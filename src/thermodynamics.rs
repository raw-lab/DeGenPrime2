//! Thermodynamics and secondary-structure checks.
//!
//! Tm uses SantaLucia (1998) unified nearest-neighbour parameters with the
//! SantaLucia (1998) entropy salt correction, ΔS += 0.368·(N−1)·ln[Na⁺], and
//! C_T/4 for non-self-complementary primers. No initiation terms are added,
//! matching the original C++ DeGenPrime so the documented 58–62 °C window keeps
//! its meaning. Gap columns are removed before stacking, and ambiguous bases use
//! the mean of the stacks they stand for.
//!
//! Units: ΔH kcal/mol, ΔS cal/(mol·K), primer nM, salt mM.

use crate::models::{iupac, nn_params, DataNode, DataSequence};
use crate::settings::Settings;
use serde::Serialize;

/// Gas constant in cal/(mol·K).
const R_CAL: f32 = 1.987;

/// Result of a Primer3-style complementarity scan between two oligos.
#[derive(Debug, Clone, Copy, Default, Serialize)]
pub struct Complementarity {
    /// Best local alignment score anywhere in the duplex.
    pub any: f32,
    /// Best alignment score anchored at a 3' end.
    pub end: f32,
    /// Alignment behind `end`: (offset k, first a-index, last a-index).
    #[serde(skip)]
    pub end_alignment: Option<(isize, usize, usize)>,
}

/// The parts of the Tm formula that depend only on the reaction:
/// R·ln(C_T/4) and ln[Na⁺]. Worked out once per run instead of once per
/// candidate window; the values are the same either way.
#[derive(Debug, Clone, Copy)]
pub struct TmConsts {
    /// `None` when the primer concentration is not positive.
    rlnk: Option<f32>,
    ln_na: f32,
}

impl TmConsts {
    pub fn new(settings: &Settings) -> Self {
        let ct = settings.primer_concentration * 1e-9;
        let rlnk = if ct <= 0.0 { None } else { Some(R_CAL * (ct / 4.0).ln()) };
        let na_m = (settings.salt_concentration / 1000.0).max(1e-6);
        TmConsts { rlnk, ln_na: na_m.ln() }
    }

    /// Tm in °C of an `n`-base oligo whose stacks sum to `dh`, `ds`.
    pub fn tm(&self, n: usize, dh: f32, ds: f32) -> f32 {
        if n < 2 {
            return 0.0;
        }
        let Some(rlnk) = self.rlnk else {
            return 0.0;
        };
        let salt = 0.368 * (n as f32 - 1.0) * self.ln_na;
        let denom = (ds + rlnk + salt) / 1000.0;
        if denom.abs() < 1e-9 {
            return 0.0;
        }
        dh / denom - 273.15
    }
}

pub struct Thermodynamics;

impl Thermodynamics {
    /// Melting temperature in °C.
    pub fn calculate_tm(seq: &DataSequence, settings: &Settings) -> f32 {
        Self::tm_of_nodes(seq.nodes(), settings)
    }

    /// Tm of the consensus bases of `nodes`, gap columns skipped (the
    /// allocation-free form of [`Self::calculate_tm`]).
    pub fn tm_of_nodes(nodes: &[DataNode], settings: &Settings) -> f32 {
        Self::tm_of_nodes_with(nodes, &TmConsts::new(settings))
    }

    /// [`Self::tm_of_nodes`] with the reaction terms worked out beforehand.
    pub fn tm_of_nodes_with(nodes: &[DataNode], k: &TmConsts) -> f32 {
        let (mut n, mut dh, mut ds) = (0usize, 0.0f32, 0.0f32);
        let mut prev: Option<char> = None;
        for node in nodes.iter().filter(|x| !x.is_gap()) {
            let c = node.most_common();
            if let Some(p) = prev {
                let (h, s) = nn_params(p, c);
                dh += h;
                ds += s;
            }
            prev = Some(c);
            n += 1;
        }
        k.tm(n, dh, ds)
    }

    /// Tm of a plain base string (no gaps), e.g. one variant of a degenerate primer.
    pub fn tm_of_bases(bases: &[char], settings: &Settings) -> f32 {
        let (mut dh, mut ds) = (0.0f32, 0.0f32);
        for w in bases.windows(2) {
            let (h, s) = nn_params(w[0], w[1]);
            dh += h;
            ds += s;
        }
        TmConsts::new(settings).tm(bases.len(), dh, ds)
    }

    /// Lowest and highest Tm across every concrete variant of a degenerate
    /// primer. Returns `None` when the primer stands for more than `limit`
    /// variants (the enumeration would be too large to be useful).
    pub fn tm_range(seq: &DataSequence, settings: &Settings, limit: u64) -> Option<(f32, f32)> {
        let masks: Vec<u8> =
            seq.nodes().iter().filter(|n| !n.is_gap()).map(|n| iupac::mask(n.code()).unwrap_or(iupac::ANY)).collect();
        let fold = masks.iter().fold(1u64, |acc, m| acc.saturating_mul(iupac::fold(*m) as u64));
        if fold > limit || masks.len() < 2 {
            return None;
        }
        let choices: Vec<Vec<char>> = masks
            .iter()
            .map(|m| iupac::BITS.iter().zip(iupac::BASES).filter(|(b, _)| m & **b != 0).map(|(_, c)| c).collect())
            .collect();
        let mut idx = vec![0usize; choices.len()];
        let mut buf: Vec<char> = choices.iter().map(|c| c[0]).collect();
        let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
        loop {
            let tm = Self::tm_of_bases(&buf, settings);
            lo = lo.min(tm);
            hi = hi.max(tm);
            // odometer increment
            let mut k = 0;
            loop {
                if k == idx.len() {
                    return Some((lo, hi));
                }
                idx[k] += 1;
                if idx[k] < choices[k].len() {
                    buf[k] = choices[k][idx[k]];
                    break;
                }
                idx[k] = 0;
                buf[k] = choices[k][0];
                k += 1;
            }
        }
    }

    /// ΔG of the stacks at `temp_k` (kcal/mol); gaps skipped.
    pub fn calculate_gibbs(seq: &DataSequence, temp_k: f32) -> f32 {
        let bases: Vec<char> = seq.nodes().iter().filter(|n| !n.is_gap()).map(|n| n.most_common()).collect();
        let (mut dh, mut ds) = (0.0f32, 0.0f32);
        for w in bases.windows(2) {
            let (h, s) = nn_params(w[0], w[1]);
            dh += h;
            ds += s;
        }
        dh - (temp_k * ds) / 1000.0
    }

    /// Stability of the last five 3' bases: −ΔG(37 °C) of their four stacks,
    /// in kcal/mol (bigger = more stable 3' end, as in Primer3's report).
    pub fn end_stability(seq: &DataSequence) -> f32 {
        let un = seq.ungapped();
        let n = un.len();
        if n < 2 {
            return 0.0;
        }
        let tail = un.sub_seq(n.saturating_sub(5), 5);
        -Self::calculate_gibbs(&tail, 310.15)
    }

    /// Legacy hairpin heuristic kept for continuity with 2.x output: the
    /// best number of base pairs (0–3) between a 3-base stem and a
    /// downstream 3-base stem separated by a loop of at least 5 bases.
    pub fn hairpin_score(seq: &DataSequence) -> usize {
        let bases: Vec<char> = seq.consensus_oligo().chars().collect();
        let n = bases.len();
        if n < 11 {
            return 0;
        }
        let mut best = 0;
        for i in 0..n - 10 {
            for j in (i + 8)..=(n - 3) {
                let matches = (0..3).filter(|&k| DataNode::complement_base(bases[i + k]) == bases[j + 2 - k]).count();
                best = best.max(matches);
            }
        }
        best
    }

    /// Legacy 3'-dimer score: longest run of base pairs when the 3' ends of
    /// the two primers overlap (up to 8 bases). 2.1 compared the forward 3'
    /// end with the reverse primer's 5' end, which cannot form an extendable
    /// dimer; this pairs the two 3' ends as intended.
    pub fn dimer_score(fwd: &DataSequence, rev: &DataSequence) -> usize {
        let f: Vec<char> = fwd.consensus_oligo().chars().collect();
        let r: Vec<char> = rev.consensus_oligo().chars().collect();
        if f.is_empty() || r.is_empty() {
            return 0;
        }
        let mut best = 0;
        for k in 1..=f.len().min(r.len()).min(8) {
            let f_tail = &f[f.len() - k..];
            let r_tail = &r[r.len() - k..];
            let matches =
                f_tail.iter().zip(r_tail.iter().rev()).filter(|(a, b)| DataNode::complement_base(**a) == **b).count();
            best = best.max(matches);
        }
        best
    }

    /// Primer3-style complementarity between `a` and `b` (both 5'→3').
    ///
    /// The two oligos are aligned antiparallel at every offset without gaps.
    /// Each position scores +1 for a Watson–Crick pair, −1 for a mismatch and
    /// −0.25 when an `N` is involved (Primer3's conventions). `any` is the best
    /// local segment; `end` is the best segment that includes the 3'-terminal
    /// base of either oligo — the dimers polymerase can extend. Gapped Primer3
    /// scores are at least as high, so these are a close lower bound.
    pub fn complementarity(a: &DataSequence, b: &DataSequence) -> Complementarity {
        let am = masks_of(a);
        let bm = masks_of(b);
        complementarity_masks(&am, &bm)
    }

    /// Self-complementarity of one oligo (self-dimer / hairpin propensity).
    pub fn self_complementarity(seq: &DataSequence) -> Complementarity {
        let m = masks_of(seq);
        complementarity_masks(&m, &m)
    }

    /// [`Self::self_complementarity`] of the consensus bases of `nodes`,
    /// reusing `buf` for the base masks.
    pub fn self_complementarity_nodes(nodes: &[DataNode], buf: &mut Vec<u8>) -> Complementarity {
        buf.clear();
        buf.extend(nodes.iter().filter(|n| !n.is_gap()).map(|n| mask_of_char(n.most_common())));
        complementarity_masks(buf, buf)
    }

    /// A text drawing of the 3'-anchored duplex behind `c.end`, e.g.
    ///
    /// ```text
    /// 5' ATGCAAGGCC 3'
    ///          ||||
    ///       3' CCGGAACGTA 5'
    /// ```
    pub fn duplex_diagram(a: &DataSequence, b: &DataSequence, c: &Complementarity) -> Option<String> {
        let (k, from, to) = c.end_alignment?;
        let a_s: Vec<char> = a.consensus_oligo().chars().collect();
        let b_s: Vec<char> = b.consensus_oligo().chars().collect();
        let b_rev: Vec<char> = b_s.iter().rev().copied().collect();
        // Column of a[0] and of b_rev[0] so that a[i] sits above b_rev[i - k].
        let (a_col, b_col) = if k >= 0 { (0usize, k as usize) } else { ((-k) as usize, 0usize) };
        let pad = 3; // width of "5' "
        let mut top = String::from("5' ");
        top.push_str(&" ".repeat(a_col));
        top.extend(a_s.iter());
        top.push_str(" 3'");
        let mut mid = " ".repeat(pad + a_col);
        for (i, &c) in a_s.iter().enumerate() {
            let t = i as isize - k;
            let paired = i >= from
                && i <= to
                && t >= 0
                && (t as usize) < b_rev.len()
                && pairs(mask_of_char(c), mask_of_char(b_rev[t as usize]));
            mid.push(if paired { '|' } else { ' ' });
        }
        let mid = mid.trim_end().to_string();
        let mut bot = String::from("3' ");
        bot.push_str(&" ".repeat(b_col));
        bot.extend(b_rev.iter());
        bot.push_str(" 5'");
        Some(format!("{top}\n{mid}\n{bot}"))
    }
}

fn mask_of_char(c: char) -> u8 {
    iupac::mask(c).unwrap_or(iupac::ANY)
}

fn masks_of(seq: &DataSequence) -> Vec<u8> {
    seq.nodes().iter().filter(|n| !n.is_gap()).map(|n| mask_of_char(n.most_common())).collect()
}

const fn pairs(x: u8, y: u8) -> bool {
    x & iupac::complement(y) != 0
}

fn pair_score(x: u8, y: u8) -> f32 {
    PAIR_SCORES[(x & 15) as usize][(y & 15) as usize]
}

/// [`pair_score_slow`] for every pair of masks.
static PAIR_SCORES: [[f32; 16]; 16] = {
    let mut t = [[0.0f32; 16]; 16];
    let mut x = 0;
    while x < 16 {
        let mut y = 0;
        while y < 16 {
            t[x][y] = pair_score_slow(x as u8, y as u8);
            y += 1;
        }
        x += 1;
    }
    t
};

const fn pair_score_slow(x: u8, y: u8) -> f32 {
    if x == iupac::ANY || y == iupac::ANY {
        -0.25
    } else if pairs(x, y) {
        1.0
    } else {
        -1.0
    }
}

fn complementarity_masks(a: &[u8], b: &[u8]) -> Complementarity {
    let (la, lb) = (a.len() as isize, b.len() as isize);
    let mut out = Complementarity::default();
    if la == 0 || lb == 0 {
        return out;
    }
    // b reversed: br[t] = b[lb-1-t]; a[i] pairs with br[i-k].
    // b reversed without copying: br[t] = b[lb - 1 - t].
    let score_at = |i: isize, k: isize| pair_score(a[i as usize], b[(lb - 1 - (i - k)) as usize]);
    for k in -(lb - 1)..=(la - 1) {
        let i0 = k.max(0);
        let i1 = (la - 1).min(k + lb - 1);
        // Best local segment (Kadane).
        let mut cur = 0.0f32;
        for i in i0..=i1 {
            let s = score_at(i, k);
            cur = if cur > 0.0 { cur + s } else { s };
            out.any = out.any.max(cur);
        }
        // Anchored at a's 3' end (a[la-1]).
        if i1 == la - 1 {
            let mut sum = 0.0f32;
            for i in (i0..=i1).rev() {
                sum += score_at(i, k);
                if sum > out.end {
                    out.end = sum;
                    out.end_alignment = Some((k, i as usize, i1 as usize));
                }
            }
        }
        // Anchored at b's 3' end (br[0], i.e. i == k).
        if i0 == k {
            let mut sum = 0.0f32;
            for i in i0..=i1 {
                sum += score_at(i, k);
                if sum > out.end {
                    out.end = sum;
                    out.end_alignment = Some((k, i0 as usize, i as usize));
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::DataSequence;

    #[test]
    fn tm_reasonable_for_typical_primer() {
        let s = DataSequence::from_string("ATGCAATTGGCCGGAATTCC");
        let tm = Thermodynamics::calculate_tm(&s, &Settings::default());
        assert!(tm > 40.0 && tm < 80.0, "tm = {tm}");
    }

    #[test]
    fn tm_unchanged_from_2_1_for_plain_primers() {
        // Regression anchor: 2.1 reported 60.93 °C for this anfD primer.
        let s = DataSequence::from_string("GCTACCTGAACACCATCCCG");
        let tm = Thermodynamics::calculate_tm(&s, &Settings::default());
        assert!((tm - 60.93).abs() < 0.02, "tm = {tm}");
    }

    #[test]
    fn gc_rich_has_higher_tm_than_at_rich() {
        let settings = Settings::default();
        let gc = Thermodynamics::calculate_tm(&DataSequence::from_string("GGCGGCGGCGGCGGCGGCGG"), &settings);
        let at = Thermodynamics::calculate_tm(&DataSequence::from_string("ATATATATATATATATATAT"), &settings);
        assert!(gc > at);
    }

    #[test]
    fn salt_correction_changes_tm() {
        let s = DataSequence::from_string("ATGCAATTGGCCGGAATTCC");
        let low = Thermodynamics::calculate_tm(&s, &Settings::default().with_salt_concentration(10.0));
        let high = Thermodynamics::calculate_tm(&s, &Settings::default().with_salt_concentration(200.0));
        assert!(high > low);
    }

    #[test]
    fn degenerate_tm_lies_inside_variant_range() {
        let settings = Settings::default();
        let s = DataSequence::from_string("ATGCRYCCGGAATTCCAGTC");
        let tm = Thermodynamics::calculate_tm(&s, &settings);
        let (lo, hi) = Thermodynamics::tm_range(&s, &settings, 4096).unwrap();
        assert!(lo <= tm + 0.5 && tm <= hi + 0.5, "{lo} <= {tm} <= {hi}");
        assert!(tm > 55.0, "IUPAC bases must not collapse to N (2.1 gave 52.8): {tm}");
    }

    #[test]
    fn dimer_score_pairs_three_prime_ends() {
        let f = DataSequence::from_string("ATATATATATATATGGGCCC");
        let r3 = DataSequence::from_string("TATATATATATATAGGGCCC");
        let r5 = DataSequence::from_string("GGGCCCTATATATATATATA");
        assert!(Thermodynamics::dimer_score(&f, &r3) >= 6);
        assert!(Thermodynamics::dimer_score(&f, &r3) > Thermodynamics::dimer_score(&f, &r5));
    }

    #[test]
    fn complementarity_scores_palindromes_high() {
        // GAATTC is its own reverse complement: a 3'-terminal palindrome.
        let pal = DataSequence::from_string("ACGTAGCTAGGAATTC");
        let c = Thermodynamics::self_complementarity(&pal);
        assert!(c.end >= 6.0, "end = {}", c.end);
        assert!(c.any >= c.end);
        let plain = DataSequence::from_string("AAAACCCCAAAACCCC");
        assert!(Thermodynamics::self_complementarity(&plain).end < 2.0);
    }

    #[test]
    fn duplex_diagram_draws_pairs() {
        let pal = DataSequence::from_string("ACGTAGCTAGGAATTC");
        let c = Thermodynamics::self_complementarity(&pal);
        let d = Thermodynamics::duplex_diagram(&pal, &pal, &c).unwrap();
        assert!(d.contains("||||||"), "{d}");
        assert_eq!(d.lines().count(), 3);
    }
}
