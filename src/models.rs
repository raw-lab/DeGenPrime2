//! Core data models for DeGenPrime.
//!
//! * [`DataNode`] — one alignment column (or one primer base): its IUPAC code,
//!   its most common base, and the fraction of sequences that agree with it.
//! * [`DataSequence`] — a run of nodes: a consensus built from an alignment, a
//!   candidate primer window, or a primer typed in by the user.
//! * [`Primer`] / [`PrimerPair`] — candidate windows and scored pairs.
//!
//! Nucleotides are handled internally as 4-bit masks (A=1, C=2, G=4, T=8), so
//! IUPAC unions, complements and degeneracy folds are bit operations and every
//! result is deterministic (2.1 broke ties by iterating a `HashMap`, whose order
//! is randomised per process, so the same input could give different primers).
//!
//! Nearest-neighbour parameters follow SantaLucia (1998) PNAS 95(4):1460.

use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

/// IUPAC nucleotide codes as 4-bit masks.
pub mod iupac {
    pub const A: u8 = 0b0001;
    pub const C: u8 = 0b0010;
    pub const G: u8 = 0b0100;
    pub const T: u8 = 0b1000;
    pub const ANY: u8 = 0b1111;

    /// The four concrete bases in mask order.
    pub const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    pub const BITS: [u8; 4] = [A, C, G, T];

    /// Mask for an IUPAC character (case-insensitive; `U` reads as `T`,
    /// `X` as `N`). Returns `None` for gaps and non-nucleotide characters.
    pub fn mask(c: char) -> Option<u8> {
        if c.is_ascii() {
            let m = MASKS[c as usize];
            return if m == 0 { None } else { Some(m) };
        }
        None
    }

    /// `mask` for every ASCII character (0 = not a nucleotide code).
    static MASKS: [u8; 128] = {
        let mut t = [0u8; 128];
        let mut i = 0;
        while i < 128 {
            t[i] = match mask_slow(i as u8 as char) {
                Some(m) => m,
                None => 0,
            };
            i += 1;
        }
        t
    };

    const fn mask_slow(c: char) -> Option<u8> {
        Some(match c.to_ascii_uppercase() {
            'A' => A,
            'C' => C,
            'G' => G,
            'T' | 'U' => T,
            'R' => A | G,
            'Y' => C | T,
            'S' => C | G,
            'W' => A | T,
            'K' => G | T,
            'M' => A | C,
            'B' => C | G | T,
            'D' => A | G | T,
            'H' => A | C | T,
            'V' => A | C | G,
            'N' | 'X' => ANY,
            _ => return None,
        })
    }

    /// IUPAC letter for a mask. An empty mask has no letter and maps to `N`.
    pub fn letter(mask: u8) -> char {
        match mask & ANY {
            0b0001 => 'A',
            0b0010 => 'C',
            0b0011 => 'M',
            0b0100 => 'G',
            0b0101 => 'R',
            0b0110 => 'S',
            0b0111 => 'V',
            0b1000 => 'T',
            0b1001 => 'W',
            0b1010 => 'Y',
            0b1011 => 'H',
            0b1100 => 'K',
            0b1101 => 'D',
            0b1110 => 'B',
            _ => 'N',
        }
    }

    /// Watson–Crick complement of a mask (A<->T, C<->G).
    pub const fn complement(m: u8) -> u8 {
        ((m & A) << 3) | ((m & T) >> 3) | ((m & C) << 1) | ((m & G) >> 1)
    }

    /// Number of concrete bases a code stands for (1 for A/C/G/T, 4 for N).
    pub fn fold(m: u8) -> u32 {
        (m & ANY).count_ones()
    }

    /// Expected G+C fraction of a (possibly ambiguous) base.
    pub fn gc_fraction(m: u8) -> f32 {
        GC_FRACTION[(m & ANY) as usize]
    }

    /// G+C fraction of every mask (0 for the empty mask).
    static GC_FRACTION: [f32; 16] = {
        let mut t = [0.0f32; 16];
        let mut m = 1;
        while m < 16 {
            t[m] = (m as u8 & (C | G)).count_ones() as f32 / (m as u8).count_ones() as f32;
            m += 1;
        }
        t
    };

    /// True for characters accepted in sequences: IUPAC letters plus gaps.
    pub fn is_sequence_char(c: char) -> bool {
        matches!(c, '-' | '.') || mask(c).is_some()
    }
}

/// SantaLucia (1998) unified ΔH (kcal/mol), indexed `[5' base][3' base]` in A,C,G,T order.
const NN_DH: [[f32; 4]; 4] =
    [[-7.9, -8.4, -7.8, -7.2], [-8.5, -8.0, -10.6, -7.8], [-8.2, -9.8, -8.0, -8.4], [-7.2, -8.2, -8.5, -7.9]];

/// SantaLucia (1998) unified ΔS (cal/(mol·K)), same indexing as [`NN_DH`].
const NN_DS: [[f32; 4]; 4] = [
    [-22.2, -22.4, -21.0, -20.4],
    [-22.7, -19.9, -27.2, -21.0],
    [-22.2, -24.4, -19.9, -22.4],
    [-21.3, -22.2, -22.7, -22.2],
];

/// Nearest-neighbour (ΔH, ΔS) for the stack `b1 -> b2`. Ambiguous codes are
/// averaged over every base they stand for (so `R` contributes the mean of
/// the A and G stacks). Gaps contribute nothing.
pub fn nn_params(b1: char, b2: char) -> (f32, f32) {
    match (iupac::mask(b1), iupac::mask(b2)) {
        (Some(x), Some(y)) => nn_table()[x as usize][y as usize],
        _ => (0.0, 0.0),
    }
}

/// [`nn_params`] for every pair of IUPAC masks, computed once. Tm is
/// evaluated for millions of windows on long alignments; the table gives
/// exactly the numbers the averaging loop gives, without repeating it.
fn nn_table() -> &'static [[(f32, f32); 16]; 16] {
    static TABLE: std::sync::OnceLock<[[(f32, f32); 16]; 16]> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        let mut t = [[(0.0f32, 0.0f32); 16]; 16];
        for (m1, row) in t.iter_mut().enumerate() {
            for (m2, cell) in row.iter_mut().enumerate() {
                *cell = nn_params_masks(m1 as u8, m2 as u8);
            }
        }
        t
    })
}

fn nn_params_masks(m1: u8, m2: u8) -> (f32, f32) {
    let (mut dh, mut ds, mut n) = (0.0f32, 0.0f32, 0.0f32);
    for (i, bi) in iupac::BITS.iter().enumerate() {
        if m1 & bi == 0 {
            continue;
        }
        for (j, bj) in iupac::BITS.iter().enumerate() {
            if m2 & bj == 0 {
                continue;
            }
            dh += NN_DH[i][j];
            ds += NN_DS[i][j];
            n += 1.0;
        }
    }
    if n == 0.0 {
        (0.0, 0.0)
    } else {
        (dh / n, ds / n)
    }
}

/// A single aligned column (or a single primer base).
///
/// * `code` — IUPAC code covering every base observed in the column, or `-`
///   when most sequences have a gap here.
/// * `most_common` — the consensus base (ties resolve in A, C, G, T order).
/// * `ratio` — fraction of sequences that carry the consensus state. Gaps
///   count against conservation: a column that is `A` in one sequence and a
///   gap in four is 20 % conserved, not 100 % (the 2.1 behaviour).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DataNode {
    code: char,
    most_common: char,
    ratio: f32,
}

impl DataNode {
    /// Build a node from the characters observed in one alignment column.
    pub fn new(chars: Vec<char>) -> Self {
        Self::from_column(chars)
    }

    /// Build a node from an iterator over one alignment column.
    ///
    /// Ambiguous calls in the input (e.g. `R`) add fractional counts to each
    /// base they stand for and widen the IUPAC code. `N` carries no
    /// information and is ignored. A column becomes a gap node (`-`) when
    /// more than half of the informative sequences have a gap.
    pub fn from_column<I: IntoIterator<Item = char>>(column: I) -> Self {
        let mut acc = ColumnCounts::default();
        for c in column {
            acc.add(if c.is_ascii() { c as u8 } else { 0 });
        }
        acc.node()
    }
}

/// Running tally of one alignment column, filled one sequence at a time.
/// [`crate::io::SequenceReader::alignment_to_consensus`] keeps one per
/// column and reads the alignment row by row, which is cache-friendly and
/// never copies the sequences; [`DataNode::from_column`] uses the same
/// arithmetic, in the same order, so both give identical nodes.
#[derive(Debug, Clone, Copy, Default)]
pub struct ColumnCounts {
    counts: [f32; 4],
    union: u8,
    gaps: u32,
    called: u32,
}

impl ColumnCounts {
    /// Add one sequence's character (ASCII; case-insensitive).
    #[inline]
    pub fn add(&mut self, c: u8) {
        if c == b'-' || c == b'.' {
            self.gaps += 1;
            return;
        }
        match iupac::mask(c as char) {
            Some(iupac::ANY) | None => {}
            Some(m) => {
                let share = 1.0 / iupac::fold(m) as f32;
                for (i, bit) in iupac::BITS.iter().enumerate() {
                    if m & bit != 0 {
                        self.counts[i] += share;
                    }
                }
                self.union |= m;
                self.called += 1;
            }
        }
    }

    /// The node for this column.
    pub fn node(&self) -> DataNode {
        let (counts, union, gaps, called) = (self.counts, self.union, self.gaps, self.called);
        let informative = called + gaps;
        if informative == 0 {
            return DataNode { code: 'N', most_common: 'N', ratio: 0.0 };
        }
        if gaps * 2 > informative {
            return DataNode { code: '-', most_common: '-', ratio: gaps as f32 / informative as f32 };
        }

        let mut best = 0usize;
        for i in 1..4 {
            if counts[i] > counts[best] {
                best = i;
            }
        }
        DataNode {
            code: iupac::letter(union),
            most_common: iupac::BASES[best],
            ratio: counts[best] / informative as f32,
        }
    }
}

impl DataNode {
    /// A node for one base typed by the user. IUPAC codes are preserved
    /// (2.1 turned every ambiguous base into `N`). `ratio` is the chance a
    /// random variant of the degenerate base matches: 1 for A/C/G/T, ½ for R.
    pub fn from_base(base: char) -> Self {
        if base == '-' || base == '.' {
            return DataNode { code: '-', most_common: '-', ratio: 1.0 };
        }
        match iupac::mask(base) {
            Some(m) => {
                let letter = iupac::letter(m);
                DataNode { code: letter, most_common: letter, ratio: 1.0 / iupac::fold(m) as f32 }
            }
            None => DataNode { code: 'N', most_common: 'N', ratio: 0.0 },
        }
    }

    pub fn code(&self) -> char {
        self.code
    }
    pub fn most_common(&self) -> char {
        self.most_common
    }
    pub fn ratio(&self) -> f32 {
        self.ratio
    }
    pub fn is_gap(&self) -> bool {
        self.code == '-'
    }

    /// IUPAC complement of a single character (gaps stay gaps).
    pub fn complement_base(base: char) -> char {
        if base == '-' {
            return '-';
        }
        match iupac::mask(base) {
            Some(m) => iupac::letter(iupac::complement(m)),
            None => base,
        }
    }

    /// The complementary node (same conservation ratio).
    pub fn complement(&self) -> Self {
        DataNode {
            code: Self::complement_base(self.code),
            most_common: Self::complement_base(self.most_common),
            ratio: self.ratio,
        }
    }

    /// ΔH contribution of the stack `self -> other` (kcal/mol).
    pub fn enthalpy_to(&self, other: &DataNode) -> f32 {
        nn_params(self.most_common, other.most_common).0
    }

    /// ΔS contribution of the stack `self -> other` (cal/(mol·K)).
    pub fn entropy_to(&self, other: &DataNode) -> f32 {
        nn_params(self.most_common, other.most_common).1
    }
}

/// A sequence of [`DataNode`]s: a consensus, a primer window, or a primer.
#[derive(Debug, Clone)]
pub struct DataSequence {
    nodes: Vec<DataNode>,
    start_pos: usize,
}

impl DataSequence {
    pub fn new(nodes: Vec<DataNode>) -> Self {
        DataSequence { nodes, start_pos: 0 }
    }

    pub fn with_start(nodes: Vec<DataNode>, start_pos: usize) -> Self {
        DataSequence { nodes, start_pos }
    }

    /// Build from a plain primer string; IUPAC codes are kept as typed.
    pub fn from_string(seq: &str) -> Self {
        let nodes = seq.chars().filter(|c| !c.is_whitespace()).map(DataNode::from_base).collect();
        DataSequence::new(nodes)
    }

    pub fn len(&self) -> usize {
        self.nodes.len()
    }
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
    pub fn nodes(&self) -> &[DataNode] {
        &self.nodes
    }
    pub fn start_pos(&self) -> usize {
        self.start_pos
    }

    /// Slice `length` nodes starting at `start`. Out-of-range requests are
    /// clipped instead of panicking.
    pub fn sub_seq(&self, start: usize, length: usize) -> DataSequence {
        let start = start.min(self.nodes.len());
        let end = start.saturating_add(length).min(self.nodes.len());
        DataSequence::with_start(self.nodes[start..end].to_vec(), self.start_pos + start)
    }

    /// The same sequence with gap columns removed — what you actually order.
    pub fn ungapped(&self) -> DataSequence {
        DataSequence::with_start(self.nodes.iter().copied().filter(|n| !n.is_gap()).collect(), self.start_pos)
    }

    /// Number of non-gap positions.
    pub fn actual_size(&self) -> usize {
        self.nodes.iter().filter(|n| !n.is_gap()).count()
    }

    pub fn reverse(&self) -> DataSequence {
        let mut rev = self.nodes.clone();
        rev.reverse();
        DataSequence::new(rev)
    }

    pub fn complement(&self) -> DataSequence {
        DataSequence::new(self.nodes.iter().map(|n| n.complement()).collect())
    }

    pub fn reverse_complement(&self) -> DataSequence {
        self.complement().reverse()
    }

    /// IUPAC codes, gaps included (alignment coordinates).
    pub fn codes(&self) -> String {
        self.nodes.iter().map(|n| n.code()).collect()
    }

    /// The degenerate oligo to order: IUPAC codes with gap columns removed.
    pub fn oligo(&self) -> String {
        self.nodes.iter().filter(|n| !n.is_gap()).map(|n| n.code()).collect()
    }

    /// Most-common base at each position (no degeneracy), gaps included.
    pub fn consensus_string(&self) -> String {
        self.nodes.iter().map(|n| n.most_common()).collect()
    }

    /// Most-common bases with gap columns removed.
    pub fn consensus_oligo(&self) -> String {
        self.nodes.iter().filter(|n| !n.is_gap()).map(|n| n.most_common()).collect()
    }

    /// G+C fraction over non-gap positions, using the consensus base.
    /// Ambiguous bases count fractionally (S = 1, R = ½, N = ½, W = 0).
    pub fn gc_ratio(&self) -> f32 {
        measure::gc_ratio(&self.nodes)
    }

    /// Conservation penalty: Σ (1 − ratio) × 100. Lower is better.
    pub fn penalty(&self) -> f32 {
        measure::penalty(&self.nodes)
    }

    /// Mean conservation of the window (0–1).
    pub fn mean_conservation(&self) -> f32 {
        if self.nodes.is_empty() {
            return 0.0;
        }
        self.nodes.iter().map(|n| n.ratio()).sum::<f32>() / self.nodes.len() as f32
    }

    /// Number of degenerate (non-A/C/G/T) positions, gaps excluded.
    pub fn degeneracy_count(&self) -> usize {
        self.nodes.iter().filter(|n| !matches!(n.code(), 'A' | 'T' | 'G' | 'C' | '-')).count()
    }

    /// Degeneracy fold: how many distinct sequences the degenerate oligo
    /// stands for (product of per-position folds, saturating).
    pub fn fold(&self) -> u64 {
        self.nodes
            .iter()
            .filter(|n| !n.is_gap())
            .map(|n| iupac::mask(n.code()).map(iupac::fold).unwrap_or(1) as u64)
            .fold(1u64, |acc, f| acc.saturating_mul(f))
    }

    /// Number of gap columns.
    pub fn gap_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_gap()).count()
    }

    /// Quick hairpin heuristic: are the first three and last three bases of
    /// the oligo complementary (5'-ABC…C'B'A'-3')?
    pub fn has_complementary_ends(&self) -> bool {
        measure::complementary_ends(&self.nodes)
    }
}

/// The [`DataSequence`] measures on plain node slices. The design pipeline
/// screens hundreds of thousands of windows with these without copying
/// them into a `DataSequence` first; the arithmetic, and its order, is the
/// same, so both paths give bit-identical numbers.
pub mod measure {
    use super::{iupac, DataNode};

    /// G+C fraction over non-gap positions, using the consensus base.
    pub fn gc_ratio(nodes: &[DataNode]) -> f32 {
        let mut total = 0usize;
        let mut gc = 0.0f32;
        for n in nodes.iter().filter(|n| !n.is_gap()) {
            total += 1;
            if let Some(m) = iupac::mask(n.most_common()) {
                gc += iupac::gc_fraction(m);
            }
        }
        if total == 0 {
            0.0
        } else {
            gc / total as f32
        }
    }

    /// Conservation penalty: Σ (1 − ratio) × 100, gap columns included.
    pub fn penalty(nodes: &[DataNode]) -> f32 {
        nodes.iter().map(|n| (1.0 - n.ratio()) * 100.0).sum::<f32>()
    }

    /// Are the first three and last three non-gap consensus bases
    /// complementary (5'-ABC…C'B'A'-3')?
    pub fn complementary_ends(nodes: &[DataNode]) -> bool {
        if nodes.iter().filter(|n| !n.is_gap()).count() < 6 {
            return false;
        }
        let mut front = nodes.iter().filter(|n| !n.is_gap());
        let mut back = nodes.iter().rev().filter(|n| !n.is_gap());
        (0..3).all(|_| match (front.next(), back.next()) {
            (Some(a), Some(b)) => DataNode::complement_base(a.most_common()) == b.most_common(),
            _ => false,
        })
    }
}

/// A candidate primer: a window `[index, index + length)` of a parent sequence.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Primer {
    pub index: usize,
    pub length: usize,
    pub penalty: f32,
}

impl Primer {
    pub fn new(index: usize, length: usize) -> Self {
        Primer { index, length, penalty: 0.0 }
    }

    pub fn with_penalty(mut self, penalty: f32) -> Self {
        self.penalty = penalty;
        self
    }

    /// Last position covered (inclusive). Zero-length primers return `index`.
    pub fn end_index(&self) -> usize {
        (self.index + self.length).saturating_sub(1).max(self.index)
    }
}

impl Eq for Primer {}

impl Ord for Primer {
    /// Penalty first, then position and length — consistent with `Eq`
    /// (2.1 compared penalty only, which breaks `Ord`'s contract).
    fn cmp(&self, other: &Self) -> Ordering {
        self.penalty
            .total_cmp(&other.penalty)
            .then_with(|| self.index.cmp(&other.index))
            .then_with(|| self.length.cmp(&other.length))
    }
}

impl PartialOrd for Primer {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A forward + reverse primer pair.
///
/// `reverse.index` is an offset on the **reverse-complement** strand; use
/// [`PrimerPair::reverse_fwd_strand_position`] for forward-strand coordinates.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PrimerPair {
    pub forward: Primer,
    pub reverse: Primer,
    pub amplicon_size: usize,
    pub tm_diff: f32,
    pub fwd_tm: f32,
    pub rev_tm: f32,
    /// Pair objective used for ranking (lower is better); see `calculator`.
    #[serde(default)]
    pub score: f32,
}

impl PrimerPair {
    pub fn new(forward: Primer, reverse: Primer, amplicon_size: usize, tm_diff: f32, fwd_tm: f32, rev_tm: f32) -> Self {
        let score = forward.penalty + reverse.penalty + tm_diff;
        PrimerPair { forward, reverse, amplicon_size, tm_diff, fwd_tm, rev_tm, score }
    }

    pub fn with_score(mut self, score: f32) -> Self {
        self.score = score;
        self
    }

    /// Left-most forward-strand position (0-based) of the reverse primer's
    /// footprint, i.e. where its 3' end sits on the template.
    pub fn reverse_fwd_strand_position(&self, template_len: usize) -> usize {
        template_len.saturating_sub(self.reverse.index).saturating_sub(self.reverse.length)
    }

    fn key_cmp(&self, other: &Self) -> Ordering {
        self.forward
            .index
            .cmp(&other.forward.index)
            .then_with(|| self.forward.length.cmp(&other.forward.length))
            .then_with(|| self.reverse.index.cmp(&other.reverse.index))
            .then_with(|| self.reverse.length.cmp(&other.reverse.length))
    }
}

impl PartialEq for PrimerPair {
    fn eq(&self, other: &Self) -> bool {
        self.forward == other.forward && self.reverse == other.reverse
    }
}
impl Eq for PrimerPair {}

impl Ord for PrimerPair {
    /// Best score first; ties broken by position so ordering is total.
    fn cmp(&self, other: &Self) -> Ordering {
        self.score.total_cmp(&other.score).then_with(|| self.key_cmp(other))
    }
}

impl PartialOrd for PrimerPair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A raw biological sequence read from a file (FASTA / CLUSTAL).
#[derive(Debug, Clone)]
pub struct Sequence {
    pub id: String,
    pub sequence: String,
}

impl Sequence {
    pub fn new(id: String, sequence: String) -> Self {
        Sequence { id, sequence }
    }

    pub fn len(&self) -> usize {
        self.sequence.len()
    }
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Uppercase, strip whitespace.
    pub fn normalize(&self) -> String {
        self.sequence.chars().filter(|c| !c.is_whitespace()).map(|c| c.to_ascii_uppercase()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_base_works_for_ambig() {
        assert_eq!(DataNode::complement_base('A'), 'T');
        assert_eq!(DataNode::complement_base('R'), 'Y');
        assert_eq!(DataNode::complement_base('N'), 'N');
        assert_eq!(DataNode::complement_base('-'), '-');
        for c in "ACGTRYSWKMBDHVN".chars() {
            let back = DataNode::complement_base(DataNode::complement_base(c));
            assert_eq!(back, c, "complement must be an involution for {c}");
        }
    }

    #[test]
    fn datanode_picks_iupac_pair() {
        let node = DataNode::new(vec!['A', 'A', 'A', 'G']);
        assert_eq!(node.code(), 'R');
        assert_eq!(node.most_common(), 'A');
        assert!((node.ratio() - 0.75).abs() < 1e-6);
    }

    #[test]
    fn ties_resolve_deterministically() {
        // 2.1 used HashMap iteration order here, so this flipped between runs.
        for _ in 0..50 {
            let node = DataNode::new(vec!['G', 'A', 'G', 'A']);
            assert_eq!(node.most_common(), 'A');
            assert_eq!(node.code(), 'R');
        }
    }

    #[test]
    fn gaps_count_against_conservation() {
        let mostly_gap = DataNode::new(vec!['A', '-', '-', '-', '-']);
        assert!(mostly_gap.is_gap());
        assert_eq!(DataSequence::new(vec![mostly_gap]).gap_count(), 1);
        let some_gap = DataNode::new(vec!['A', 'A', 'A', '-', '-']);
        assert_eq!(some_gap.code(), 'A');
        assert!((some_gap.ratio() - 0.6).abs() < 1e-6);
    }

    #[test]
    fn ambiguous_input_calls_widen_the_code() {
        let node = DataNode::new(vec!['A', 'A', 'R']);
        assert_eq!(node.code(), 'R');
        assert_eq!(node.most_common(), 'A');
    }

    #[test]
    fn user_primer_keeps_iupac_codes() {
        let s = DataSequence::from_string("ATGCRYN");
        assert_eq!(s.codes(), "ATGCRYN");
        assert_eq!(s.degeneracy_count(), 3);
        assert_eq!(s.fold(), 2 * 2 * 4);
    }

    #[test]
    fn reverse_complement_is_involutive() {
        let s = DataSequence::from_string("ATGCAATTGCRY");
        let back = s.reverse_complement().reverse_complement();
        assert_eq!(s.codes(), back.codes());
    }

    #[test]
    fn gc_ratio_correct() {
        assert!((DataSequence::from_string("GGCC").gc_ratio() - 1.0).abs() < 1e-6);
        assert!((DataSequence::from_string("ATAT").gc_ratio()).abs() < 1e-6);
        assert!((DataSequence::from_string("ATGC").gc_ratio() - 0.5).abs() < 1e-6);
        assert!((DataSequence::from_string("SSWW").gc_ratio() - 0.5).abs() < 1e-6);
    }

    #[test]
    fn sub_seq_never_panics() {
        let s = DataSequence::from_string("ACGT");
        assert_eq!(s.sub_seq(10, 5).len(), 0);
        assert_eq!(s.sub_seq(2, 99).codes(), "GT");
    }

    #[test]
    fn ambiguous_stacks_are_averaged() {
        let (dh_r, _) = nn_params('R', 'C');
        let (dh_a, _) = nn_params('A', 'C');
        let (dh_g, _) = nn_params('G', 'C');
        assert!((dh_r - (dh_a + dh_g) / 2.0).abs() < 1e-5);
        assert_eq!(nn_params('-', 'A'), (0.0, 0.0));
    }

    #[test]
    fn primer_ordering_is_total() {
        let a = Primer::new(1, 20).with_penalty(5.0);
        let b = Primer::new(2, 20).with_penalty(5.0);
        assert_ne!(a.cmp(&b), Ordering::Equal);
    }
}
