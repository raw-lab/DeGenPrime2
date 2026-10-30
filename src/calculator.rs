//! Primer design pipeline.
//!
//! 1. Every window of `min..=max` primer length inside the design region, on
//!    both strands, is a candidate. The region is mirrored onto the reverse
//!    complement. Windows are screened in parallel, in chunks of the
//!    alignment, straight from the consensus (no per-window copies).
//! 2. Each window goes through [`filters::screen_nodes`]. Windows that start
//!    or end on a gap column are skipped: they duplicate a shorter primer.
//! 3. Forward and reverse survivors are paired and ranked by
//!
//!    `score = Pf + Pr + |Tm_f − Tm_opt| + |Tm_r − Tm_opt| + |ΔTm|`
//!
//!    where `P` is the conservation penalty Σ(1 − ratio)·100 and `Tm_opt` is the
//!    middle of the Tm window. One column at 80 % conservation costs 20 points,
//!    so conservation dominates and the Tm terms break ties. `RankBy::TmDiff`
//!    restores the 2.1 ordering (smallest ΔTm first).
//! 4. The pair search never looks at every forward × reverse combination.
//!    Reverse candidates are indexed by where they end ([`BlockTree`]), so a
//!    forward primer only meets partners that give a product of the allowed
//!    size; inside that window, partners come in order of cost (or of Tm, for
//!    `TmDiff`), and the search stops as soon as no remaining partner can beat
//!    the pairs already held. Only the best few thousand pairs are held, in a
//!    bounded heap. The final list skips near-duplicate pairs and pairs whose
//!    3' ends can form a primer-dimer; if too many of the held pairs were
//!    skipped, the search runs again with a larger pool.
//!
//! The result is exactly the best pairs under the ranking (ties broken by
//! position), whatever order the search visits them in, and it does not
//! depend on the number of threads.

use crate::filters::{self, FilterResult, Rejection, Scratch};
use crate::models::{measure, DataSequence, Primer, PrimerPair};
use crate::settings::{RankBy, Settings};
use crate::thermodynamics::Thermodynamics;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Forward,
    Reverse,
}

/// Why candidate primers on one strand were rejected (first failure counts).
#[derive(Debug, Clone, Default, Serialize)]
pub struct StrandStats {
    pub considered: usize,
    pub degeneracy: usize,
    pub deletions: usize,
    pub gc_content: usize,
    pub gc_3prime: usize,
    pub repeats: usize,
    pub complementary_ends: usize,
    pub tm_low: usize,
    pub tm_high: usize,
    pub self_complementarity: usize,
    pub ok: usize,
}

impl StrandStats {
    fn record(&mut self, r: Rejection) {
        match r {
            Rejection::Degeneracy => self.degeneracy += 1,
            Rejection::Deletions => self.deletions += 1,
            Rejection::GcContent => self.gc_content += 1,
            Rejection::Gc3Prime => self.gc_3prime += 1,
            Rejection::Repeats => self.repeats += 1,
            Rejection::ComplementaryEnds => self.complementary_ends += 1,
            Rejection::TmLow => self.tm_low += 1,
            Rejection::TmHigh => self.tm_high += 1,
            Rejection::SelfComplementarity => self.self_complementarity += 1,
        }
    }

    fn add(&mut self, o: &StrandStats) {
        self.considered += o.considered;
        self.degeneracy += o.degeneracy;
        self.deletions += o.deletions;
        self.gc_content += o.gc_content;
        self.gc_3prime += o.gc_3prime;
        self.repeats += o.repeats;
        self.complementary_ends += o.complementary_ends;
        self.tm_low += o.tm_low;
        self.tm_high += o.tm_high;
        self.self_complementarity += o.self_complementarity;
        self.ok += o.ok;
    }

    /// Rejection reasons with counts, most common first.
    pub fn top_reasons(&self) -> Vec<(&'static str, usize)> {
        let mut v = vec![
            ("too degenerate", self.degeneracy),
            ("too many gaps", self.deletions),
            ("GC content out of range", self.gc_content),
            ("3' end too GC-rich", self.gc_3prime),
            ("repeats", self.repeats),
            ("complementary ends", self.complementary_ends),
            ("Tm too low", self.tm_low),
            ("Tm too high", self.tm_high),
            ("self-complementary", self.self_complementarity),
        ];
        v.retain(|(_, n)| *n > 0);
        v.sort_by_key(|&(_, n)| std::cmp::Reverse(n));
        v
    }
}

/// What happened to candidate pairs.
///
/// `combinations`, `overlap`, `product_too_short` and `product_too_long`
/// are exact counts over every forward × reverse combination (each
/// combination is counted under the first of these it fails). The other
/// counts cover the pairs the search examined in detail, which stops early
/// once no remaining pair can beat the ones already held.
#[derive(Debug, Clone, Default, Serialize)]
pub struct PairStats {
    pub combinations: u64,
    pub examined: u64,
    pub overlap: u64,
    pub product_too_short: u64,
    pub product_too_long: u64,
    pub tm_diff: u64,
    pub pair_complementarity: u64,
    pub too_similar: u64,
    pub ok: u64,
}

impl PairStats {
    /// Combinations whose primers do not overlap and give a product of the
    /// allowed size.
    pub fn position_valid(&self) -> u64 {
        self.combinations
            .saturating_sub(self.overlap)
            .saturating_sub(self.product_too_short)
            .saturating_sub(self.product_too_long)
    }
}

/// Primer3-style statistics explaining a design run.
#[derive(Debug, Clone, Default, Serialize)]
pub struct Explain {
    pub forward: StrandStats,
    pub reverse: StrandStats,
    pub pairs: PairStats,
}

#[derive(Debug, Clone)]
pub struct DesignOutcome {
    pub pairs: Vec<PrimerPair>,
    pub explain: Explain,
}

/// A screened candidate with everything pairing needs.
#[derive(Debug, Clone, Copy)]
struct Candidate {
    primer: Primer,
    tm: f32,
    /// Conservation penalty + |Tm − Tm_opt|.
    cost: f32,
    /// Footprint on the forward strand, `[from, to)`.
    from: usize,
    to: usize,
    /// Forward-strand coordinate of the primer's 3' end.
    end3: usize,
}

#[derive(Debug, Clone, Copy)]
struct PoolItem {
    key1: f32,
    key2: f32,
    f_from: usize,
    f_to: usize,
    r_from: usize,
    r_to: usize,
    fi: usize,
    ri: usize,
    score: f32,
    tm_diff: f32,
    amplicon: usize,
}

impl PoolItem {
    fn order(&self, o: &Self) -> Ordering {
        self.key1
            .total_cmp(&o.key1)
            .then_with(|| self.key2.total_cmp(&o.key2))
            .then_with(|| self.f_from.cmp(&o.f_from))
            .then_with(|| self.f_to.cmp(&o.f_to))
            .then_with(|| self.r_from.cmp(&o.r_from))
            .then_with(|| self.r_to.cmp(&o.r_to))
    }
}

impl PartialEq for PoolItem {
    fn eq(&self, o: &Self) -> bool {
        self.order(o) == Ordering::Equal
    }
}
impl Eq for PoolItem {}
impl PartialOrd for PoolItem {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for PoolItem {
    fn cmp(&self, o: &Self) -> Ordering {
        self.order(o)
    }
}

/// Smallest block of the [`BlockTree`]; shorter runs are scanned directly.
const LEAF: usize = 32;

/// Reverse candidates in the order of where their footprint ends, plus,
/// for every aligned block of that order (32, 64, 128, ... candidates), the
/// same candidates sorted by a ranking key. Any range of end positions is
/// covered by O(log n) such blocks, so the partners of a forward primer that
/// give a product of the allowed size can be visited best first.
struct BlockTree {
    order: Vec<u32>,
    ends: Vec<usize>,
    levels: Vec<Vec<u32>>,
}

impl BlockTree {
    fn new(rev: &[Candidate], key: impl Fn(&Candidate, &Candidate) -> Ordering) -> Self {
        let mut order: Vec<u32> = (0..rev.len() as u32).collect();
        order.sort_by(|&a, &b| {
            let (x, y) = (&rev[a as usize], &rev[b as usize]);
            x.to.cmp(&y.to).then(x.from.cmp(&y.from))
        });
        let ends = order.iter().map(|&i| rev[i as usize].to).collect();
        let n = order.len();
        let cmp = |a: &u32, b: &u32| key(&rev[*a as usize], &rev[*b as usize]);
        let mut first = order.clone();
        for chunk in first.chunks_mut(LEAF) {
            chunk.sort_by(cmp);
        }
        let mut levels = vec![first];
        let mut size = LEAF;
        while size < n {
            let prev = levels.last().expect("at least one level");
            let mut next = Vec::with_capacity(n);
            let mut i = 0;
            while i < n {
                let mid = (i + size).min(n);
                let end = (i + 2 * size).min(n);
                let (mut a, mut b) = (i, mid);
                while a < mid && b < end {
                    if cmp(&prev[b], &prev[a]) == Ordering::Less {
                        next.push(prev[b]);
                        b += 1;
                    } else {
                        next.push(prev[a]);
                        a += 1;
                    }
                }
                next.extend_from_slice(&prev[a..mid]);
                next.extend_from_slice(&prev[b..end]);
                i = end;
            }
            levels.push(next);
            size *= 2;
        }
        BlockTree { order, ends, levels }
    }

    /// Positions in `order` of the candidates ending in `[lo, hi]`.
    fn range(&self, lo: usize, hi: Option<usize>) -> (usize, usize) {
        let a = self.ends.partition_point(|&e| e < lo);
        let b = match hi {
            Some(h) => self.ends.partition_point(|&e| e <= h),
            None => self.ends.len(),
        };
        (a, b.max(a))
    }

    /// Call `f` on pieces that together cover `order[lo..hi]`: blocks sorted
    /// by the key (`true`) and short unsorted runs at the edges (`false`).
    fn visit(&self, lo: usize, hi: usize, mut f: impl FnMut(&[u32], bool)) {
        if lo >= hi {
            return;
        }
        let mut a = lo;
        let aligned = (a.div_ceil(LEAF) * LEAF).min(hi);
        if a < aligned {
            f(&self.order[a..aligned], false);
            a = aligned;
        }
        while a + LEAF <= hi {
            let mut k = 0;
            while k + 1 < self.levels.len() {
                let s = LEAF << (k + 1);
                if a % s == 0 && a + s <= hi {
                    k += 1;
                } else {
                    break;
                }
            }
            let s = LEAF << k;
            f(&self.levels[k][a..a + s], true);
            a += s;
        }
        if a < hi {
            f(&self.order[a..hi], false);
        }
    }
}

/// Bounded pool of the best pairs: a max-heap with the worst pair on top.
struct Pool {
    heap: BinaryHeap<PoolItem>,
    cap: usize,
}

impl Pool {
    fn new(cap: usize) -> Self {
        Pool { heap: BinaryHeap::with_capacity(cap + 1), cap }
    }
    fn full(&self) -> bool {
        self.heap.len() >= self.cap
    }
    /// The key a new pair must not exceed to enter, once the pool is full.
    fn bar(&self) -> Option<f32> {
        if self.full() {
            self.heap.peek().map(|w| w.key1)
        } else {
            None
        }
    }
    fn offer(&mut self, item: PoolItem) {
        if self.heap.len() < self.cap {
            self.heap.push(item);
        } else if let Some(worst) = self.heap.peek() {
            if item < *worst {
                self.heap.pop();
                self.heap.push(item);
            }
        }
    }
}

pub struct PrimerCalculator;

impl PrimerCalculator {
    /// Candidate windows on the forward strand (2.x API).
    pub fn generate_candidates(seq: &DataSequence, settings: &Settings) -> Vec<Primer> {
        Self::generate_candidates_on(seq, settings, Strand::Forward)
    }

    /// The first and last-plus-one window start on `strand`, and the
    /// lengths to try, or `None` when no window fits. For the reverse strand
    /// `seq` is the reverse complement and the design region `[b, e)` maps
    /// to `[n−e, n−b)`.
    fn window_bounds(n: usize, settings: &Settings, strand: Strand) -> Option<(usize, usize, usize, usize)> {
        let min_len = settings.min_primer_len.max(1);
        let max_len = settings.max_primer_len.min(n);
        if min_len > max_len {
            return None;
        }
        let b = settings.begin_nucleotide.unwrap_or(0).min(n);
        let e = settings.end_nucleotide.unwrap_or(n).min(n);
        if e <= b {
            return None;
        }
        let (start, end) = match strand {
            Strand::Forward => (b, e),
            Strand::Reverse => (n - e, n - b),
        };
        if end < start + min_len {
            return None;
        }
        Some((start, end, min_len, max_len))
    }

    /// Candidate windows on `strand`. For the reverse strand `seq` is the
    /// reverse complement and the design region `[b, e)` maps to `[n−e, n−b)`.
    pub fn generate_candidates_on(seq: &DataSequence, settings: &Settings, strand: Strand) -> Vec<Primer> {
        let Some((start, end, min_len, max_len)) = Self::window_bounds(seq.len(), settings, strand) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for len in min_len..=max_len {
            if end < start + len {
                continue;
            }
            for s in start..=(end - len) {
                out.push(Primer::new(s, len));
            }
        }
        out
    }

    /// Screen candidates and attach their conservation penalty (2.x API).
    pub fn filter_and_score(candidates: &[Primer], seq: &DataSequence, settings: &Settings) -> Vec<Primer> {
        candidates
            .iter()
            .filter_map(|p| {
                let w = seq.sub_seq(p.index, p.length);
                filters::screen(&w, settings).ok().map(|_| p.with_penalty(w.penalty()))
            })
            .collect()
    }

    /// Screen every window on one strand, in parallel chunks of window
    /// starts. Chunks are merged in order, so the result is the same for any
    /// number of threads.
    fn screen_strand(seq: &DataSequence, settings: &Settings, strand: Strand) -> (Vec<Candidate>, StrandStats) {
        let n = seq.len();
        let Some((start, end, min_len, max_len)) = Self::window_bounds(n, settings, strand) else {
            return (Vec::new(), StrandStats::default());
        };
        const CHUNK: usize = 2048;
        let last = end - min_len; // last window start, inclusive
        let chunks: Vec<(usize, usize)> =
            (start..=last).step_by(CHUNK).map(|a| (a, (a + CHUNK).min(last + 1))).collect();
        let opt = settings.tm_optimum();
        let nodes = seq.nodes();
        let parts: Vec<(Vec<Candidate>, StrandStats)> = chunks
            .par_iter()
            .map(|&(a, z)| {
                let mut sc = Scratch::default();
                let mut stats = StrandStats::default();
                let mut out = Vec::new();
                for s in a..z {
                    for len in min_len..=max_len {
                        if s + len > end {
                            break;
                        }
                        stats.considered += 1;
                        let w = &nodes[s..s + len];
                        if w[0].is_gap() || w[len - 1].is_gap() {
                            stats.deletions += 1;
                            continue;
                        }
                        match filters::screen_nodes(w, settings, &mut sc) {
                            Err(r) => stats.record(r),
                            Ok(tm) => {
                                stats.ok += 1;
                                let penalty = measure::penalty(w);
                                let (from, to) = match strand {
                                    Strand::Forward => (s, s + len),
                                    Strand::Reverse => (n - s - len, n - s),
                                };
                                let end3 = match strand {
                                    Strand::Forward => to - 1,
                                    Strand::Reverse => from,
                                };
                                out.push(Candidate {
                                    primer: Primer::new(s, len).with_penalty(penalty),
                                    tm,
                                    cost: penalty + (tm - opt).abs(),
                                    from,
                                    to,
                                    end3,
                                });
                            }
                        }
                    }
                }
                (out, stats)
            })
            .collect();
        let mut all = Vec::with_capacity(parts.iter().map(|p| p.0.len()).sum());
        let mut stats = StrandStats::default();
        for (c, st) in parts {
            all.extend(c);
            stats.add(&st);
        }
        (all, stats)
    }

    /// Full pipeline with statistics. Returns `Ok` with an empty pair list
    /// (and the statistics explaining why) when nothing passes.
    pub fn design(seq: &DataSequence, rev_seq: &DataSequence, settings: &Settings) -> Result<DesignOutcome, String> {
        settings.validate()?;
        if seq.len() != rev_seq.len() {
            return Err("forward and reverse templates must have the same length".into());
        }
        let ((fwd, fs), (rev, rs)) = rayon::join(
            || Self::screen_strand(seq, settings, Strand::Forward),
            || Self::screen_strand(rev_seq, settings, Strand::Reverse),
        );
        let mut explain = Explain { forward: fs, reverse: rs, pairs: PairStats::default() };
        let pairs = if fwd.is_empty() || rev.is_empty() {
            Vec::new()
        } else {
            Self::select_pairs(&fwd, &rev, seq, rev_seq, settings, &mut explain.pairs)
        };
        Ok(DesignOutcome { pairs, explain })
    }

    /// Exact counts of the combinations that overlap, or give a product
    /// that is too short or too long (in that order of precedence).
    fn count_position_failures(fwd: &[Candidate], rev: &[Candidate], settings: &Settings, stats: &mut PairStats) {
        stats.combinations = fwd.len() as u64 * rev.len() as u64;
        let mut froms: Vec<usize> = rev.iter().map(|r| r.from).collect();
        froms.sort_unstable();
        let mut by_len: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for r in rev {
            by_len.entry(r.to - r.from).or_default().push(r.from);
        }
        for v in by_len.values_mut() {
            v.sort_unstable();
        }
        let min_amp = settings.min_amplicon_len as i64;
        let (mut overlap, mut short, mut long) = (0u64, 0u64, 0u64);
        for f in fwd {
            overlap += froms.partition_point(|&x| x < f.to) as u64;
            for (&len, fr) in &by_len {
                let len = len as i64;
                let lo = f.to as i64; // non-overlapping partners start here
                let total_from_lo = fr.len() - fr.partition_point(|&x| (x as i64) < lo);
                // too short: r.to < f.from + min, i.e. r.from < f.from + min - len
                let short_hi = f.from as i64 + min_amp - len;
                let n_short = if short_hi > lo {
                    fr.partition_point(|&x| (x as i64) < short_hi) - fr.partition_point(|&x| (x as i64) < lo)
                } else {
                    0
                };
                short += n_short as u64;
                if let Some(max) = settings.max_amplicon_len {
                    // too long: r.to > f.from + max, i.e. r.from > f.from + max - len
                    let long_lo = (f.from as i64 + max as i64 - len + 1).max(lo).max(short_hi);
                    let n_long = fr.len() - fr.partition_point(|&x| (x as i64) < long_lo);
                    long += n_long.min(total_from_lo) as u64;
                }
            }
        }
        stats.overlap = overlap;
        stats.product_too_short = short;
        stats.product_too_long = long;
    }

    /// The best `cap` pairs that satisfy every pair constraint except the
    /// final dimer and near-duplicate checks, and whether the pool filled up.
    fn search_pairs(
        fwd: &[Candidate],
        rev: &[Candidate],
        settings: &Settings,
        cap: usize,
        stats: &mut PairStats,
    ) -> (Vec<PoolItem>, bool) {
        let mut pool = Pool::new(cap);
        let (mut examined, mut tm_fail) = (0u64, 0u64);
        let min_amp = settings.min_amplicon_len;
        let max_amp = settings.max_amplicon_len;
        let max_dtm = settings.max_tm_diff;
        // Partners must end inside [f.from + min, f.from + max] and after f.
        let window = |f: &Candidate| ((f.from + min_amp).max(f.to + 1), max_amp.map(|m| f.from + m));
        let item = |fi: usize, ri: usize, key1: f32, key2: f32, score: f32, tm_diff: f32| {
            let (f, r) = (&fwd[fi], &rev[ri]);
            PoolItem {
                key1,
                key2,
                f_from: f.from,
                f_to: f.to,
                r_from: r.from,
                r_to: r.to,
                fi,
                ri,
                score,
                tm_diff,
                amplicon: r.to - f.from,
            }
        };
        match settings.rank_by {
            RankBy::Score => {
                let tree = BlockTree::new(rev, |a, b| {
                    a.cost.total_cmp(&b.cost).then(a.from.cmp(&b.from)).then(a.to.cmp(&b.to))
                });
                let mut f_order: Vec<usize> = (0..fwd.len()).collect();
                f_order.sort_by(|&a, &b| {
                    let (x, y) = (&fwd[a], &fwd[b]);
                    x.cost.total_cmp(&y.cost).then(x.from.cmp(&y.from)).then(x.to.cmp(&y.to))
                });
                let min_r_cost = rev.iter().map(|r| r.cost).fold(f32::INFINITY, f32::min);
                for &fi in &f_order {
                    let f = &fwd[fi];
                    if let Some(bar) = pool.bar() {
                        if f.cost + min_r_cost > bar {
                            break;
                        }
                    }
                    let (lo_end, hi_end) = window(f);
                    let (lo, hi) = tree.range(lo_end, hi_end);
                    tree.visit(lo, hi, |slice, sorted| {
                        for &ri in slice {
                            let ri = ri as usize;
                            let r = &rev[ri];
                            let lower = f.cost + r.cost;
                            if let Some(bar) = pool.bar() {
                                if lower > bar {
                                    if sorted {
                                        break;
                                    }
                                    continue;
                                }
                            }
                            if r.from < f.to {
                                continue;
                            }
                            examined += 1;
                            let tm_diff = (f.tm - r.tm).abs();
                            if tm_diff > max_dtm {
                                tm_fail += 1;
                                continue;
                            }
                            let score = lower + tm_diff;
                            pool.offer(item(fi, ri, score, 0.0, score, tm_diff));
                        }
                    });
                }
            }
            RankBy::TmDiff => {
                let tree =
                    BlockTree::new(rev, |a, b| a.tm.total_cmp(&b.tm).then(a.from.cmp(&b.from)).then(a.to.cmp(&b.to)));
                for (fi, f) in fwd.iter().enumerate() {
                    let (lo_end, hi_end) = window(f);
                    let (lo, hi) = tree.range(lo_end, hi_end);
                    tree.visit(lo, hi, |slice, sorted| {
                        let theta0 = pool.bar().unwrap_or(max_dtm).min(max_dtm);
                        let first =
                            if sorted { slice.partition_point(|&ri| f.tm - rev[ri as usize].tm > theta0) } else { 0 };
                        for &ri in &slice[first..] {
                            let ri = ri as usize;
                            let r = &rev[ri];
                            let theta = pool.bar().unwrap_or(max_dtm).min(max_dtm);
                            if sorted && r.tm - f.tm > theta {
                                break;
                            }
                            let tm_diff = (f.tm - r.tm).abs();
                            if tm_diff > theta || r.from < f.to {
                                continue;
                            }
                            examined += 1;
                            let score = f.cost + r.cost + tm_diff;
                            pool.offer(item(fi, ri, tm_diff, f.primer.penalty + r.primer.penalty, score, tm_diff));
                        }
                    });
                }
            }
        }
        let saturated = pool.full();
        let mut items = pool.heap.into_vec();
        items.sort();
        stats.examined = examined;
        // With nothing found the search was exhaustive, so every combination
        // of the right size and position failed on ΔTm.
        stats.tm_diff = if items.is_empty() { stats.position_valid() } else { tm_fail };
        (items, saturated)
    }

    fn select_pairs(
        fwd: &[Candidate],
        rev: &[Candidate],
        seq: &DataSequence,
        rev_seq: &DataSequence,
        settings: &Settings,
        stats: &mut PairStats,
    ) -> Vec<PrimerPair> {
        let k = settings.max_primers_output.max(1);
        let mut cap = (k * 200).clamp(2000, 20000);
        Self::count_position_failures(fwd, rev, settings, stats);
        loop {
            let (pool, saturated) = Self::search_pairs(fwd, rev, settings, cap, stats);
            let d = settings.min_three_prime_distance;
            let mut chosen: Vec<(usize, usize)> = Vec::new();
            let mut out = Vec::new();
            let (mut similar, mut dimer) = (0u64, 0u64);
            for it in pool {
                if out.len() >= k {
                    break;
                }
                let (f, r) = (&fwd[it.fi], &rev[it.ri]);
                if d > 0 && chosen.iter().any(|&(fe, re)| fe.abs_diff(f.end3) < d && re.abs_diff(r.end3) < d) {
                    similar += 1;
                    continue;
                }
                let fw = seq.sub_seq(f.primer.index, f.primer.length);
                let rw = rev_seq.sub_seq(r.primer.index, r.primer.length);
                let c = Thermodynamics::complementarity(&fw, &rw);
                if c.any > settings.max_pair_any || c.end > settings.max_pair_end {
                    dimer += 1;
                    continue;
                }
                chosen.push((f.end3, r.end3));
                out.push(PrimerPair::new(f.primer, r.primer, it.amplicon, it.tm_diff, f.tm, r.tm).with_score(it.score));
            }
            // Too many held pairs were near-duplicates or dimers: look again
            // with a bigger pool (rare; the pool holds 2000+ pairs).
            if out.len() < k && saturated && cap < (1 << 20) {
                cap *= 8;
                continue;
            }
            stats.too_similar = similar;
            stats.pair_complementarity = dimer;
            stats.ok = out.len() as u64;
            return out;
        }
    }

    /// Pair already-chosen primers (2.x API). Unlike [`Self::design`] this
    /// returns every valid combination, so keep the inputs small.
    pub fn pair_primers(
        fwd_primers: &[Primer],
        rev_primers: &[Primer],
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
        settings: &Settings,
    ) -> Vec<PrimerPair> {
        let n = fwd_seq.len();
        let opt = settings.tm_optimum();
        let tm = |s: &DataSequence, p: &Primer| Thermodynamics::calculate_tm(&s.sub_seq(p.index, p.length), settings);
        let fwd_tms: Vec<f32> = fwd_primers.iter().map(|p| tm(fwd_seq, p)).collect();
        let rev_tms: Vec<f32> = rev_primers.iter().map(|p| tm(rev_seq, p)).collect();
        let mut pairs = Vec::new();
        for (i, f) in fwd_primers.iter().enumerate() {
            let f_to = f.index + f.length;
            for (j, r) in rev_primers.iter().enumerate() {
                if r.index + r.length > n {
                    continue;
                }
                let r_from = n - r.index - r.length;
                let r_to = n - r.index;
                if r_from < f_to {
                    continue;
                }
                let amplicon = r_to - f.index;
                if amplicon < settings.min_amplicon_len {
                    continue;
                }
                if let Some(max) = settings.max_amplicon_len {
                    if amplicon > max {
                        continue;
                    }
                }
                let tm_diff = (fwd_tms[i] - rev_tms[j]).abs();
                if tm_diff > settings.max_tm_diff {
                    continue;
                }
                let score = match settings.rank_by {
                    RankBy::Score => {
                        f.penalty + r.penalty + (fwd_tms[i] - opt).abs() + (rev_tms[j] - opt).abs() + tm_diff
                    }
                    RankBy::TmDiff => tm_diff,
                };
                pairs.push(PrimerPair::new(*f, *r, amplicon, tm_diff, fwd_tms[i], rev_tms[j]).with_score(score));
            }
        }
        pairs.sort();
        pairs
    }

    /// Full pipeline (2.x API): errors when no pair survives.
    pub fn design_primers(
        seq: &DataSequence,
        rev_seq: &DataSequence,
        settings: &Settings,
    ) -> Result<Vec<PrimerPair>, String> {
        let outcome = Self::design(seq, rev_seq, settings)?;
        if outcome.pairs.is_empty() {
            let e = &outcome.explain;
            return Err(if e.forward.ok == 0 {
                "No forward primers passed filters".into()
            } else if e.reverse.ok == 0 {
                "No reverse primers passed filters".into()
            } else {
                "Candidate primers exist but no valid pairs (check amplicon size / Tm diff)".into()
            });
        }
        Ok(outcome.pairs)
    }

    /// Diagnostic analysis of one primer (`degenprime test`, "Check a primer").
    pub fn analyze_primer(seq: &DataSequence, settings: &Settings) -> PrimerAnalysis {
        let tm = Thermodynamics::calculate_tm(seq, settings);
        let range = Thermodynamics::tm_range(seq, settings, 4096);
        let sc = Thermodynamics::self_complementarity(seq);
        PrimerAnalysis {
            sequence: seq.codes(),
            length: seq.len(),
            actual_size: seq.actual_size(),
            tm,
            gc_content: seq.gc_ratio() * 100.0,
            degeneracy_count: seq.degeneracy_count(),
            gap_count: seq.gap_count(),
            penalty: seq.penalty(),
            hairpin_score: Thermodynamics::hairpin_score(seq),
            filters: filters::apply_all(seq, settings),
            fold: seq.fold(),
            tm_min: range.map(|r| r.0),
            tm_max: range.map(|r| r.1),
            self_any: sc.any,
            self_end: sc.end,
            end_stability: Thermodynamics::end_stability(seq),
            self_dimer: if sc.end > 0.0 { Thermodynamics::duplex_diagram(seq, seq, &sc) } else { None },
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerAnalysis {
    pub sequence: String,
    pub length: usize,
    pub actual_size: usize,
    pub tm: f32,
    pub gc_content: f32,
    pub degeneracy_count: usize,
    pub gap_count: usize,
    pub penalty: f32,
    pub hairpin_score: usize,
    pub filters: Vec<FilterResult>,
    /// Number of distinct sequences the degenerate primer stands for.
    #[serde(default)]
    pub fold: u64,
    /// Lowest / highest Tm across all variants (when ≤ 4096 variants).
    #[serde(default)]
    pub tm_min: Option<f32>,
    #[serde(default)]
    pub tm_max: Option<f32>,
    /// Primer3-style self-complementarity scores.
    #[serde(default)]
    pub self_any: f32,
    #[serde(default)]
    pub self_end: f32,
    /// −ΔG(37 °C) of the last five 3' bases, kcal/mol.
    #[serde(default)]
    pub end_stability: f32,
    /// Drawing of the strongest 3' self-dimer, if any.
    #[serde(default)]
    pub self_dimer: Option<String>,
}

impl PrimerAnalysis {
    pub fn all_passed(&self) -> bool {
        self.filters.iter().all(|f| f.passed)
    }

    /// Pretty-print for the CLI.
    pub fn report(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!("Primer:     5'-{}-3'\n", self.sequence));
        s.push_str(&format!("Length:     {} nt (non-gap: {})\n", self.length, self.actual_size));
        match (self.tm_min, self.tm_max) {
            (Some(lo), Some(hi)) if self.fold > 1 => {
                s.push_str(&format!("Tm:         {:.2} °C (variants {:.1}–{:.1} °C)\n", self.tm, lo, hi))
            }
            _ => s.push_str(&format!("Tm:         {:.2} °C\n", self.tm)),
        }
        s.push_str(&format!("GC:         {:.1}%\n", self.gc_content));
        s.push_str(&format!(
            "Degenerate: {} position{}, {}-fold\n",
            self.degeneracy_count,
            if self.degeneracy_count == 1 { "" } else { "s" },
            self.fold
        ));
        s.push_str(&format!("Self-compl: any {:.2}, 3' {:.2} (Primer3-style scores)\n", self.self_any, self.self_end));
        s.push_str(&format!("3' stability: {:.1} kcal/mol (last 5 bases)\n", self.end_stability));
        s.push_str(&format!("Hairpin:    {} (legacy 0–3 heuristic)\n", self.hairpin_score));
        s.push_str(&format!("Penalty:    {:.2}\n", self.penalty));
        s.push_str("\nFilters:\n");
        for f in &self.filters {
            let mark = if f.passed { "PASS" } else { "FAIL" };
            s.push_str(&format!("  [{mark}] {}", f.name));
            if let Some(r) = &f.reason {
                s.push_str(&format!(" — {r}"));
            }
            s.push('\n');
        }
        if let Some(d) = &self.self_dimer {
            s.push_str("\nStrongest 3' self-dimer:\n");
            s.push_str(d);
            s.push('\n');
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::SequenceReader;
    use crate::models::DataNode;

    fn ds_from(s: &str) -> DataSequence {
        DataSequence::new(s.chars().map(DataNode::from_base).collect())
    }

    fn example() -> (DataSequence, DataSequence) {
        let text = include_str!("../data/example_anfD.clust");
        let (seqs, _) = SequenceReader::parse(text).unwrap();
        let cons = SequenceReader::alignment_to_consensus(&seqs).unwrap();
        let rc = cons.reverse_complement();
        (cons, rc)
    }

    #[test]
    fn candidate_count_matches_expectation() {
        let seq = ds_from(&"ATGC".repeat(25));
        let cands = PrimerCalculator::generate_candidates(&seq, &Settings::default());
        assert_eq!(cands.len(), 81 + 80 + 79);
    }

    #[test]
    fn pair_primers_uses_rev_complement_coords() {
        let fwd = ds_from(&"A".repeat(200));
        let rev = fwd.reverse_complement();
        let settings = Settings { min_tm: 0.0, max_tm: 100.0, ..Settings::default() };
        let pairs = PrimerCalculator::pair_primers(&[Primer::new(0, 20)], &[Primer::new(0, 20)], &fwd, &rev, &settings);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].amplicon_size, 200);
    }

    #[test]
    fn pair_rejects_overlap() {
        let fwd = ds_from(&"A".repeat(50));
        let rev = fwd.reverse_complement();
        let pairs = PrimerCalculator::pair_primers(
            &[Primer::new(20, 20)],
            &[Primer::new(20, 20)],
            &fwd,
            &rev,
            &Settings::default(),
        );
        assert!(pairs.is_empty());
    }

    #[test]
    fn reverse_candidates_respect_the_region() {
        let seq = ds_from(&"ACGTTGCAAGCTTGCA".repeat(20));
        let n = seq.len();
        let s = Settings::default().with_region(0, 150);
        let cands = PrimerCalculator::generate_candidates_on(&seq.reverse_complement(), &s, Strand::Reverse);
        assert!(!cands.is_empty());
        for p in cands {
            let to = n - p.index;
            assert!(to <= 150, "reverse primer footprint ends at {to}, outside [0,150)");
        }
    }

    #[test]
    fn example_design_is_deterministic() {
        let (cons, rc) = example();
        let a = PrimerCalculator::design(&cons, &rc, &Settings::default()).unwrap();
        assert!(!a.pairs.is_empty(), "explain: {:?}", a.explain);
        for _ in 0..3 {
            let b = PrimerCalculator::design(&cons, &rc, &Settings::default()).unwrap();
            assert_eq!(a.pairs, b.pairs);
        }
        let scores: Vec<f32> = a.pairs.iter().map(|p| p.score).collect();
        assert!(scores.windows(2).all(|w| w[0] <= w[1]), "pairs must be sorted by score");
    }

    #[test]
    fn example_design_respects_region() {
        let (cons, rc) = example();
        let s = Settings::default().with_region(0, 700);
        let out = PrimerCalculator::design(&cons, &rc, &s).unwrap();
        assert!(!out.pairs.is_empty());
        for p in &out.pairs {
            let rev_left = p.reverse_fwd_strand_position(cons.len());
            assert!(rev_left + p.reverse.length <= 700, "reverse primer outside region");
            assert!(p.forward.index + p.forward.length <= 700);
        }
    }

    #[test]
    fn returned_pairs_are_not_near_duplicates() {
        let (cons, rc) = example();
        let out = PrimerCalculator::design(&cons, &rc, &Settings::default()).unwrap();
        let n = cons.len();
        for (i, a) in out.pairs.iter().enumerate() {
            for b in out.pairs.iter().skip(i + 1) {
                let fa = a.forward.index + a.forward.length;
                let fb = b.forward.index + b.forward.length;
                let ra = a.reverse_fwd_strand_position(n);
                let rb = b.reverse_fwd_strand_position(n);
                assert!(fa.abs_diff(fb) >= 3 || ra.abs_diff(rb) >= 3);
            }
        }
    }

    #[test]
    fn legacy_ranking_still_available() {
        let (cons, rc) = example();
        let s = Settings::default().with_rank_by(RankBy::TmDiff);
        let out = PrimerCalculator::design(&cons, &rc, &s).unwrap();
        assert!(!out.pairs.is_empty());
        let diffs: Vec<f32> = out.pairs.iter().map(|p| p.tm_diff).collect();
        assert!(diffs.windows(2).all(|w| w[0] <= w[1] + 1e-6));
    }

    /// Every valid pair, by brute force, in ranking order.
    fn all_pairs(fwd: &[Candidate], rev: &[Candidate], s: &Settings) -> Vec<PoolItem> {
        let mut v = Vec::new();
        for (fi, f) in fwd.iter().enumerate() {
            for (ri, r) in rev.iter().enumerate() {
                if r.from < f.to {
                    continue;
                }
                let amplicon = r.to - f.from;
                if amplicon < s.min_amplicon_len || s.max_amplicon_len.is_some_and(|m| amplicon > m) {
                    continue;
                }
                let tm_diff = (f.tm - r.tm).abs();
                if tm_diff > s.max_tm_diff {
                    continue;
                }
                let score = f.cost + r.cost + tm_diff;
                let (key1, key2) = match s.rank_by {
                    RankBy::Score => (score, 0.0),
                    RankBy::TmDiff => (tm_diff, f.primer.penalty + r.primer.penalty),
                };
                v.push(PoolItem {
                    key1,
                    key2,
                    f_from: f.from,
                    f_to: f.to,
                    r_from: r.from,
                    r_to: r.to,
                    fi,
                    ri,
                    score,
                    tm_diff,
                    amplicon,
                });
            }
        }
        v.sort();
        v
    }

    fn random_candidates(n: usize, len: usize, seed: u64) -> Vec<Candidate> {
        let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        (0..n)
            .map(|_| {
                let l = 18 + (next() % 5) as usize;
                let from = (next() % (len - l) as u64) as usize;
                // Coarse values make ties common, which is what breaks a
                // careless pruning rule.
                let tm = 58.0 + (next() % 9) as f32 * 0.5;
                let penalty = (next() % 4) as f32 * 10.0;
                Candidate {
                    primer: Primer::new(from, l).with_penalty(penalty),
                    tm,
                    cost: penalty + (tm - 60.0).abs(),
                    from,
                    to: from + l,
                    end3: from,
                }
            })
            .collect()
    }

    #[test]
    fn pair_search_matches_brute_force() {
        for seed in 1..40u64 {
            let len = 300 + (seed as usize * 37) % 900;
            let fwd = random_candidates(20 + (seed as usize * 7) % 90, len, seed);
            let rev = random_candidates(20 + (seed as usize * 11) % 90, len, seed + 1000);
            for rank in [RankBy::Score, RankBy::TmDiff] {
                for (min_amp, max_amp) in [(100, None), (50, Some(200)), (150, Some(151)), (10, Some(5000))] {
                    for cap in [1usize, 7, 50, 100_000] {
                        let mut s = Settings::default().with_rank_by(rank);
                        s.min_amplicon_len = min_amp;
                        s.max_amplicon_len = max_amp;
                        s.max_tm_diff = 1.5;
                        let mut stats = PairStats::default();
                        PrimerCalculator::count_position_failures(&fwd, &rev, &s, &mut stats);
                        let (got, _) = PrimerCalculator::search_pairs(&fwd, &rev, &s, cap, &mut stats);
                        let mut want = all_pairs(&fwd, &rev, &s);
                        want.truncate(cap);
                        let key =
                            |v: &[PoolItem]| v.iter().map(|p| (p.f_from, p.f_to, p.r_from, p.r_to)).collect::<Vec<_>>();
                        assert_eq!(key(&got), key(&want), "seed {seed} {rank:?} {min_amp}-{max_amp:?} cap {cap}");
                    }
                }
            }
        }
    }

    #[test]
    fn position_failure_counts_are_exact() {
        for seed in 1..25u64 {
            let fwd = random_candidates(60, 800, seed);
            let rev = random_candidates(70, 800, seed + 77);
            for (min_amp, max_amp) in [(100, None), (60, Some(250)), (300, Some(300))] {
                let s = Settings { min_amplicon_len: min_amp, max_amplicon_len: max_amp, ..Settings::default() };
                let mut stats = PairStats::default();
                PrimerCalculator::count_position_failures(&fwd, &rev, &s, &mut stats);
                let (mut overlap, mut short, mut long) = (0u64, 0u64, 0u64);
                for f in &fwd {
                    for r in &rev {
                        let amp = r.to as i64 - f.from as i64;
                        if r.from < f.to {
                            overlap += 1;
                        } else if amp < min_amp as i64 {
                            short += 1;
                        } else if max_amp.is_some_and(|m| amp > m as i64) {
                            long += 1;
                        }
                    }
                }
                assert_eq!(stats.combinations, (fwd.len() * rev.len()) as u64);
                assert_eq!((stats.overlap, stats.product_too_short, stats.product_too_long), (overlap, short, long));
            }
        }
    }

    #[test]
    fn results_do_not_depend_on_thread_count() {
        let (cons, rc) = example();
        let one = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        let many = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        let s = Settings::default().with_max_primers(20);
        let a = one.install(|| PrimerCalculator::design(&cons, &rc, &s).unwrap());
        let b = many.install(|| PrimerCalculator::design(&cons, &rc, &s).unwrap());
        assert_eq!(a.pairs, b.pairs);
        assert_eq!(a.explain.forward.ok, b.explain.forward.ok);
    }
}
