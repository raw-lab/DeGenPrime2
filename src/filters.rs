//! Primer quality filters.
//!
//! [`screen`] is the fast path used by the design pipeline: it runs every
//! check in a fixed order, stops at the first failure and returns the Tm it
//! computed, so the pipeline never computes it twice. The `filter_*`
//! functions build human-readable [`FilterResult`]s for `degenprime test` and
//! the "Check a primer" page. All checks work on the oligo with gap columns
//! removed, and all of them read their limits from [`Settings`].

use crate::models::{measure, DataNode, DataSequence};
use crate::settings::Settings;
use crate::thermodynamics::{Thermodynamics, TmConsts};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FilterResult {
    pub name: String,
    pub passed: bool,
    pub reason: Option<String>,
}

impl FilterResult {
    pub fn pass(name: &str) -> Self {
        FilterResult { name: name.into(), passed: true, reason: None }
    }

    pub fn fail(name: &str, reason: impl Into<String>) -> Self {
        FilterResult { name: name.into(), passed: false, reason: Some(reason.into()) }
    }
}

/// The first check a candidate window failed (pipeline order).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rejection {
    Degeneracy,
    Deletions,
    GcContent,
    Gc3Prime,
    Repeats,
    ComplementaryEnds,
    TmLow,
    TmHigh,
    SelfComplementarity,
}

/// Reusable buffers for [`screen_nodes`]; keep one per thread.
#[derive(Debug, Default)]
pub struct Scratch {
    oligo: Vec<DataNode>,
    masks: Vec<u8>,
    /// Tm reaction terms, with the (primer, salt) concentrations they were
    /// worked out for.
    tm: Option<(u32, u32, TmConsts)>,
}

/// Run every check on a candidate window. `Ok(tm)` when it passes.
pub fn screen(window: &DataSequence, settings: &Settings) -> Result<f32, Rejection> {
    screen_nodes(window.nodes(), settings, &mut Scratch::default())
}

/// [`screen`] on a slice of the template, without allocating (the design
/// pipeline calls this for every window of the alignment). Same checks,
/// same order, same arithmetic.
pub fn screen_nodes(window: &[DataNode], settings: &Settings, sc: &mut Scratch) -> Result<f32, Rejection> {
    let (mut actual, mut degenerate, mut gaps) = (0usize, 0usize, 0usize);
    for n in window {
        if n.is_gap() {
            gaps += 1;
        } else {
            actual += 1;
            if !matches!(n.code(), 'A' | 'C' | 'G' | 'T') {
                degenerate += 1;
            }
        }
    }
    if degenerate as f32 / actual.max(1) as f32 > settings.max_degeneracy_ratio {
        return Err(Rejection::Degeneracy);
    }
    if gaps as f32 / window.len().max(1) as f32 > settings.max_deletion_ratio {
        return Err(Rejection::Deletions);
    }
    sc.oligo.clear();
    sc.oligo.extend(window.iter().copied().filter(|n| !n.is_gap()));
    let oligo = &sc.oligo[..];
    let gc = measure::gc_ratio(oligo);
    if gc < settings.min_gc_content || gc > settings.max_gc_content {
        return Err(Rejection::GcContent);
    }
    if oligo.len() >= 5 && measure::gc_ratio(&oligo[oligo.len() - 5..]) > settings.max_gc_end {
        return Err(Rejection::Gc3Prime);
    }
    if repeat_in(oligo.iter().map(|n| n.most_common()), settings.max_poly_x).is_some() {
        return Err(Rejection::Repeats);
    }
    if measure::complementary_ends(oligo) {
        return Err(Rejection::ComplementaryEnds);
    }
    let key = (settings.primer_concentration.to_bits(), settings.salt_concentration.to_bits());
    let consts = match sc.tm {
        Some((p, q, k)) if (p, q) == key => k,
        _ => {
            let k = TmConsts::new(settings);
            sc.tm = Some((key.0, key.1, k));
            k
        }
    };
    let tm = Thermodynamics::tm_of_nodes_with(oligo, &consts);
    if tm < settings.min_tm {
        return Err(Rejection::TmLow);
    }
    if tm > settings.max_tm {
        return Err(Rejection::TmHigh);
    }
    let c = Thermodynamics::self_complementarity_nodes(oligo, &mut sc.masks);
    if c.any > settings.max_self_any || c.end > settings.max_self_end {
        return Err(Rejection::SelfComplementarity);
    }
    Ok(tm)
}

/// The first repeat problem in a run of bases: a single-base run longer
/// than `max_poly_x`, or three tandem copies of a dinucleotide.
enum Repeat {
    Run(usize, char),
    Dinucleotide(char, char),
}

fn repeat_in(bases: impl Iterator<Item = char>, max_poly_x: usize) -> Option<Repeat> {
    // A sliding window of the last six bases is all either check needs.
    let mut last = ['\0'; 6];
    let mut seen = 0usize;
    let mut run = 1usize;
    let mut dinuc: Option<Repeat> = None;
    for c in bases {
        if seen > 0 {
            if c == last[5] {
                run += 1;
                if run > max_poly_x {
                    return Some(Repeat::Run(run, c));
                }
            } else {
                run = 1;
            }
        }
        last.rotate_left(1);
        last[5] = c;
        seen += 1;
        if dinuc.is_none() && seen >= 6 {
            let (a, b) = (last[0], last[1]);
            if a != b && last[2] == a && last[3] == b && last[4] == a && last[5] == b {
                dinuc = Some(Repeat::Dinucleotide(a, b));
            }
        }
    }
    dinuc
}

/// Homopolymer runs longer than `max_poly_x`, or three tandem copies of a
/// dinucleotide (ATATAT). A single isolated ATAT is fine.
pub fn repeat_problem(seq: &DataSequence, settings: &Settings) -> Option<String> {
    let bases = seq.nodes().iter().filter(|n| !n.is_gap()).map(|n| n.most_common());
    repeat_in(bases, settings.max_poly_x).map(|r| match r {
        Repeat::Run(run, c) => format!("homopolymer run of {run} '{c}'s (max {})", settings.max_poly_x),
        Repeat::Dinucleotide(a, b) => format!("tandem dinucleotide repeat {a}{b}{a}{b}{a}{b}"),
    })
}

pub fn filter_degeneracy(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let degen = seq.degeneracy_count();
    let denom = seq.actual_size().max(1);
    let ratio = degen as f32 / denom as f32;
    if ratio > settings.max_degeneracy_ratio {
        FilterResult::fail(
            "degeneracy",
            format!(
                "{degen} degenerate positions of {denom} ({:.1}% > {:.1}%)",
                ratio * 100.0,
                settings.max_degeneracy_ratio * 100.0
            ),
        )
    } else {
        FilterResult::pass("degeneracy")
    }
}

pub fn filter_deletions(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let gaps = seq.gap_count();
    let ratio = gaps as f32 / seq.len().max(1) as f32;
    if ratio > settings.max_deletion_ratio {
        FilterResult::fail(
            "deletions",
            format!(
                "{gaps} gap columns of {} ({:.1}% > {:.1}%)",
                seq.len(),
                ratio * 100.0,
                settings.max_deletion_ratio * 100.0
            ),
        )
    } else {
        FilterResult::pass("deletions")
    }
}

pub fn filter_gc_content(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let oligo = seq.ungapped();
    let gc = oligo.gc_ratio();
    if gc < settings.min_gc_content || gc > settings.max_gc_content {
        return FilterResult::fail(
            "gc_content",
            format!(
                "GC {:.1}% outside {:.0}–{:.0}%",
                gc * 100.0,
                settings.min_gc_content * 100.0,
                settings.max_gc_content * 100.0
            ),
        );
    }
    if oligo.len() >= 5 {
        let end_gc = oligo.sub_seq(oligo.len() - 5, 5).gc_ratio();
        if end_gc > settings.max_gc_end {
            return FilterResult::fail(
                "gc_content_3prime",
                format!("3' end GC {:.0}% > {:.0}% (last 5 bases)", end_gc * 100.0, settings.max_gc_end * 100.0),
            );
        }
    }
    FilterResult::pass("gc_content")
}

pub fn filter_repeats(seq: &DataSequence, settings: &Settings) -> FilterResult {
    match repeat_problem(&seq.ungapped(), settings) {
        Some(why) => FilterResult::fail("repeats", why),
        None => FilterResult::pass("repeats"),
    }
}

/// Complementary-ends filter: a quick heuristic for 5'/3' self-folding.
pub fn filter_complementary_ends(seq: &DataSequence, _settings: &Settings) -> FilterResult {
    if seq.has_complementary_ends() {
        FilterResult::fail("complementary_ends", "first 3 and last 3 bases are complementary")
    } else {
        FilterResult::pass("complementary_ends")
    }
}

pub fn filter_temperature(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let tm = Thermodynamics::calculate_tm(seq, settings);
    if tm < settings.min_tm || tm > settings.max_tm {
        FilterResult::fail(
            "temperature",
            format!("Tm {:.1} °C outside {:.1}–{:.1} °C", tm, settings.min_tm, settings.max_tm),
        )
    } else {
        FilterResult::pass("temperature")
    }
}

/// Primer3-style self-complementarity (self-dimers and hairpins).
pub fn filter_self_complementarity(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let c = Thermodynamics::self_complementarity(seq);
    if c.end > settings.max_self_end {
        FilterResult::fail(
            "self_complementarity",
            format!("3' self-complementarity {:.1} > {:.1}", c.end, settings.max_self_end),
        )
    } else if c.any > settings.max_self_any {
        FilterResult::fail(
            "self_complementarity",
            format!("self-complementarity {:.1} > {:.1}", c.any, settings.max_self_any),
        )
    } else {
        FilterResult::pass("self_complementarity")
    }
}

/// Every filter, passes included, for diagnostic reports.
pub fn apply_all(seq: &DataSequence, settings: &Settings) -> Vec<FilterResult> {
    vec![
        filter_degeneracy(seq, settings),
        filter_deletions(seq, settings),
        filter_gc_content(seq, settings),
        filter_repeats(seq, settings),
        filter_complementary_ends(seq, settings),
        filter_temperature(seq, settings),
        filter_self_complementarity(seq, settings),
    ]
}

/// Yes/no version of [`apply_all`].
pub fn passes_all(seq: &DataSequence, settings: &Settings) -> bool {
    screen(seq, settings).is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::DataSequence;

    #[test]
    fn good_primer_passes_gc_filter() {
        let s = DataSequence::from_string("ATGCAATTGGCC");
        assert!(filter_gc_content(&s, &Settings::default()).passed);
    }

    #[test]
    fn homopolymer_fails_repeats() {
        let s = DataSequence::from_string("ATGAAAAATGC");
        assert!(!filter_repeats(&s, &Settings::default()).passed);
    }

    #[test]
    fn poly_x_limit_is_configurable() {
        let s = DataSequence::from_string("ATGAAAAATGC");
        assert!(filter_repeats(&s, &Settings::default().with_max_poly_x(5)).passed);
    }

    #[test]
    fn isolated_dinuc_passes_repeats() {
        let s = DataSequence::from_string("ATATGCAATTGGCC");
        assert!(filter_repeats(&s, &Settings::default()).passed);
    }

    #[test]
    fn triple_dinuc_fails_repeats() {
        let s = DataSequence::from_string("CGATATATGGCC");
        assert!(!filter_repeats(&s, &Settings::default()).passed);
    }

    #[test]
    fn temperature_filter_uses_settings() {
        let s = DataSequence::from_string("ATGCAATTGGCCGGAATTCC");
        let strict = Settings::default().with_tm_range(60.0, 60.5);
        let result = filter_temperature(&s, &strict);
        let tm = Thermodynamics::calculate_tm(&s, &strict);
        assert_eq!(result.passed, (60.0..=60.5).contains(&tm));
    }

    #[test]
    fn palindromic_three_prime_end_is_rejected() {
        let s = DataSequence::from_string("ACGTAGCTAGGAATTC");
        assert!(!filter_self_complementarity(&s, &Settings::default()).passed);
    }

    #[test]
    fn deletions_filter_now_sees_gap_columns() {
        use crate::models::DataNode;
        let mut nodes: Vec<DataNode> = "ACGTACGTACGTACGTACGT".chars().map(DataNode::from_base).collect();
        for n in nodes.iter_mut().take(5) {
            *n = DataNode::new(vec!['A', '-', '-', '-']);
        }
        let s = DataSequence::new(nodes);
        assert_eq!(s.gap_count(), 5);
        assert!(!filter_deletions(&s, &Settings::default()).passed);
    }
}
