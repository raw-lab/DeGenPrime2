//! Primer quality filters.
//!
//! Each filter returns a `FilterResult` so the caller can build diagnostic
//! reports. All filters consult `Settings` directly so user parameters are
//! never silently ignored — this was a bug in v2 where `filter_temperature`
//! used hard-coded 50 nM / 50 mM regardless of what the user passed.

use crate::models::DataSequence;
use crate::settings::Settings;
use crate::thermodynamics::Thermodynamics;
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
        FilterResult {
            name: name.into(),
            passed: false,
            reason: Some(reason.into()),
        }
    }
}

/// Degeneracy filter: too many ambiguous IUPAC positions => fail.
pub fn filter_degeneracy(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let degen = seq.degeneracy_count();
    let denom = seq.actual_size().max(1);
    let ratio = degen as f32 / denom as f32;
    if ratio > settings.max_degeneracy_ratio {
        FilterResult::fail(
            "degeneracy",
            format!("{degen} degenerate positions / {denom} ({:.1}% > {:.1}%)",
                ratio * 100.0,
                settings.max_degeneracy_ratio * 100.0),
        )
    } else {
        FilterResult::pass("degeneracy")
    }
}

/// Deletion filter: too many gaps => fail.
pub fn filter_deletions(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let gaps = seq.gap_count();
    let ratio = gaps as f32 / seq.len().max(1) as f32;
    if ratio > settings.max_deletion_ratio {
        FilterResult::fail(
            "deletions",
            format!("{gaps} gaps / {} ({:.1}% > {:.1}%)",
                seq.len(),
                ratio * 100.0,
                settings.max_deletion_ratio * 100.0),
        )
    } else {
        FilterResult::pass("deletions")
    }
}

/// GC content filter: enforces the overall window and a 3' end cap.
pub fn filter_gc_content(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let gc = seq.gc_ratio();
    if gc < settings.min_gc_content || gc > settings.max_gc_content {
        return FilterResult::fail(
            "gc_content",
            format!(
                "GC {:.1}% outside [{:.0}%, {:.0}%]",
                gc * 100.0,
                settings.min_gc_content * 100.0,
                settings.max_gc_content * 100.0
            ),
        );
    }

    if seq.len() >= 5 {
        let end = seq.sub_seq(seq.len() - 5, 5);
        let end_gc = end.gc_ratio();
        if end_gc > settings.max_gc_end {
            return FilterResult::fail(
                "gc_content_3prime",
                format!(
                    "3' end GC {:.0}% > {:.0}%",
                    end_gc * 100.0,
                    settings.max_gc_end * 100.0
                ),
            );
        }
    }
    FilterResult::pass("gc_content")
}

/// Repetition filter: bans long homopolymer runs and obviously bad tandem
/// repeats (>= 3 dinucleotide repeats in a row, e.g. ATATAT).
///
/// The v2 version flagged ANY occurrence of `w[0]==w[2] && w[1]==w[3]`,
/// which rejected most primers because `ATGCAT` contains a window with that
/// pattern (`TGCATG` -> T==T, G==G means w[0]==w[2], w[1]==w[3]... actually
/// in `TGAT`, T==A? no. But many real primers do have one alternating
/// 4-mer somewhere). Here we require 3 repeats of the dimer in a row.
pub fn filter_repeats(seq: &DataSequence, _settings: &Settings) -> FilterResult {
    let codes: Vec<char> = seq.nodes().iter().map(|n| n.most_common()).collect();

    // 1) Homopolymer ≥5 fails.
    let mut run_len = 1;
    let mut run_char = codes.first().copied().unwrap_or(' ');
    for c in codes.iter().skip(1) {
        if *c == run_char {
            run_len += 1;
            if run_len >= 5 {
                return FilterResult::fail(
                    "repeats",
                    format!("homopolymer run of {run_len} '{c}'s"),
                );
            }
        } else {
            run_char = *c;
            run_len = 1;
        }
    }

    // 2) Tandem dinucleotide repeats: require ≥3 copies in a row (so 6 bases).
    if codes.len() >= 6 {
        for i in 0..(codes.len() - 5) {
            let a = codes[i];
            let b = codes[i + 1];
            if a == b { continue; }
            if codes[i + 2] == a && codes[i + 3] == b
               && codes[i + 4] == a && codes[i + 5] == b
            {
                return FilterResult::fail(
                    "repeats",
                    format!("tandem dinucleotide repeat {}{}", a, b),
                );
            }
        }
    }

    FilterResult::pass("repeats")
}

/// Complementary-ends filter: a quick heuristic for 3'-3' self-folding.
pub fn filter_complementary_ends(seq: &DataSequence, _settings: &Settings) -> FilterResult {
    if seq.has_complementary_ends() {
        FilterResult::fail(
            "complementary_ends",
            "first 3 and last 3 bases are complementary".to_string(),
        )
    } else {
        FilterResult::pass("complementary_ends")
    }
}

/// Melting-temperature filter: now properly uses the settings the user supplied.
pub fn filter_temperature(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let tm = Thermodynamics::calculate_tm(seq, settings);
    if tm < settings.min_tm || tm > settings.max_tm {
        FilterResult::fail(
            "temperature",
            format!(
                "Tm {:.1}°C outside [{:.1}, {:.1}]",
                tm, settings.min_tm, settings.max_tm
            ),
        )
    } else {
        FilterResult::pass("temperature")
    }
}

/// Apply every standard filter. Returns the full list, including passes,
/// so the caller can build a diagnostic report.
pub fn apply_all(seq: &DataSequence, settings: &Settings) -> Vec<FilterResult> {
    vec![
        filter_degeneracy(seq, settings),
        filter_deletions(seq, settings),
        filter_gc_content(seq, settings),
        filter_repeats(seq, settings),
        filter_complementary_ends(seq, settings),
        filter_temperature(seq, settings),
    ]
}

/// Fast yes/no version for the design pipeline.
pub fn passes_all(seq: &DataSequence, settings: &Settings) -> bool {
    [filter_degeneracy, filter_deletions, filter_gc_content,
     filter_repeats, filter_complementary_ends, filter_temperature]
        .iter()
        .all(|f| f(seq, settings).passed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::DataSequence;

    #[test]
    fn good_primer_passes_gc_filter() {
        let s = DataSequence::from_string("ATGCAATTGGCC"); // 50% GC
        assert!(filter_gc_content(&s, &Settings::default()).passed);
    }

    #[test]
    fn homopolymer_fails_repeats() {
        let s = DataSequence::from_string("ATGAAAAATGC"); // AAAAA
        assert!(!filter_repeats(&s, &Settings::default()).passed);
    }

    #[test]
    fn isolated_dinuc_passes_repeats() {
        // "ATAT" in a larger primer should NOT fail (this was the v2 bug)
        let s = DataSequence::from_string("ATATGCAATTGGCC");
        assert!(filter_repeats(&s, &Settings::default()).passed);
    }

    #[test]
    fn triple_dinuc_fails_repeats() {
        // 3 repeats of "AT" (=6 bases ATATAT) SHOULD fail
        let s = DataSequence::from_string("CGATATATGGCC");
        assert!(!filter_repeats(&s, &Settings::default()).passed);
    }

    #[test]
    fn temperature_filter_uses_settings() {
        // A primer that's OK at default settings becomes unacceptable if we
        // require a very narrow high range.
        let s = DataSequence::from_string("ATGCAATTGGCCGGAATTCC");
        let strict = Settings::default().with_tm_range(60.0, 60.5);
        let result = filter_temperature(&s, &strict);
        // Either passes or fails — we just need to confirm the filter ACTUALLY
        // reads the settings (it would have always passed in v2's broken version)
        let tm = Thermodynamics::calculate_tm(&s, &strict);
        assert_eq!(result.passed, tm >= 60.0 && tm <= 60.5);
    }
}
