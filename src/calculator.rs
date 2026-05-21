//! Primer design pipeline.
//!
//! High-level flow:
//!   1. Generate all candidate primer windows (forward + reverse-complement).
//!   2. Run each candidate through `filters::passes_all`.
//!   3. Rank survivors by combined penalty (low conservation hurts).
//!   4. Pair forward with reverse, score by Tm difference.
//!
//! IMPORTANT FIX vs DeGenPrime2 v2: the reverse primer's `index` is on the
//! reverse-complement strand. To compute amplicon size and forward-strand
//! coordinates we use `template_len - rev.index - rev.length` as the 5' position
//! of the reverse primer on the forward strand.

use crate::filters::{self, FilterResult};
use crate::models::{DataSequence, Primer, PrimerPair};
use crate::settings::Settings;
use crate::thermodynamics::Thermodynamics;

pub struct PrimerCalculator;

impl PrimerCalculator {
    /// Generate every candidate primer in `[min_len, max_len]` across the
    /// usable span of `seq`. Respects begin/end region constraints in `settings`.
    pub fn generate_candidates(seq: &DataSequence, settings: &Settings) -> Vec<Primer> {
        let n = seq.len();
        let min_len = settings.min_primer_len;
        let max_len = settings.max_primer_len.min(n);
        if min_len > max_len || n < min_len {
            return Vec::new();
        }

        let region_start = settings.begin_nucleotide.unwrap_or(0).min(n);
        let region_end = settings.end_nucleotide.unwrap_or(n).min(n);
        if region_end <= region_start {
            return Vec::new();
        }

        let mut out = Vec::new();
        for len in min_len..=max_len {
            // Candidate starting positions: anywhere the whole primer fits
            // inside the requested region.
            let last_start = region_end.saturating_sub(len);
            for start in region_start..=last_start {
                out.push(Primer::new(start, len));
            }
        }
        out
    }

    /// Filter candidates through every quality check and attach a penalty.
    pub fn filter_and_score(
        candidates: &[Primer],
        seq: &DataSequence,
        settings: &Settings,
    ) -> Vec<Primer> {
        candidates
            .iter()
            .filter_map(|p| {
                let sub = seq.sub_seq(p.index, p.length);
                if !filters::passes_all(&sub, settings) {
                    return None;
                }
                Some(p.with_penalty(sub.penalty()))
            })
            .collect()
    }

    /// Pair forward and reverse primers respecting amplicon size and Tm-diff
    /// constraints. Sorted ascending by Tm difference (best first).
    pub fn pair_primers(
        fwd_primers: &[Primer],
        rev_primers: &[Primer],
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
        settings: &Settings,
    ) -> Vec<PrimerPair> {
        let template_len = fwd_seq.len();
        let min_amp = settings.min_amplicon_len;
        let max_amp = settings.max_amplicon_len;
        let max_tm_diff = settings.max_tm_diff;

        // Precompute Tms so we don't redo them in the inner loop.
        let fwd_tms: Vec<f32> = fwd_primers
            .iter()
            .map(|p| Thermodynamics::calculate_tm(&fwd_seq.sub_seq(p.index, p.length), settings))
            .collect();
        let rev_tms: Vec<f32> = rev_primers
            .iter()
            .map(|p| Thermodynamics::calculate_tm(&rev_seq.sub_seq(p.index, p.length), settings))
            .collect();

        let mut pairs = Vec::new();
        for (i, fwd) in fwd_primers.iter().enumerate() {
            // Forward primer occupies [fwd.index, fwd.index + fwd.length) on forward strand.
            let fwd_5p = fwd.index;

            for (j, rev) in rev_primers.iter().enumerate() {
                // Reverse primer's index is on the RC strand. Its footprint
                // on the FORWARD strand starts at:
                let rev_5p_on_fwd = template_len
                    .saturating_sub(rev.index)
                    .saturating_sub(rev.length);

                // For a valid PCR pair the reverse must be downstream of forward.
                if rev_5p_on_fwd <= fwd_5p {
                    continue;
                }

                // Amplicon size = distance from forward 5' to reverse 5' end + reverse length.
                let amplicon = rev_5p_on_fwd + rev.length - fwd_5p;
                if amplicon < min_amp {
                    continue;
                }
                if let Some(max) = max_amp {
                    if amplicon > max { continue; }
                }

                let tm_diff = (fwd_tms[i] - rev_tms[j]).abs();
                if tm_diff > max_tm_diff {
                    continue;
                }

                pairs.push(PrimerPair::new(
                    *fwd, *rev, amplicon, tm_diff, fwd_tms[i], rev_tms[j],
                ));
            }
        }

        pairs.sort();
        pairs
    }

    /// Full pipeline: generate, filter, pair, truncate.
    pub fn design_primers(
        seq: &DataSequence,
        rev_seq: &DataSequence,
        settings: &Settings,
    ) -> Result<Vec<PrimerPair>, String> {
        settings.validate()?;

        let fwd_candidates = Self::generate_candidates(seq, settings);
        let rev_candidates = Self::generate_candidates(rev_seq, settings);

        let fwd_good = Self::filter_and_score(&fwd_candidates, seq, settings);
        let rev_good = Self::filter_and_score(&rev_candidates, rev_seq, settings);

        if fwd_good.is_empty() {
            return Err("No forward primers passed filters".into());
        }
        if rev_good.is_empty() {
            return Err("No reverse primers passed filters".into());
        }

        let mut pairs = Self::pair_primers(&fwd_good, &rev_good, seq, rev_seq, settings);
        pairs.truncate(settings.max_primers_output);

        if pairs.is_empty() {
            return Err("Candidate primers exist but no valid pairs (check amplicon size / Tm diff)".into());
        }
        Ok(pairs)
    }

    /// Diagnostic analysis of a user-provided primer string. This is what
    /// the GUI's "test primer" button calls.
    pub fn analyze_primer(seq: &DataSequence, settings: &Settings) -> PrimerAnalysis {
        let tm = Thermodynamics::calculate_tm(seq, settings);
        let filters = filters::apply_all(seq, settings);
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
            filters,
        }
    }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
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
}

impl PrimerAnalysis {
    pub fn all_passed(&self) -> bool {
        self.filters.iter().all(|f| f.passed)
    }

    /// Pretty-print for CLI display.
    pub fn report(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!("Primer:    {}\n", self.sequence));
        s.push_str(&format!("Length:    {} bp  (non-gap: {})\n", self.length, self.actual_size));
        s.push_str(&format!("Tm:        {:.2} °C\n", self.tm));
        s.push_str(&format!("GC:        {:.1}%\n", self.gc_content));
        s.push_str(&format!("Degen pos: {}\n", self.degeneracy_count));
        s.push_str(&format!("Gaps:      {}\n", self.gap_count));
        s.push_str(&format!("Penalty:   {:.2}\n", self.penalty));
        s.push_str(&format!("Hairpin:   {} (matches)\n", self.hairpin_score));
        s.push_str("\nFilters:\n");
        for f in &self.filters {
            let mark = if f.passed { "PASS" } else { "FAIL" };
            s.push_str(&format!("  [{mark}] {}", f.name));
            if let Some(r) = &f.reason {
                s.push_str(&format!(" — {}", r));
            }
            s.push('\n');
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::DataNode;

    fn ds_from(s: &str) -> DataSequence {
        DataSequence::new(s.chars().map(DataNode::from_base).collect())
    }

    #[test]
    fn candidate_count_matches_expectation() {
        // 100 bp, min=20, max=22 ⇒ for each len: (100 - len + 1) starts
        // Total: 81 + 80 + 79 = 240
        let seq = ds_from(&"ATGC".repeat(25));
        let cands = PrimerCalculator::generate_candidates(&seq, &Settings::default());
        assert_eq!(cands.len(), 81 + 80 + 79);
    }

    #[test]
    fn pair_primers_uses_rev_complement_coords() {
        // Construct a 200-bp template; place a fwd primer at the start and
        // a reverse primer at the end on the RC strand. Amplicon should
        // span the full template (200), not some negative or absurd value.
        let template = "A".repeat(200);
        let fwd = ds_from(&template);
        let rev = fwd.reverse_complement();
        let mut settings = Settings::default();
        // Relax everything because all-A primers fail GC etc.
        settings.min_gc_content = 0.0;
        settings.max_gc_content = 1.0;
        settings.max_degeneracy_ratio = 1.0;
        settings.min_tm = 0.0;
        settings.max_tm = 100.0;
        // Skip repeats/comp_ends filters by hand-picking the primers:
        let fwd_primer = Primer::new(0, 20);
        let rev_primer = Primer::new(0, 20);
        // amplicon should equal 200
        let pairs = PrimerCalculator::pair_primers(
            &[fwd_primer], &[rev_primer], &fwd, &rev, &settings,
        );
        // We didn't filter so we got at least the pair.
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].amplicon_size, 200);
    }

    #[test]
    fn pair_rejects_overlap() {
        // Reverse primer on RC strand at a position that would map BEFORE the
        // forward primer on the forward strand → must be skipped.
        let template = "A".repeat(50);
        let fwd = ds_from(&template);
        let rev = fwd.reverse_complement();
        let settings = Settings::default();
        let fwd_primer = Primer::new(20, 20); // forward at position 20
        let rev_primer = Primer::new(20, 20); // RC at 20 → forward strand pos 50-20-20 = 10
        let pairs = PrimerCalculator::pair_primers(
            &[fwd_primer], &[rev_primer], &fwd, &rev, &settings,
        );
        assert!(pairs.is_empty(), "should reject pairs where reverse precedes forward");
    }
}
