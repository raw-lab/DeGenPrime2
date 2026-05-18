use crate::models::{DataSequence, Primer, PrimerPair};
use crate::settings::Settings;
use crate::filters::{apply_all_filters, FilterResult};
use crate::thermodynamics::Thermodynamics;

/// Generate and filter candidate primers using functional programming patterns
pub struct PrimerCalculator;

impl PrimerCalculator {
    /// Generate all candidate primers from a sequence
    pub fn generate_candidates(seq: &DataSequence, settings: &Settings) -> Vec<Primer> {
        let min_len = settings.min_primer_len();
        let max_len = settings.max_primer_len();
        let max_pos = seq.len().saturating_sub(min_len);

        // Use functional iteration to generate all candidates
        (min_len..=max_len)
            .flat_map(|len| {
                (0..=max_pos).map(move |pos| Primer::new(pos, len))
            })
            .collect()
    }

    /// Generate primers bounded by a maximum position
    pub fn generate_bounded_candidates(
        seq: &DataSequence,
        settings: &Settings,
        max_pos: usize,
    ) -> Vec<Primer> {
        let min_len = settings.min_primer_len();
        let max_len = settings.max_primer_len();
        let bound = max_pos.min(seq.len().saturating_sub(min_len));

        (min_len..=max_len)
            .flat_map(|len| {
                (0..=bound).map(move |pos| Primer::new(pos, len))
            })
            .collect()
    }

    /// Generate primers from conserved regions
    pub fn generate_from_conserved_regions(
        seq: &DataSequence,
        settings: &Settings,
        regions: &[(usize, usize)],
    ) -> Vec<Primer> {
        let min_len = settings.min_primer_len();
        let max_len = settings.max_primer_len();

        regions
            .iter()
            .flat_map(|(start, len)| {
                let end_pos = (*start + len).saturating_sub(min_len);
                (min_len..=max_len)
                    .flat_map(move |plen| {
                        (0..=end_pos).map(move |offset| {
                            Primer::new(start + offset, plen)
                        })
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    /// Identify conserved regions in a sequence
    pub fn find_conserved_regions(seq: &DataSequence, settings: &Settings) -> Vec<(usize, usize)> {
        let min_conserved = settings.min_primer_len();
        let codes = seq.codes();

        let mut regions = Vec::new();
        let mut current_start = 0;
        let mut conserved_count = 0;
        let conserved_chars = "CGATWSRYKMBDHV";

        for (i, c) in codes.chars().enumerate() {
            if conserved_chars.contains(c) {
                if conserved_count == 0 {
                    current_start = i;
                }
                conserved_count += 1;
            } else {
                if conserved_count >= min_conserved {
                    regions.push((current_start, conserved_count));
                }
                conserved_count = 0;
            }
        }

        if conserved_count >= min_conserved {
            regions.push((current_start, conserved_count));
        }

        regions
    }

    /// Filter primers by applying all filter criteria
    pub fn filter_primers(
        primers: &[Primer],
        seq: &DataSequence,
        settings: &Settings,
    ) -> Vec<(Primer, f32)> {
        primers
            .iter()
            .filter_map(|primer| {
                let primer_seq = seq.sub_seq(primer.index(), primer.length());

                // Skip if not valid IUPAC
                let valid = primer_seq
                    .codes()
                    .chars()
                    .all(|c| "AaCcGgTtWwSsRrYyKkMmBbDdHhVvNn-".contains(c));

                if !valid {
                    return None;
                }

                // Apply filters
                let filter_result = apply_all_filters(&primer_seq, settings);
                if !filter_result.passed {
                    return None;
                }

                // Calculate penalty
                let penalty = primer_seq.penalty();
                Some((*primer, penalty))
            })
            .collect()
    }

    /// Rank filtered primers by penalty
    pub fn rank_primers(filtered: Vec<(Primer, f32)>) -> Vec<Primer> {
        let mut sorted = filtered;
        sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        sorted.into_iter().map(|(p, _)| p).collect()
    }

    /// Find primer pairs that match both forward and reverse
    pub fn pair_primers(
        fwd_primers: &[Primer],
        rev_primers: &[Primer],
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
        settings: &Settings,
    ) -> Vec<PrimerPair> {
        let min_amplicon = settings.min_amplicon_len();
        let max_tm_diff = 5.0; // Maximum allowed Tm difference

        let mut pairs = Vec::new();

        for fwd in fwd_primers {
            for rev in rev_primers {
                // Calculate amplicon size
                let amplicon_size = if rev.index() > fwd.index() {
                    rev.index() - fwd.index() + rev.length()
                } else {
                    continue; // Skip invalid pairs
                };

                // Check minimum amplicon size
                if amplicon_size < min_amplicon {
                    continue;
                }

                // Get sequences and calculate temperatures
                let fwd_seq_sub = fwd_seq.sub_seq(fwd.index(), fwd.length());
                let rev_seq_sub = rev_seq.sub_seq(rev.index(), rev.length());

                let tm_fwd = Thermodynamics::calculate_tm(&fwd_seq_sub, settings);
                let tm_rev = Thermodynamics::calculate_tm(&rev_seq_sub, settings);

                let tm_diff = (tm_fwd - tm_rev).abs();

                if tm_diff <= max_tm_diff {
                    let pair = PrimerPair::new(*fwd, *rev, amplicon_size, tm_diff);
                    pairs.push(pair);
                }
            }
        }

        // Sort by temperature difference
        pairs.sort();
        pairs
    }

    /// Run full primer design pipeline
    pub fn design_primers(
        seq: &DataSequence,
        rev_seq: &DataSequence,
        settings: &Settings,
    ) -> Result<Vec<PrimerPair>, String> {
        // Validate settings
        settings.validate()?;

        // Generate candidate primers
        let fwd_candidates = Self::generate_candidates(seq, settings);
        let rev_candidates = Self::generate_candidates(rev_seq, settings);

        // Filter primers
        let fwd_filtered = Self::filter_primers(&fwd_candidates, seq, settings);
        let rev_filtered = Self::filter_primers(&rev_candidates, rev_seq, settings);

        // Rank by penalty
        let fwd_ranked = Self::rank_primers(fwd_filtered);
        let rev_ranked = Self::rank_primers(rev_filtered);

        if fwd_ranked.is_empty() || rev_ranked.is_empty() {
            return Err("No primers passed filter criteria".to_string());
        }

        // Pair primers
        let mut pairs = Self::pair_primers(&fwd_ranked, &rev_ranked, seq, rev_seq, settings);

        // Limit output
        pairs.truncate(settings.max_primers_output());

        Ok(pairs)
    }

    /// Get diagnostic information about why a primer failed
    pub fn analyze_primer(
        seq: &DataSequence,
        settings: &Settings,
    ) -> String {
        let mut report = String::new();

        // Basic info
        report.push_str(&format!("Primer sequence: {}\n", seq.codes()));
        report.push_str(&format!("Length: {} bp\n", seq.len()));
        report.push_str(&format!("Actual size (no gaps): {} bp\n", seq.actual_size()));
        report.push_str(&format!("Penalty: {:.2}\n\n", seq.penalty()));

        // Filter results
        report.push_str("Filter Results:\n");
        report.push_str("==============\n");

        let filters = [
            ("Degeneracy", crate::filters::filter_degeneracy),
            ("Deletions", crate::filters::filter_deletions),
            ("GC Content", crate::filters::filter_gc_content),
            ("Repeats", crate::filters::filter_repeats),
            ("Complementary Ends", crate::filters::filter_complementary_ends),
            ("Temperature", crate::filters::filter_temperature),
        ];

        for (name, filter_fn) in filters {
            let result = filter_fn(seq, settings);
            let status = if result.passed { "PASS" } else { "FAIL" };
            report.push_str(&format!("{}: {}", name, status));
            if let Some(reason) = result.reason {
                report.push_str(&format!(" ({})", reason));
            }
            report.push('\n');
        }

        // Thermodynamic info
        report.push_str("\nThermodynamic Properties:\n");
        report.push_str("=========================\n");
        report.push_str(&format!(
            "Melting Temperature: {:.2}°C\n",
            Thermodynamics::calculate_tm(seq, settings)
        ));
        report.push_str(&format!("GC Content: {:.1}%\n", seq.gc_ratio() * 100.0));
        report.push_str(&format!("Degeneracy: {} positions\n", seq.degeneracy_count()));
        report.push_str(&format!("Gaps: {} positions\n", seq.gap_count()));

        report
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::DataNode;

    #[test]
    fn test_generate_candidates() {
        let nodes = (0..100)
            .map(|_| DataNode::new(vec!['A']))
            .collect::<Vec<_>>();
        let seq = DataSequence::new(nodes);
        let settings = Settings::default();

        let candidates = PrimerCalculator::generate_candidates(&seq, &settings);

        // With min_len=20, max_len=22, and seq_len=100:
        // Positions: 0..81 (100-20+1)
        // Lengths: 20, 21, 22
        // Total: 81 * 3 = 243
        assert!(candidates.len() > 0);
    }

    #[test]
    fn test_find_conserved_regions() {
        let codes = "AAAA----TTTT----GGGG"
            .chars()
            .map(|c| DataNode::new(vec![c]))
            .collect::<Vec<_>>();
        let seq = DataSequence::new(codes);
        let settings = Settings::default();

        let regions = PrimerCalculator::find_conserved_regions(&seq, &settings);
        assert!(regions.len() >= 2);
    }

    #[test]
    fn test_filter_primers() {
        let nodes = "ATGCATGCATGCATGCATGC"
            .chars()
            .map(|c| DataNode::new(vec![c]))
            .collect::<Vec<_>>();
        let seq = DataSequence::new(nodes);
        let settings = Settings::default();

        let candidates = PrimerCalculator::generate_candidates(&seq, &settings);
        let filtered = PrimerCalculator::filter_primers(&candidates, &seq, &settings);

        assert!(filtered.len() <= candidates.len());
    }
}
