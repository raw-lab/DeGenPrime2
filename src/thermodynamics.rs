use crate::models::{DataSequence, DataNode};
use crate::settings::Settings;

/// Thermodynamic calculations for primer melting temperature and secondary structure prediction
pub struct Thermodynamics;

impl Thermodynamics {
    /// Calculate melting temperature using nearest-neighbor thermodynamics
    /// Based on SantaLucia (1998) and Owczarzy et al. (2008)
    pub fn calculate_tm(seq: &DataSequence, settings: &Settings) -> f32 {
        let nodes = seq.nodes();

        if nodes.is_empty() {
            return 0.0;
        }

        // Calculate ΔH and ΔS using nearest-neighbor parameters
        let (delta_h, delta_s) = Self::calculate_thermodynamic_params(nodes);

        // Apply salt correction
        let na_corr = Self::salt_correction(settings.salt_concentration());

        // Calculate Tm = ΔH / (ΔS + R*ln(primer_conc/2)) + salt_correction
        let r = 1.987; // cal/(mol·K)
        let primer_conc = settings.primer_concentration() / 1_000_000_000.0; // Convert nM to M

        let denominator = delta_s + r * (primer_conc / 2.0).ln();

        if denominator.abs() < 1e-10 {
            return 0.0;
        }

        let tm = (delta_h * 1000.0 / denominator) - 273.15 + na_corr;
        tm.max(0.0)
    }

    /// Calculate ΔH and ΔS for a sequence using nearest-neighbor parameters
    fn calculate_thermodynamic_params(nodes: &[DataNode]) -> (f32, f32) {
        // Initial values for DNA-DNA hybrid at 37°C
        let mut delta_h = -3.1; // kcal/mol - initial (terminal AT penalty)
        let mut delta_s = -12.6; // e.u. (entropy units)

        if nodes.len() > 1 {
            // Add contribution from each nearest-neighbor pair
            for window in nodes.windows(2) {
                delta_h += window[0].enthalpy_to(&window[1]);
                delta_s += window[0].entropy_to(&window[1]);
            }
        } else if nodes.len() == 1 {
            // Single base - use self-interaction
            delta_h += nodes[0].enthalpy_to(&nodes[0]);
            delta_s += nodes[0].entropy_to(&nodes[0]);
        }

        (delta_h, delta_s)
    }

    /// Salt correction based on Owczarzy et al. (2008)
    fn salt_correction(na_conc_mm: f32) -> f32 {
        if na_conc_mm <= 0.0 {
            return 0.0;
        }

        // Convert mM to M
        let na_m = na_conc_mm / 1000.0;

        // Correction formula: 12.0 * ln([Na+])
        // Works well for [Na+] between 0.01 and 1.0 M
        if na_m < 0.01 || na_m > 1.0 {
            // Use linear approximation for out-of-range values
            16.6 * (na_m.log10())
        } else {
            12.0 * na_m.ln()
        }
    }

    /// Calculate Gibbs free energy for sequence stability
    pub fn calculate_gibbs(seq: &DataSequence, temp_k: f32) -> f32 {
        let nodes = seq.nodes();

        if nodes.is_empty() {
            return 0.0;
        }

        let (delta_h, delta_s) = Self::calculate_thermodynamic_params(nodes);

        // ΔG = ΔH - T*ΔS
        // Convert units: ΔH in kcal/mol, ΔS in cal/(mol·K)
        let delta_g = (delta_h * 1000.0) - (temp_k * delta_s);

        delta_g
    }

    /// Check for potential hairpin (intramolecular) secondary structures
    pub fn hairpin_check(seq: &DataSequence) -> Option<(usize, usize)> {
        let codes = seq.codes();
        let n = codes.len();

        // Check for complementary regions that could fold back
        // Look for GC-rich complementary sequences separated by loop region
        for loop_size in 3..=10 {
            for i in 0..(n.saturating_sub(2 * 5 + loop_size)) {
                let stem_size = 5; // Minimum stem length
                
                if i + stem_size + loop_size + stem_size > n {
                    continue;
                }

                let start_seq: String = codes.chars().skip(i).take(stem_size).collect();
                let loop_end = i + stem_size + loop_size;
                let end_seq: String = codes
                    .chars()
                    .skip(loop_end)
                    .take(stem_size)
                    .collect();

                if Self::is_complement_pair(&start_seq, &end_seq) {
                    return Some((i, loop_end));
                }
            }
        }

        None
    }

    /// Check for dimer formation potential
    pub fn dimer_check(fwd: &DataSequence, rev: &DataSequence) -> f32 {
        let fwd_codes = fwd.codes();
        let rev_codes = rev.codes();

        // Check 3' end complementarity (3 bp minimum)
        let mut max_match = 0;
        let min_match = 3;

        for fwd_end in 1..=fwd_codes.len().min(10) {
            for rev_start in 0..=rev_codes.len().saturating_sub(min_match) {
                let fwd_slice: String = fwd_codes
                    .chars()
                    .rev()
                    .take(fwd_end)
                    .collect::<Vec<_>>()
                    .iter()
                    .rev()
                    .collect();

                let rev_slice: String = rev_codes.chars().skip(rev_start).take(fwd_end).collect();

                let matches = fwd_slice
                    .chars()
                    .zip(rev_slice.chars())
                    .filter(|(a, b)| Self::is_complement(*a, *b))
                    .count();

                if matches >= min_match {
                    max_match = max_match.max(matches);
                }
            }
        }

        max_match as f32
    }

    /// Check if two nucleotides are complementary
    fn is_complement(base1: char, base2: char) -> bool {
        matches!(
            (base1.to_ascii_uppercase(), base2.to_ascii_uppercase()),
            ('A', 'T') | ('T', 'A') | ('G', 'C') | ('C', 'G')
        )
    }

    /// Check if two sequences are reverse complements
    fn is_complement_pair(seq1: &str, seq2: &str) -> bool {
        if seq1.len() != seq2.len() {
            return false;
        }

        seq1.chars()
            .zip(seq2.chars().rev())
            .all(|(a, b)| Self::is_complement(a, b))
    }

    /// Calculate annealing temperature based on Tm
    pub fn calculate_annealing_temp(tm: f32) -> f32 {
        // Typical PCR annealing temperature is 3-5°C below the Tm
        (tm - 5.0).max(4.0)
    }

    /// Estimate extension time based on amplicon length
    pub fn estimate_extension_time(amplicon_len: usize) -> f32 {
        // Assume ~1000 bp per second for Taq polymerase
        let seconds = (amplicon_len as f32) / 1000.0;
        seconds.max(0.5) // Minimum 0.5 seconds
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::{DataNode, DataSequence};

    #[test]
    fn test_salt_correction() {
        // Salt correction should increase with salt concentration
        let corr_low = Thermodynamics::salt_correction(50.0);
        let corr_high = Thermodynamics::salt_correction(100.0);
        assert!(corr_high > corr_low);
    }

    #[test]
    fn test_calculate_tm() {
        let nodes = vec![
            DataNode::new(vec!['A']),
            DataNode::new(vec!['T']),
            DataNode::new(vec!['G']),
            DataNode::new(vec!['C']),
        ];
        let seq = DataSequence::new(nodes);
        let settings = Settings::default();

        let tm = Thermodynamics::calculate_tm(&seq, &settings);
        assert!(tm > 0.0 && tm < 100.0);
    }

    #[test]
    fn test_dimer_check() {
        let fwd_nodes = vec![
            DataNode::new(vec!['A']),
            DataNode::new(vec!['T']),
            DataNode::new(vec!['G']),
        ];
        let rev_nodes = vec![
            DataNode::new(vec!['C']),
            DataNode::new(vec!['A']),
            DataNode::new(vec!['T']),
        ];

        let fwd = DataSequence::new(fwd_nodes);
        let rev = DataSequence::new(rev_nodes);

        let score = Thermodynamics::dimer_check(&fwd, &rev);
        assert!(score >= 0.0);
    }
}
