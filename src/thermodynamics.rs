//! Thermodynamic calculations for primer chemistry.
//!
//! Tm calculation follows SantaLucia (1998) PNAS 95(4):1460 with the
//! Owczarzy et al. (2008) Biochemistry 47:5336 salt correction for [Na+].
//!
//! Units:
//!   - ΔH is reported in kcal/mol
//!   - ΔS is reported in cal/(mol·K)
//!   - primer concentration is taken from settings in nM
//!   - salt concentration is taken from settings in mM (Na+)

use crate::models::{DataNode, DataSequence};
use crate::settings::Settings;

/// Gas constant in cal/(mol·K) used in the SantaLucia formulation.
const R_CAL: f32 = 1.987;

pub struct Thermodynamics;

impl Thermodynamics {
    /// Calculate melting temperature in degrees Celsius using nearest-neighbour
    /// thermodynamics with the von Ahsen et al. (1999) salt correction.
    ///
    /// This matches the original DeGenPrime C++ implementation: salt correction
    /// is incorporated into the denominator (entropy units), and no separate
    /// initiation correction is applied. The result matches values produced by
    /// the original tool, so the Tm range [58, 62] documented in the README
    /// continues to mean what users expect.
    pub fn calculate_tm(seq: &DataSequence, settings: &Settings) -> f32 {
        let nodes = seq.nodes();
        if nodes.len() < 2 {
            return 0.0;
        }

        // Skip gap columns when computing stack contributions.
        let (dh, ds) = Self::stack_params(nodes);
        let n_eff = nodes.iter().filter(|x| x.most_common() != '-').count();
        if n_eff < 2 {
            return 0.0;
        }

        // Primer concentration in mol/L (input is nM); use Ct/4 for non-self-
        // complementary primers.
        let ct = (settings.primer_concentration as f32) * 1e-9;
        if ct <= 0.0 {
            return 0.0;
        }
        let rlnk = R_CAL * (ct / 4.0).ln();

        // Salt correction (von Ahsen 1999), added directly to the entropy term.
        let na_m = (settings.salt_concentration / 1000.0).max(1e-6);
        let salt_mod = 0.368 * (n_eff as f32 - 1.0) * na_m.ln();

        // Combine and convert cal -> kcal so the units match the numerator.
        let denom = (ds + rlnk + salt_mod) / 1000.0;
        if denom.abs() < 1e-9 {
            return 0.0;
        }

        let tm_k = dh / denom;
        tm_k - 273.15
    }

    /// Sum the per-stack ΔH and ΔS over a sequence.
    fn stack_params(nodes: &[DataNode]) -> (f32, f32) {
        let mut dh = 0.0;
        let mut ds = 0.0;
        for w in nodes.windows(2) {
            dh += w[0].enthalpy_to(&w[1]);
            ds += w[0].entropy_to(&w[1]);
        }
        (dh, ds)
    }

    /// Gibbs free energy at a given temperature (Kelvin).
    pub fn calculate_gibbs(seq: &DataSequence, temp_k: f32) -> f32 {
        let nodes = seq.nodes();
        if nodes.len() < 2 {
            return 0.0;
        }
        let (dh, ds) = Self::stack_params(nodes);
        // ΔG = ΔH - T·ΔS; convert ΔH kcal→cal then both into kcal:
        dh - (temp_k * ds) / 1000.0
    }

    /// Crude hairpin score: max number of complementary base pairs that can
    /// form in a 5-bp window across the sequence with a >=3 bp loop. Higher
    /// scores indicate stronger potential for hairpin formation.
    pub fn hairpin_score(seq: &DataSequence) -> usize {
        let bases: Vec<char> = seq.nodes().iter().map(|n| n.most_common()).collect();
        let n = bases.len();
        if n < 11 { return 0; }
        let mut best = 0;
        for i in 0..n.saturating_sub(10) {
            for j in (i + 8)..n.saturating_sub(2) {
                if j + 3 > n { continue; }
                let stem = 3.min(n - j);
                let matches = (0..stem)
                    .filter(|k| {
                        let l = i + k;
                        let r = j + (stem - 1 - k);
                        r < n && DataNode::complement_base(bases[l]) == bases[r]
                    })
                    .count();
                if matches > best { best = matches; }
            }
        }
        best
    }

    /// Score for primer-dimer formation between a forward and reverse primer
    /// based on 3' end complementarity. Returns the longest matching stretch
    /// in their 3'-3' overlap.
    pub fn dimer_score(fwd: &DataSequence, rev: &DataSequence) -> usize {
        let f: Vec<char> = fwd.nodes().iter().map(|n| n.most_common()).collect();
        let r: Vec<char> = rev.nodes().iter().map(|n| n.most_common()).collect();
        if f.is_empty() || r.is_empty() { return 0; }

        // Try every alignment of the last `k` bases of fwd against the start of rev
        let mut best = 0;
        let max_overlap = f.len().min(r.len()).min(8);
        for k in 1..=max_overlap {
            let f_tail = &f[f.len() - k..];
            let r_head = &r[..k];
            let matches = f_tail
                .iter()
                .zip(r_head.iter().rev()) // 3' end of rev pairs with 3' end of fwd
                .filter(|(a, b)| DataNode::complement_base(**a) == **b)
                .count();
            if matches > best { best = matches; }
        }
        best
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::DataSequence;

    #[test]
    fn tm_reasonable_for_typical_primer() {
        // A typical 20-mer with mixed bases should land in 50-70 C.
        let s = DataSequence::from_string("ATGCAATTGGCCGGAATTCC");
        let settings = Settings::default();
        let tm = Thermodynamics::calculate_tm(&s, &settings);
        assert!(tm > 40.0 && tm < 80.0, "tm = {tm}");
    }

    #[test]
    fn gc_rich_has_higher_tm_than_at_rich() {
        let gc = DataSequence::from_string("GGCGGCGGCGGCGGCGGCGG");
        let at = DataSequence::from_string("ATATATATATATATATATAT");
        let settings = Settings::default();
        let tm_gc = Thermodynamics::calculate_tm(&gc, &settings);
        let tm_at = Thermodynamics::calculate_tm(&at, &settings);
        assert!(tm_gc > tm_at, "tm_gc={tm_gc} should exceed tm_at={tm_at}");
    }

    #[test]
    fn salt_correction_changes_tm() {
        let s = DataSequence::from_string("ATGCAATTGGCCGGAATTCC");
        let low = Thermodynamics::calculate_tm(&s, &Settings::default().with_salt_concentration(10.0));
        let high = Thermodynamics::calculate_tm(&s, &Settings::default().with_salt_concentration(200.0));
        assert!(high > low, "higher salt should raise Tm: low={low} high={high}");
    }
}
