//! Settings for a DeGenPrime run.
//!
//! Defaults match the original C++ implementation where one existed; the
//! Primer3-style complementarity limits use Primer3's own defaults. Builder
//! setters clamp to the documented safe ranges; callers that want to tell the
//! user about a clamp should go through [`crate::api::SettingsInput`], which
//! reports every adjustment as a warning instead of applying it silently.

use serde::{Deserialize, Serialize};

/// Documented safe ranges (shared by the CLI help, the web form and the API).
pub mod limits {
    pub const PRIMER_LEN: (usize, usize) = (18, 25);
    pub const TM: (f32, f32) = (50.0, 65.0);
    pub const MAX_PRIMERS: (usize, usize) = (1, 50);
    pub const MIN_PRIMER_CONC_NM: f32 = 50.0;
    pub const MIN_SALT_MM: f32 = 1.0;
}

/// How primer pairs are ranked.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum RankBy {
    /// Conservation penalties + Tm terms (see `calculator` docs). Default.
    #[default]
    Score,
    /// 2.1 behaviour: smallest ΔTm first, conservation as tie-break.
    TmDiff,
}

impl RankBy {
    pub fn parse(s: &str) -> Option<RankBy> {
        match s.trim().to_ascii_lowercase().replace('-', "_").as_str() {
            "score" | "penalty" | "combined" => Some(RankBy::Score),
            "tm_diff" | "tmdiff" | "tm" | "legacy" => Some(RankBy::TmDiff),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct Settings {
    pub min_primer_len: usize,
    pub max_primer_len: usize,
    pub min_amplicon_len: usize,
    pub max_amplicon_len: Option<usize>,

    pub min_tm: f32,
    pub max_tm: f32,

    /// GC fractions (0.0 – 1.0).
    pub min_gc_content: f32,
    pub max_gc_content: f32,
    /// Max GC fraction in the last five 3' bases.
    pub max_gc_end: f32,

    /// Max fraction of degenerate positions in a primer.
    pub max_degeneracy_ratio: f32,
    /// Max fraction of gap (deletion) columns in a primer window.
    pub max_deletion_ratio: f32,
    /// Longest allowed mononucleotide run (Primer3's PRIMER_MAX_POLY_X).
    pub max_poly_x: usize,

    /// Primer3-style complementarity limits (alignment scores, 1 per base pair).
    pub max_self_any: f32,
    pub max_self_end: f32,
    pub max_pair_any: f32,
    pub max_pair_end: f32,

    pub primer_concentration: f32, // nM
    pub salt_concentration: f32,   // mM monovalent cation

    pub max_primers_output: usize,
    pub max_tm_diff: f32,
    /// Returned pairs must differ by at least this many bases at the 3' end
    /// of the forward or the reverse primer. 0 allows near-duplicates.
    pub min_three_prime_distance: usize,
    pub rank_by: RankBy,

    pub input_file: Option<String>,
    pub output_file: Option<String>,

    /// Optional design region in alignment columns: `[begin, end)`, 0-based.
    pub begin_nucleotide: Option<usize>,
    pub end_nucleotide: Option<usize>,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            min_primer_len: 20,
            max_primer_len: 22,
            min_amplicon_len: 100,
            max_amplicon_len: None,
            min_tm: 58.0,
            max_tm: 62.0,
            min_gc_content: 0.40,
            max_gc_content: 0.60,
            max_gc_end: 0.80,
            max_degeneracy_ratio: 0.15,
            max_deletion_ratio: 0.15,
            max_poly_x: 4,
            max_self_any: 8.0,
            max_self_end: 3.0,
            max_pair_any: 8.0,
            max_pair_end: 3.0,
            primer_concentration: 200.0,
            salt_concentration: 50.0,
            max_primers_output: 5,
            max_tm_diff: 5.0,
            min_three_prime_distance: 3,
            rank_by: RankBy::Score,
            input_file: None,
            output_file: None,
            begin_nucleotide: None,
            end_nucleotide: None,
        }
    }
}

impl Settings {
    pub fn new() -> Self {
        Self::default()
    }

    /// Optimal Tm used by the ranking score: the middle of the Tm window.
    pub fn tm_optimum(&self) -> f32 {
        (self.min_tm + self.max_tm) / 2.0
    }

    // --- builder setters: clamp to safe ranges -----------------------------

    pub fn with_primer_len(mut self, min: usize, max: usize) -> Self {
        let (lo, hi) = limits::PRIMER_LEN;
        self.min_primer_len = min.clamp(lo, hi);
        self.max_primer_len = max.clamp(lo, hi);
        if self.min_primer_len > self.max_primer_len {
            std::mem::swap(&mut self.min_primer_len, &mut self.max_primer_len);
        }
        self
    }

    pub fn with_amplicon_len(mut self, min: usize) -> Self {
        self.min_amplicon_len = min;
        self
    }

    pub fn with_max_amplicon_len(mut self, max: usize) -> Self {
        self.max_amplicon_len = Some(max);
        self
    }

    pub fn with_tm_range(mut self, min: f32, max: f32) -> Self {
        let (lo, hi) = limits::TM;
        self.min_tm = min.clamp(lo, hi);
        self.max_tm = max.clamp(lo, hi);
        if self.min_tm > self.max_tm {
            std::mem::swap(&mut self.min_tm, &mut self.max_tm);
        }
        self
    }

    /// Set GC content from percentages (0–100).
    pub fn with_gc_content(mut self, min_pct: f32, max_pct: f32) -> Self {
        self.min_gc_content = (min_pct / 100.0).clamp(0.0, 1.0);
        self.max_gc_content = (max_pct / 100.0).clamp(0.0, 1.0);
        if self.min_gc_content > self.max_gc_content {
            std::mem::swap(&mut self.min_gc_content, &mut self.max_gc_content);
        }
        self
    }

    /// Max GC in the last five 3' bases, as a percentage.
    pub fn with_max_gc_end(mut self, pct: f32) -> Self {
        self.max_gc_end = (pct / 100.0).clamp(0.0, 1.0);
        self
    }

    pub fn with_max_degeneracy(mut self, pct: f32) -> Self {
        self.max_degeneracy_ratio = (pct / 100.0).clamp(0.0, 1.0);
        self
    }

    pub fn with_max_deletions(mut self, pct: f32) -> Self {
        self.max_deletion_ratio = (pct / 100.0).clamp(0.0, 1.0);
        self
    }

    pub fn with_max_poly_x(mut self, n: usize) -> Self {
        self.max_poly_x = n.max(1);
        self
    }

    pub fn with_complementarity(mut self, self_any: f32, self_end: f32, pair_any: f32, pair_end: f32) -> Self {
        self.max_self_any = self_any.max(0.0);
        self.max_self_end = self_end.max(0.0);
        self.max_pair_any = pair_any.max(0.0);
        self.max_pair_end = pair_end.max(0.0);
        self
    }

    pub fn with_primer_concentration(mut self, c_nm: f32) -> Self {
        self.primer_concentration = c_nm.max(limits::MIN_PRIMER_CONC_NM);
        self
    }

    pub fn with_salt_concentration(mut self, c_mm: f32) -> Self {
        self.salt_concentration = c_mm.max(limits::MIN_SALT_MM);
        self
    }

    pub fn with_max_primers(mut self, n: usize) -> Self {
        let (lo, hi) = limits::MAX_PRIMERS;
        self.max_primers_output = n.clamp(lo, hi);
        self
    }

    pub fn with_max_tm_diff(mut self, d: f32) -> Self {
        self.max_tm_diff = d.max(0.0);
        self
    }

    pub fn with_min_three_prime_distance(mut self, d: usize) -> Self {
        self.min_three_prime_distance = d;
        self
    }

    pub fn with_rank_by(mut self, r: RankBy) -> Self {
        self.rank_by = r;
        self
    }

    pub fn with_input_file(mut self, p: String) -> Self {
        self.input_file = Some(p);
        self
    }

    pub fn with_output_file(mut self, p: String) -> Self {
        self.output_file = Some(p);
        self
    }

    /// Restrict design to alignment columns `[begin, end)` (0-based).
    pub fn with_region(mut self, begin: usize, end: usize) -> Self {
        let (b, e) = if begin <= end { (begin, end) } else { (end, begin) };
        self.begin_nucleotide = Some(b);
        self.end_nucleotide = Some(e);
        self
    }

    /// Set either end of the region independently (2.1 ignored `--begin`
    /// unless `--end` was also given).
    pub fn with_region_bounds(mut self, begin: Option<usize>, end: Option<usize>) -> Self {
        match (begin, end) {
            (Some(b), Some(e)) => return self.with_region(b, e),
            (b, e) => {
                self.begin_nucleotide = b;
                self.end_nucleotide = e;
            }
        }
        self
    }

    /// Cross-field sanity checks.
    pub fn validate(&self) -> Result<(), String> {
        let numbers = [
            ("Tm minimum", self.min_tm),
            ("Tm maximum", self.max_tm),
            ("GC minimum", self.min_gc_content),
            ("GC maximum", self.max_gc_content),
            ("3' end GC", self.max_gc_end),
            ("degenerate positions", self.max_degeneracy_ratio),
            ("gap columns", self.max_deletion_ratio),
            ("self-complementarity", self.max_self_any),
            ("3' self-complementarity", self.max_self_end),
            ("pair complementarity", self.max_pair_any),
            ("3' pair complementarity", self.max_pair_end),
            ("primer concentration", self.primer_concentration),
            ("salt concentration", self.salt_concentration),
            ("maximum ΔTm", self.max_tm_diff),
        ];
        if let Some((name, _)) = numbers.iter().find(|(_, v)| !v.is_finite()) {
            return Err(format!("The {name} must be a finite number."));
        }
        if self.min_primer_len == 0 || self.min_primer_len > self.max_primer_len {
            return Err("min_primer_len must be between 1 and max_primer_len".into());
        }
        if self.min_tm >= self.max_tm {
            return Err("min_tm must be lower than max_tm".into());
        }
        if self.min_gc_content >= self.max_gc_content {
            return Err("min_gc must be lower than max_gc".into());
        }
        if let Some(max_amp) = self.max_amplicon_len {
            if max_amp < self.min_amplicon_len {
                return Err("max_amplicon must be at least the minimum amplicon size".into());
            }
        }
        if let (Some(b), Some(e)) = (self.begin_nucleotide, self.end_nucleotide) {
            if e <= b {
                return Err("the design region is empty (end must be after begin)".into());
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_validate() {
        assert!(Settings::default().validate().is_ok());
    }

    #[test]
    fn tm_setter_clamps_and_swaps() {
        let s = Settings::default().with_tm_range(70.0, 40.0);
        assert!(s.min_tm < s.max_tm);
        assert!(s.min_tm >= 50.0 && s.max_tm <= 65.0);
    }

    #[test]
    fn gc_setter_converts_pct() {
        let s = Settings::default().with_gc_content(45.0, 55.0);
        assert!((s.min_gc_content - 0.45).abs() < 1e-6);
        assert!((s.max_gc_content - 0.55).abs() < 1e-6);
    }

    #[test]
    fn region_bounds_can_be_one_sided() {
        let s = Settings::default().with_region_bounds(Some(900), None);
        assert_eq!(s.begin_nucleotide, Some(900));
        assert_eq!(s.end_nucleotide, None);
    }

    #[test]
    fn old_json_without_new_fields_still_loads() {
        let old = r#"{"min_primer_len":20,"max_primer_len":22,"min_amplicon_len":100,
            "max_amplicon_len":null,"min_tm":58.0,"max_tm":62.0,"min_gc_content":0.4,
            "max_gc_content":0.6,"max_gc_end":0.8,"max_degeneracy_ratio":0.15,
            "max_deletion_ratio":0.15,"primer_concentration":200.0,"salt_concentration":50.0,
            "max_primers_output":5,"max_tm_diff":5.0,"input_file":null,"output_file":null,
            "begin_nucleotide":null,"end_nucleotide":null}"#;
        let s: Settings = serde_json::from_str(old).unwrap();
        assert_eq!(s.max_poly_x, 4);
        assert_eq!(s.rank_by, RankBy::Score);
    }
}
