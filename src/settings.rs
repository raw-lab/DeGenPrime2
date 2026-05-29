//! Global settings for a DeGenPrime run.
//!
//! Defaults match the original C++ implementation. All setters return `Self`
//! for builder-style chaining and clamp inputs to safe ranges so the
//! algorithm never sees nonsense values.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Settings {
    pub min_primer_len: usize,
    pub max_primer_len: usize,
    pub min_amplicon_len: usize,
    pub max_amplicon_len: Option<usize>,

    pub min_tm: f32,
    pub max_tm: f32,

    pub min_gc_content: f32,   // 0.0 - 1.0
    pub max_gc_content: f32,
    pub max_gc_end: f32,       // GC ratio in last 5 bases, 0.0 - 1.0

    pub max_degeneracy_ratio: f32,
    pub max_deletion_ratio: f32,

    pub primer_concentration: f32, // nM
    pub salt_concentration: f32,   // mM (monovalent Na+)

    pub max_primers_output: usize,
    pub max_tm_diff: f32,

    pub input_file: Option<String>,
    pub output_file: Option<String>,

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
            primer_concentration: 200.0,
            salt_concentration: 50.0,
            max_primers_output: 5,
            max_tm_diff: 5.0,
            input_file: None,
            output_file: None,
            begin_nucleotide: None,
            end_nucleotide: None,
        }
    }
}

impl Settings {
    pub fn new() -> Self { Self::default() }

    // --- builder setters: clamp to safe ranges -----------------------------

    pub fn with_primer_len(mut self, min: usize, max: usize) -> Self {
        self.min_primer_len = min.clamp(18, 25);
        self.max_primer_len = max.clamp(18, 25);
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
        self.min_tm = min.clamp(50.0, 65.0);
        self.max_tm = max.clamp(50.0, 65.0);
        if self.min_tm > self.max_tm {
            std::mem::swap(&mut self.min_tm, &mut self.max_tm);
        }
        self
    }

    /// Set GC content from percentages (0-100).
    pub fn with_gc_content(mut self, min_pct: f32, max_pct: f32) -> Self {
        self.min_gc_content = (min_pct / 100.0).clamp(0.0, 1.0);
        self.max_gc_content = (max_pct / 100.0).clamp(0.0, 1.0);
        if self.min_gc_content > self.max_gc_content {
            std::mem::swap(&mut self.min_gc_content, &mut self.max_gc_content);
        }
        self
    }

    pub fn with_primer_concentration(mut self, c_nm: f32) -> Self {
        // Lower bound 50 nM (per original DeGenPrime contract).
        self.primer_concentration = c_nm.max(50.0);
        self
    }

    pub fn with_salt_concentration(mut self, c_mm: f32) -> Self {
        // Permit lower salt for unusual habitats but keep at least 1 mM.
        self.salt_concentration = c_mm.max(1.0);
        self
    }

    pub fn with_max_primers(mut self, n: usize) -> Self {
        self.max_primers_output = n.clamp(1, 50);
        self
    }

    pub fn with_max_tm_diff(mut self, d: f32) -> Self {
        self.max_tm_diff = d.max(0.0);
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

    pub fn with_region(mut self, begin: usize, end: usize) -> Self {
        let (b, e) = if begin <= end { (begin, end) } else { (end, begin) };
        self.begin_nucleotide = Some(b);
        self.end_nucleotide = Some(e);
        self
    }

    /// Final sanity check. Most invariants are enforced by the setters; this
    /// catches the few that cross-cut multiple fields.
    pub fn validate(&self) -> Result<(), String> {
        if self.min_primer_len >= self.max_primer_len + 1 && self.min_primer_len > self.max_primer_len {
            return Err("min_primer_len must be <= max_primer_len".into());
        }
        if self.min_tm >= self.max_tm {
            return Err("min_tm must be < max_tm".into());
        }
        if self.min_gc_content >= self.max_gc_content {
            return Err("min_gc must be < max_gc".into());
        }
        if let Some(max_amp) = self.max_amplicon_len {
            if max_amp < self.min_amplicon_len {
                return Err("max_amplicon_len must be >= min_amplicon_len".into());
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
        // both clamped to [50, 65], then min/max swap so min < max
        assert!(s.min_tm < s.max_tm);
        assert!(s.min_tm >= 50.0 && s.max_tm <= 65.0);
    }

    #[test]
    fn gc_setter_converts_pct() {
        let s = Settings::default().with_gc_content(45.0, 55.0);
        assert!((s.min_gc_content - 0.45).abs() < 1e-6);
        assert!((s.max_gc_content - 0.55).abs() < 1e-6);
    }
}
