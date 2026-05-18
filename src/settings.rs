/// Global settings for DeGenPrime
#[derive(Debug, Clone)]
pub struct Settings {
    // Primer length constraints
    min_primer_len: usize,
    max_primer_len: usize,

    // Amplicon constraints
    min_amplicon_len: usize,

    // Temperature constraints (°C)
    min_tm: f32,
    max_tm: f32,

    // GC content constraints
    min_gc_content: f32,
    max_gc_content: f32,
    max_gc_end: f32,

    // Filter tolerances
    max_degeneracy_ratio: f32,
    max_deletion_ratio: f32,

    // Thermodynamic parameters
    primer_concentration: f32,        // nM
    salt_concentration: f32,          // mM
    thermodynamic_temp: f32,          // K

    // Output parameters
    max_primers_output: usize,
    min_temp_diff: f32,

    // Algorithm flags
    measure_by_amplicon: bool,
    non_degenerate: bool,

    // File paths
    input_file: Option<String>,
    output_file: Option<String>,

    // Region selection
    begin_nucleotide: Option<usize>,
    end_nucleotide: Option<usize>,

    // Delta G threshold for repetition filter
    delta_g_min: f32,
}

impl Settings {
    /// Create a new Settings struct with default values
    pub fn new() -> Self {
        Settings::default()
    }

    // Getters
    pub fn min_primer_len(&self) -> usize {
        self.min_primer_len
    }

    pub fn max_primer_len(&self) -> usize {
        self.max_primer_len
    }

    pub fn min_amplicon_len(&self) -> usize {
        self.min_amplicon_len
    }

    pub fn min_tm(&self) -> f32 {
        self.min_tm
    }

    pub fn max_tm(&self) -> f32 {
        self.max_tm
    }

    pub fn min_gc_content(&self) -> f32 {
        self.min_gc_content
    }

    pub fn max_gc_content(&self) -> f32 {
        self.max_gc_content
    }

    pub fn max_gc_end(&self) -> f32 {
        self.max_gc_end
    }

    pub fn max_degeneracy_ratio(&self) -> f32 {
        self.max_degeneracy_ratio
    }

    pub fn max_deletion_ratio(&self) -> f32 {
        self.max_deletion_ratio
    }

    pub fn primer_concentration(&self) -> f32 {
        self.primer_concentration
    }

    pub fn salt_concentration(&self) -> f32 {
        self.salt_concentration
    }

    pub fn thermodynamic_temp(&self) -> f32 {
        self.thermodynamic_temp
    }

    pub fn max_primers_output(&self) -> usize {
        self.max_primers_output
    }

    pub fn min_temp_diff(&self) -> f32 {
        self.min_temp_diff
    }

    pub fn measure_by_amplicon(&self) -> bool {
        self.measure_by_amplicon
    }

    pub fn non_degenerate(&self) -> bool {
        self.non_degenerate
    }

    pub fn input_file(&self) -> Option<&str> {
        self.input_file.as_deref()
    }

    pub fn output_file(&self) -> Option<&str> {
        self.output_file.as_deref()
    }

    pub fn begin_nucleotide(&self) -> Option<usize> {
        self.begin_nucleotide
    }

    pub fn end_nucleotide(&self) -> Option<usize> {
        self.end_nucleotide
    }

    pub fn delta_g_min(&self) -> f32 {
        self.delta_g_min
    }

    // Setters
    pub fn with_primer_len(mut self, min: usize, max: usize) -> Self {
        self.min_primer_len = min.max(18).min(25);
        self.max_primer_len = max.max(18).min(25);
        self
    }

    pub fn with_amplicon_len(mut self, len: usize) -> Self {
        self.min_amplicon_len = len;
        self
    }

    pub fn with_tm_range(mut self, min: f32, max: f32) -> Self {
        self.min_tm = min.max(50.0).min(65.0);
        self.max_tm = max.max(50.0).min(65.0);
        self
    }

    pub fn with_gc_content(mut self, min: f32, max: f32) -> Self {
        self.min_gc_content = (min / 100.0).max(0.0).min(1.0);
        self.max_gc_content = (max / 100.0).max(0.0).min(1.0);
        self
    }

    pub fn with_primer_concentration(mut self, conc: f32) -> Self {
        self.primer_concentration = conc.max(50.0);
        self
    }

    pub fn with_salt_concentration(mut self, conc: f32) -> Self {
        self.salt_concentration = conc.max(50.0);
        self
    }

    pub fn with_input_file(mut self, path: String) -> Self {
        self.input_file = Some(path);
        self
    }

    pub fn with_output_file(mut self, path: String) -> Self {
        self.output_file = Some(path);
        self
    }

    pub fn with_region(mut self, begin: usize, end: usize) -> Self {
        self.begin_nucleotide = Some(begin);
        self.end_nucleotide = Some(end);
        self.measure_by_amplicon = false;
        self
    }

    pub fn with_max_primers(mut self, count: usize) -> Self {
        self.max_primers_output = count.min(10);
        self
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.min_primer_len >= self.max_primer_len {
            return Err("Minimum primer length must be less than maximum".to_string());
        }

        if self.min_tm >= self.max_tm {
            return Err("Minimum Tm must be less than maximum".to_string());
        }

        if self.min_tm < 50.0 || self.max_tm > 65.0 {
            return Err("Tm must be between 50.0 and 65.0°C".to_string());
        }

        if self.min_gc_content >= self.max_gc_content {
            return Err("Minimum GC content must be less than maximum".to_string());
        }

        Ok(())
    }
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            min_primer_len: 20,
            max_primer_len: 22,
            min_amplicon_len: 100,
            min_tm: 58.0,
            max_tm: 62.0,
            min_gc_content: 0.40,
            max_gc_content: 0.60,
            max_gc_end: 0.80,
            max_degeneracy_ratio: 0.15,
            max_deletion_ratio: 0.15,
            primer_concentration: 200.0,
            salt_concentration: 50.0,
            thermodynamic_temp: 273.15 + 37.0,
            max_primers_output: 5,
            min_temp_diff: 5.0,
            measure_by_amplicon: true,
            non_degenerate: true,
            input_file: None,
            output_file: None,
            begin_nucleotide: None,
            end_nucleotide: None,
            delta_g_min: -4.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_settings() {
        let settings = Settings::default();
        assert!(settings.validate().is_ok());
    }

    #[test]
    fn test_invalid_tm_range() {
        let settings = Settings::default().with_tm_range(45.0, 70.0);
        // Should clamp to valid range
        assert!(settings.min_tm >= 50.0);
        assert!(settings.max_tm <= 65.0);
    }

    #[test]
    fn test_settings_builder_pattern() {
        let settings = Settings::new()
            .with_primer_len(19, 24)
            .with_tm_range(56.0, 63.0)
            .with_input_file("test.fasta".to_string());

        assert_eq!(settings.min_primer_len(), 19);
        assert_eq!(settings.max_primer_len(), 24);
        assert!(settings.input_file().is_some());
    }
}
