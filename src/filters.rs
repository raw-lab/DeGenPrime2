use crate::models::{DataSequence, Primer, DataNode};
use crate::settings::Settings;

/// Result of applying a filter
#[derive(Debug, Clone)]
pub struct FilterResult {
    pub passed: bool,
    pub reason: Option<String>,
}

impl FilterResult {
    pub fn pass() -> Self {
        FilterResult {
            passed: true,
            reason: None,
        }
    }

    pub fn fail(reason: impl Into<String>) -> Self {
        FilterResult {
            passed: false,
            reason: Some(reason.into()),
        }
    }
}

/// Compose multiple filters using function composition
pub struct FilterChain {
    filters: Vec<Box<dyn Fn(&DataSequence, &Settings) -> FilterResult>>,
}

impl FilterChain {
    pub fn new() -> Self {
        FilterChain {
            filters: Vec::new(),
        }
    }

    pub fn add<F>(mut self, filter: F) -> Self
    where
        F: Fn(&DataSequence, &Settings) -> FilterResult + 'static,
    {
        self.filters.push(Box::new(filter));
        self
    }

    /// Apply all filters, short-circuit on first failure
    pub fn apply_all(&self, seq: &DataSequence, settings: &Settings) -> FilterResult {
        self.filters
            .iter()
            .find_map(|filter| {
                let result = filter(seq, settings);
                if !result.passed {
                    Some(result)
                } else {
                    None
                }
            })
            .unwrap_or_else(FilterResult::pass)
    }
}

impl Default for FilterChain {
    fn default() -> Self {
        Self::new()
    }
}

/// Filter: Degeneracy - check for too many ambiguous bases
pub fn filter_degeneracy(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let degenerate_count = seq.degeneracy_count();
    let ratio = degenerate_count as f32 / seq.actual_size().max(1) as f32;

    if ratio > settings.max_degeneracy_ratio() {
        FilterResult::fail(format!(
            "Too much degeneracy: {} positions out of {}",
            degenerate_count,
            seq.actual_size()
        ))
    } else {
        FilterResult::pass()
    }
}

/// Filter: Deletions - check for too many gaps
pub fn filter_deletions(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let gap_count = seq.gap_count();
    let ratio = gap_count as f32 / seq.len() as f32;

    if ratio > settings.max_deletion_ratio() {
        FilterResult::fail(format!(
            "Too many deletions: {} gaps out of {}",
            gap_count,
            seq.len()
        ))
    } else {
        FilterResult::pass()
    }
}

/// Filter: GC Content - check for optimal GC content
pub fn filter_gc_content(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let gc = seq.gc_ratio();
    let min_gc = settings.min_gc_content();
    let max_gc = settings.max_gc_content();

    if gc < min_gc || gc > max_gc {
        FilterResult::fail(format!(
            "GC content out of range: {:.1}% (expected {:.1}%-{:.1}%)",
            gc * 100.0,
            min_gc * 100.0,
            max_gc * 100.0
        ))
    } else {
        // Check GC at the 3' end (last 5 nucleotides)
        if seq.len() >= 5 {
            let end = seq.sub_seq(seq.len() - 5, 5);
            let end_gc = end.gc_ratio();
            if end_gc > settings.max_gc_end() {
                return FilterResult::fail(format!(
                    "Too much GC at 3' end: {:.1}% (max {:.1}%)",
                    end_gc * 100.0,
                    settings.max_gc_end() * 100.0
                ));
            }
        }
        FilterResult::pass()
    }
}

/// Filter: Repetition - check for homopolymer runs and tandem repeats
pub fn filter_repeats(seq: &DataSequence, _settings: &Settings) -> FilterResult {
    let codes = seq.codes();

    // Check for homopolymer runs (4+ identical bases in a row)
    let max_homopolymer = 4;
    let mut current_char = '\0';
    let mut count = 0;

    for c in codes.chars() {
        if c == current_char {
            count += 1;
            if count > max_homopolymer {
                return FilterResult::fail(format!(
                    "Homopolymer run: {} bases of '{}'",
                    count, c
                ));
            }
        } else {
            current_char = c;
            count = 1;
        }
    }

    // Check for tandem repeats (dinucleotide repeats)
    let has_tandem_repeat = codes
        .as_bytes()
        .windows(4)
        .any(|w| w[0] == w[2] && w[1] == w[3]);

    if has_tandem_repeat {
        FilterResult::fail("Tandem repeats detected".to_string())
    } else {
        FilterResult::pass()
    }
}

/// Filter: Complementary ends - check for potential secondary structures
pub fn filter_complementary_ends(seq: &DataSequence, _settings: &Settings) -> FilterResult {
    if seq.has_complementary_ends() {
        FilterResult::fail("Complementary ends detected - potential for secondary structure".to_string())
    } else {
        FilterResult::pass()
    }
}

/// Filter: Temperature - check if melting temperature is in acceptable range
pub fn filter_temperature(seq: &DataSequence, settings: &Settings) -> FilterResult {
    let tm = calculate_tm(seq);
    let min_tm = settings.min_tm();
    let max_tm = settings.max_tm();

    if tm < min_tm || tm > max_tm {
        FilterResult::fail(format!(
            "Melting temperature out of range: {:.1}°C (expected {:.1}°C-{:.1}°C)",
            tm, min_tm, max_tm
        ))
    } else {
        FilterResult::pass()
    }
}

/// Calculate primer melting temperature using nearest-neighbor method
/// Based on SantaLucia (1998) DNA thermodynamics
pub fn calculate_tm(seq: &DataSequence) -> f32 {
    let nodes = seq.nodes();
    if nodes.is_empty() {
        return 0.0;
    }

    // Calculate thermodynamic parameters using nearest-neighbor model
    let mut dh = -3.1; // Terminal AT penalty
    let mut ds = -12.6;

    if nodes.len() > 1 {
        for window in nodes.windows(2) {
            dh += window[0].enthalpy_to(&window[1]);
            ds += window[0].entropy_to(&window[1]);
        }
    } else {
        dh += nodes[0].enthalpy_to(&nodes[0]);
        ds += nodes[0].entropy_to(&nodes[0]);
    }

    // Calculate salt correction (from Owczarzy et al. 2008)
    let na_conc = 50.0; // mM - standard PCR salt concentration
    let k_corr = 12.0 * ((na_conc / 1000.0).ln());

    // Tm = (dH / (dS + R*ln(C/2))) - 273.15
    // where R = 1.987 cal/(mol*K) and C is primer concentration
    let r = 1.987;
    let tm = dh * 1000.0 / (ds + r * ((50.0 / 1_000_000_000.0) / 2.0).ln()) - 273.15 + k_corr;

    tm.max(0.0)
}

/// Check if a primer sequence contains only valid IUPAC codes
pub fn is_valid_iupac(seq: &DataSequence) -> bool {
    let valid_chars = "AaCcGgTtWwSsRrYyKkMmBbDdHhVvNn-";
    seq.codes().chars().all(|c| valid_chars.contains(c))
}

/// Filter builder for more complex filtering scenarios
pub struct FilterBuilder {
    min_primer_len: usize,
    max_primer_len: usize,
    min_degeneracy: f32,
    min_deletions: f32,
}

impl FilterBuilder {
    pub fn new() -> Self {
        FilterBuilder {
            min_primer_len: 18,
            max_primer_len: 25,
            min_degeneracy: 0.3,
            min_deletions: 0.3,
        }
    }

    pub fn build(self) -> FilterChain {
        let min_len = self.min_primer_len;
        let max_len = self.max_primer_len;

        FilterChain::new()
            .add(move |seq, _| {
                if seq.actual_size() < min_len || seq.actual_size() > max_len {
                    FilterResult::fail(format!(
                        "Length out of range: {} (expected {}-{})",
                        seq.actual_size(),
                        min_len,
                        max_len
                    ))
                } else {
                    FilterResult::pass()
                }
            })
            .add(filter_degeneracy)
            .add(filter_deletions)
            .add(filter_gc_content)
            .add(filter_repeats)
            .add(filter_complementary_ends)
            .add(filter_temperature)
    }
}

impl Default for FilterBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Apply a set of filters to a sequence and return first failure
pub fn apply_all_filters(seq: &DataSequence, settings: &Settings) -> FilterResult {
    [
        filter_degeneracy,
        filter_deletions,
        filter_gc_content,
        filter_repeats,
        filter_complementary_ends,
        filter_temperature,
    ]
    .iter()
    .find_map(|filter| {
        let result = filter(seq, settings);
        if !result.passed {
            Some(result)
        } else {
            None
        }
    })
    .unwrap_or_else(FilterResult::pass)
}

/// Applicative filter composition - applies filter and returns transformed result
pub fn compose_filters<'a>(
    filters: Vec<Box<dyn Fn(&DataSequence, &Settings) -> FilterResult + 'a>>,
) -> impl Fn(&DataSequence, &Settings) -> Vec<FilterResult> + 'a {
    move |seq, settings| {
        filters
            .iter()
            .map(|filter| filter(seq, settings))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::DataNode;

    #[test]
    fn test_filter_gc_content() {
        let settings = Settings::default();
        let nodes = vec![
            DataNode::new(vec!['G']),
            DataNode::new(vec!['C']),
            DataNode::new(vec!['A']),
            DataNode::new(vec!['T']),
        ];
        let seq = DataSequence::new(nodes);
        let result = filter_gc_content(&seq, &settings);
        assert!(result.passed);
    }

    #[test]
    fn test_filter_repeats() {
        let settings = Settings::default();
        let codes = "AAAAA".chars().map(|c| DataNode::new(vec![c])).collect();
        let seq = DataSequence::new(codes);
        let result = filter_repeats(&seq, &settings);
        assert!(!result.passed);
    }
}
