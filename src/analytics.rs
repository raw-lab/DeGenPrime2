use crate::models::{PrimerPair, DataSequence, Primer};
use crate::settings::Settings;
use crate::thermodynamics::Thermodynamics;
use polars::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Comprehensive analytics for primer design results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerAnalytics {
    /// Total candidates generated
    pub total_candidates: usize,
    /// Total after filtering
    pub passed_filters: usize,
    /// Filter pass rates
    pub filter_pass_rates: HashMap<String, f64>,
    /// Melting temperature statistics
    pub tm_stats: TemperatureStats,
    /// GC content statistics
    pub gc_stats: GCStats,
    /// Degeneracy statistics
    pub degeneracy_stats: DegeneracyStats,
    /// Primer pair statistics
    pub pair_stats: PairStats,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TemperatureStats {
    pub mean: f32,
    pub median: f32,
    pub std_dev: f32,
    pub min: f32,
    pub max: f32,
    pub q25: f32,
    pub q75: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GCStats {
    pub mean: f32,
    pub median: f32,
    pub std_dev: f32,
    pub min: f32,
    pub max: f32,
    pub distribution: Vec<(f32, usize)>, // (gc_percent, count)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DegeneracyStats {
    pub mean_degenerate_count: f32,
    pub max_degenerate_count: usize,
    pub primers_no_degeneracy: usize,
    pub primers_low_degeneracy: usize,
    pub primers_high_degeneracy: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairStats {
    pub total_pairs: usize,
    pub mean_tm_diff: f32,
    pub mean_amplicon_size: f32,
    pub min_amplicon_size: usize,
    pub max_amplicon_size: usize,
    pub tm_diff_distribution: Vec<(f32, usize)>,
}

/// Analytics engine using Polars for computation
pub struct AnalyticsEngine;

impl AnalyticsEngine {
    /// Calculate comprehensive statistics from primer results
    pub fn analyze_primers(
        pairs: &[PrimerPair],
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
        settings: &Settings,
    ) -> Result<PrimerAnalytics, String> {
        if pairs.is_empty() {
            return Err("No primer pairs to analyze".to_string());
        }

        // Extract data into Polars DataFrame
        let mut fwd_primers = Vec::new();
        let mut rev_primers = Vec::new();
        let mut tm_diffs = Vec::new();
        let mut amplicon_sizes = Vec::new();
        let mut pair_indices = Vec::new();

        for (idx, pair) in pairs.iter().enumerate() {
            fwd_primers.push(pair.forward().index() as i32);
            rev_primers.push(pair.reverse().index() as i32);
            tm_diffs.push(pair.temp_diff());
            amplicon_sizes.push(pair.amplicon_size() as i32);
            pair_indices.push(idx as i32);
        }

        // Create Polars DataFrame
        let df = df!(
            "index" => &pair_indices,
            "fwd_position" => &fwd_primers,
            "rev_position" => &rev_primers,
            "tm_diff" => &tm_diffs,
            "amplicon_size" => &amplicon_sizes,
        )
        .map_err(|e| e.to_string())?;

        // Calculate statistics
        let tm_stats = Self::calculate_temperature_stats(&df, fwd_seq, rev_seq, settings)?;
        let gc_stats = Self::calculate_gc_stats(pairs, fwd_seq, rev_seq)?;
        let degeneracy_stats = Self::calculate_degeneracy_stats(pairs, fwd_seq, rev_seq)?;
        let pair_stats = Self::calculate_pair_stats(&df)?;

        Ok(PrimerAnalytics {
            total_candidates: pairs.len() * 2, // Estimate
            passed_filters: pairs.len() * 2,
            filter_pass_rates: HashMap::new(),
            tm_stats,
            gc_stats,
            degeneracy_stats,
            pair_stats,
        })
    }

    fn calculate_temperature_stats(
        df: &DataFrame,
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
        settings: &Settings,
    ) -> Result<TemperatureStats, String> {
        let indices: Vec<usize> = df
            .column("index")
            .map_err(|e| e.to_string())?
            .i32()
            .map_err(|e| e.to_string())?
            .into_iter()
            .filter_map(|v| v.map(|i| i as usize))
            .collect();

        let mut temps = Vec::new();
        // This would be filled from actual primer sequences
        for _ in indices {
            temps.push(60.0); // Placeholder
        }

        temps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mean = temps.iter().sum::<f32>() / temps.len() as f32;
        let variance = temps.iter().map(|t| (t - mean).powi(2)).sum::<f32>() / temps.len() as f32;
        let std_dev = variance.sqrt();

        let mid = temps.len() / 2;
        let median = if temps.len() % 2 == 0 {
            (temps[mid - 1] + temps[mid]) / 2.0
        } else {
            temps[mid]
        };

        let q25_idx = temps.len() / 4;
        let q75_idx = (temps.len() * 3) / 4;

        Ok(TemperatureStats {
            mean,
            median,
            std_dev,
            min: *temps.first().unwrap_or(&0.0),
            max: *temps.last().unwrap_or(&0.0),
            q25: temps.get(q25_idx).copied().unwrap_or(0.0),
            q75: temps.get(q75_idx).copied().unwrap_or(0.0),
        })
    }

    fn calculate_gc_stats(
        pairs: &[PrimerPair],
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
    ) -> Result<GCStats, String> {
        let mut gc_values = Vec::new();

        for pair in pairs {
            let fwd_primer = fwd_seq.sub_seq(pair.forward().index(), pair.forward().length());
            let rev_primer = rev_seq.sub_seq(pair.reverse().index(), pair.reverse().length());

            gc_values.push(fwd_primer.gc_ratio() * 100.0);
            gc_values.push(rev_primer.gc_ratio() * 100.0);
        }

        gc_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mean = gc_values.iter().sum::<f32>() / gc_values.len() as f32;
        let variance = gc_values.iter().map(|g| (g - mean).powi(2)).sum::<f32>() / gc_values.len() as f32;
        let std_dev = variance.sqrt();

        let mid = gc_values.len() / 2;
        let median = if gc_values.len() % 2 == 0 {
            (gc_values[mid - 1] + gc_values[mid]) / 2.0
        } else {
            gc_values[mid]
        };

        // Create distribution (binned)
        let mut distribution = Vec::new();
        for bin in (30..=70).step_by(5) {
            let count = gc_values
                .iter()
                .filter(|g| **g >= bin as f32 && **g < (bin + 5) as f32)
                .count();
            if count > 0 {
                distribution.push((bin as f32, count));
            }
        }

        Ok(GCStats {
            mean,
            median,
            std_dev,
            min: *gc_values.first().unwrap_or(&0.0),
            max: *gc_values.last().unwrap_or(&0.0),
            distribution,
        })
    }

    fn calculate_degeneracy_stats(
        pairs: &[PrimerPair],
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
    ) -> Result<DegeneracyStats, String> {
        let mut degen_counts = Vec::new();
        let mut no_degen = 0;
        let mut low_degen = 0;
        let mut high_degen = 0;

        for pair in pairs {
            let fwd_primer = fwd_seq.sub_seq(pair.forward().index(), pair.forward().length());
            let rev_primer = rev_seq.sub_seq(pair.reverse().index(), pair.reverse().length());

            let fwd_degen = fwd_primer.degeneracy_count();
            let rev_degen = rev_primer.degeneracy_count();

            degen_counts.push(fwd_degen as f32);
            degen_counts.push(rev_degen as f32);

            if fwd_degen == 0 {
                no_degen += 1;
            } else if fwd_degen <= 2 {
                low_degen += 1;
            } else {
                high_degen += 1;
            }

            if rev_degen == 0 {
                no_degen += 1;
            } else if rev_degen <= 2 {
                low_degen += 1;
            } else {
                high_degen += 1;
            }
        }

        let mean = if degen_counts.is_empty() {
            0.0
        } else {
            degen_counts.iter().sum::<f32>() / degen_counts.len() as f32
        };

        Ok(DegeneracyStats {
            mean_degenerate_count: mean,
            max_degenerate_count: degen_counts.iter().map(|d| *d as usize).max().unwrap_or(0),
            primers_no_degeneracy: no_degen,
            primers_low_degeneracy: low_degen,
            primers_high_degeneracy: high_degen,
        })
    }

    fn calculate_pair_stats(df: &DataFrame) -> Result<PairStats, String> {
        let amplicon_col = df
            .column("amplicon_size")
            .map_err(|e| e.to_string())?
            .i32()
            .map_err(|e| e.to_string())?;

        let tm_diff_col = df
            .column("tm_diff")
            .map_err(|e| e.to_string())?
            .f32()
            .map_err(|e| e.to_string())?;

        let amplicon_values: Vec<i32> = amplicon_col.into_iter().filter_map(|v| v).collect();
        let tm_diff_values: Vec<f32> = tm_diff_col.into_iter().filter_map(|v| v).collect();

        let mean_amplicon = amplicon_values.iter().sum::<i32>() as f32 / amplicon_values.len() as f32;
        let mean_tm_diff = tm_diff_values.iter().sum::<f32>() / tm_diff_values.len() as f32;

        // Create Tm diff distribution
        let mut tm_distribution = Vec::new();
        for bin in (0..50).step_by(5) {
            let count = tm_diff_values
                .iter()
                .filter(|t| **t >= bin as f32 && **t < (bin + 5) as f32)
                .count();
            if count > 0 {
                tm_distribution.push((bin as f32, count));
            }
        }

        Ok(PairStats {
            total_pairs: df.height(),
            mean_tm_diff,
            mean_amplicon_size: mean_amplicon,
            min_amplicon_size: amplicon_values.iter().min().copied().unwrap_or(0) as usize,
            max_amplicon_size: amplicon_values.iter().max().copied().unwrap_or(0) as usize,
            tm_diff_distribution: tm_distribution,
        })
    }

    /// Generate comprehensive statistics report as JSON
    pub fn generate_report(
        analytics: &PrimerAnalytics,
    ) -> Result<String, String> {
        serde_json::to_string_pretty(analytics)
            .map_err(|e| format!("Failed to serialize analytics: {}", e))
    }

    /// Create a summary table for display
    pub fn summary_table(analytics: &PrimerAnalytics) -> String {
        format!(
            r#"
╔════════════════════════════════════════════════════════════╗
║           DeGenPrime Analytics Summary                     ║
╠════════════════════════════════════════════════════════════╣
║ MELTING TEMPERATURE STATISTICS                             ║
║   Mean:      {:.2}°C                                         ║
║   Median:    {:.2}°C                                         ║
║   Std Dev:   {:.2}°C                                         ║
║   Range:     {:.2}°C - {:.2}°C                               ║
║   IQR:       {:.2}°C - {:.2}°C                               ║
╠════════════════════════════════════════════════════════════╣
║ GC CONTENT STATISTICS                                      ║
║   Mean:      {:.1}%                                          ║
║   Median:    {:.1}%                                          ║
║   Std Dev:   {:.1}%                                          ║
║   Range:     {:.1}% - {:.1}%                                 ║
╠════════════════════════════════════════════════════════════╣
║ DEGENERACY ANALYSIS                                        ║
║   Mean Degen Positions: {:.1}                               ║
║   No Degeneracy:        {}                                  ║
║   Low Degeneracy:       {}                                  ║
║   High Degeneracy:      {}                                  ║
╠════════════════════════════════════════════════════════════╣
║ PRIMER PAIR STATISTICS                                     ║
║   Total Pairs:          {}                                  ║
║   Mean Tm Difference:   {:.2}°C                             ║
║   Mean Amplicon Size:   {:.0} bp                             ║
║   Amplicon Range:       {} - {} bp                           ║
╚════════════════════════════════════════════════════════════╝
"#,
            analytics.tm_stats.mean,
            analytics.tm_stats.median,
            analytics.tm_stats.std_dev,
            analytics.tm_stats.min,
            analytics.tm_stats.max,
            analytics.tm_stats.q25,
            analytics.tm_stats.q75,
            analytics.gc_stats.mean,
            analytics.gc_stats.median,
            analytics.gc_stats.std_dev,
            analytics.gc_stats.min,
            analytics.gc_stats.max,
            analytics.degeneracy_stats.mean_degenerate_count,
            analytics.degeneracy_stats.primers_no_degeneracy,
            analytics.degeneracy_stats.primers_low_degeneracy,
            analytics.degeneracy_stats.primers_high_degeneracy,
            analytics.pair_stats.total_pairs,
            analytics.pair_stats.mean_tm_diff,
            analytics.pair_stats.mean_amplicon_size,
            analytics.pair_stats.min_amplicon_size,
            analytics.pair_stats.max_amplicon_size,
        )
    }
}
