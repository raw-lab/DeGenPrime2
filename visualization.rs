use crate::analytics::{PrimerAnalytics, PairStats};
use crate::models::PrimerPair;
use plotly::common::{Mode, Title};
use plotly::{
    Bar, Box, Histogram, Layout, Plot, Scatter, ScatterMode, Violin
};
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Interactive visualization generator for primer analytics
pub struct VisualizationEngine;

impl VisualizationEngine {
    /// Generate comprehensive interactive HTML report
    pub fn generate_html_report(
        path: &Path,
        analytics: &PrimerAnalytics,
        pairs: &[PrimerPair],
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(path)?;

        // Write HTML header
        let html_header = Self::generate_html_header();
        file.write_all(html_header.as_bytes())?;

        // Generate and embed plots
        file.write_all(b"<div class='content'>\n")?;
        file.write_all(b"<h1>DeGenPrime Analysis Report</h1>\n")?;

        // Summary statistics
        file.write_all(b"<section class='summary'>\n")?;
        file.write_all(b"<h2>Summary Statistics</h2>\n")?;
        file.write_all(Self::summary_stats_html(analytics).as_bytes())?;
        file.write_all(b"</section>\n")?;

        // Melting temperature distribution
        file.write_all(b"<section class='plot'>\n")?;
        file.write_all(b"<h2>Melting Temperature Distribution</h2>\n")?;
        let tm_plot = Self::generate_tm_histogram(analytics)?;
        file.write_all(tm_plot.as_bytes())?;
        file.write_all(b"</section>\n")?;

        // GC content distribution
        file.write_all(b"<section class='plot'>\n")?;
        file.write_all(b"<h2>GC Content Distribution</h2>\n")?;
        let gc_plot = Self::generate_gc_histogram(analytics)?;
        file.write_all(gc_plot.as_bytes())?;
        file.write_all(b"</section>\n")?;

        // Amplicon size distribution
        file.write_all(b"<section class='plot'>\n")?;
        file.write_all(b"<h2>Amplicon Size Distribution</h2>\n")?;
        let amplicon_plot = Self::generate_amplicon_histogram(analytics)?;
        file.write_all(amplicon_plot.as_bytes())?;
        file.write_all(b"</section>\n")?;

        // Tm difference distribution
        file.write_all(b"<section class='plot'>\n")?;
        file.write_all(b"<h2>Melting Temperature Difference (Forward vs Reverse)</h2>\n")?;
        let tm_diff_plot = Self::generate_tm_diff_histogram(analytics)?;
        file.write_all(tm_diff_plot.as_bytes())?;
        file.write_all(b"</section>\n")?;

        // Scatter plots
        file.write_all(b"<section class='plot'>\n")?;
        file.write_all(b"<h2>Primer Pair Relationships</h2>\n")?;
        let scatter_plot = Self::generate_scatter_plots(pairs)?;
        file.write_all(scatter_plot.as_bytes())?;
        file.write_all(b"</section>\n")?;

        // Box plots for statistics
        file.write_all(b"<section class='plot'>\n")?;
        file.write_all(b"<h2>Statistical Distribution (Box Plots)</h2>\n")?;
        let box_plot = Self::generate_box_plots(analytics)?;
        file.write_all(box_plot.as_bytes())?;
        file.write_all(b"</section>\n")?;

        file.write_all(b"</div>\n")?;

        // Write HTML footer
        let html_footer = Self::generate_html_footer();
        file.write_all(html_footer.as_bytes())?;

        Ok(())
    }

    fn generate_html_header() -> String {
        r#"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DeGenPrime Analysis Report</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: #333;
            line-height: 1.6;
        }
        
        .content {
            max-width: 1400px;
            margin: 0 auto;
            padding: 40px 20px;
        }
        
        h1 {
            text-align: center;
            color: white;
            margin-bottom: 40px;
            font-size: 2.5em;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
        }
        
        h2 {
            color: white;
            margin-top: 30px;
            margin-bottom: 20px;
            font-size: 1.8em;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
        }
        
        section {
            background: white;
            border-radius: 12px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            transition: transform 0.3s ease, box-shadow 0.3s ease;
        }
        
        section:hover {
            transform: translateY(-5px);
            box-shadow: 0 15px 50px rgba(0,0,0,0.3);
        }
        
        .summary {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
        }
        
        .stat-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        }
        
        .stat-label {
            font-size: 0.9em;
            opacity: 0.9;
            margin-bottom: 8px;
            text-transform: uppercase;
            letter-spacing: 1px;
            font-weight: 600;
        }
        
        .stat-value {
            font-size: 1.8em;
            font-weight: bold;
        }
        
        .stat-unit {
            font-size: 0.8em;
            opacity: 0.8;
            margin-top: 5px;
        }
        
        .plot {
            background: white;
            border-radius: 12px;
            padding: 20px;
            margin-bottom: 30px;
        }
        
        .plotly-graph-div {
            width: 100% !important;
            height: 500px !important;
        }
        
        footer {
            text-align: center;
            color: white;
            padding: 20px;
            opacity: 0.8;
            margin-top: 40px;
        }
        
        @media (max-width: 768px) {
            h1 {
                font-size: 1.8em;
            }
            
            h2 {
                font-size: 1.3em;
            }
            
            section {
                padding: 15px;
            }
            
            .plotly-graph-div {
                height: 400px !important;
            }
        }
    </style>
</head>
<body>
"#.to_string()
    }

    fn generate_html_footer() -> String {
        r#"
    <footer>
        <p>DeGenPrime Interactive Report | Generated with Rust + Plotly</p>
        <p>Robust PCR Primer Design Tool</p>
    </footer>
</body>
</html>
"#.to_string()
    }

    fn summary_stats_html(analytics: &PrimerAnalytics) -> String {
        format!(
            r#"
<div class="stat-card">
    <div class="stat-label">Mean Tm</div>
    <div class="stat-value">{:.2}</div>
    <div class="stat-unit">°C</div>
</div>

<div class="stat-card">
    <div class="stat-label">Tm Std Dev</div>
    <div class="stat-value">{:.2}</div>
    <div class="stat-unit">°C</div>
</div>

<div class="stat-card">
    <div class="stat-label">Mean GC</div>
    <div class="stat-value">{:.1}</div>
    <div class="stat-unit">%</div>
</div>

<div class="stat-card">
    <div class="stat-label">Primers</div>
    <div class="stat-value">{}</div>
    <div class="stat-unit">pairs</div>
</div>

<div class="stat-card">
    <div class="stat-label">Mean Amplicon</div>
    <div class="stat-value">{:.0}</div>
    <div class="stat-unit">bp</div>
</div>

<div class="stat-card">
    <div class="stat-label">Mean Tm Diff</div>
    <div class="stat-value">{:.2}</div>
    <div class="stat-unit">°C</div>
</div>
"#,
            analytics.tm_stats.mean,
            analytics.tm_stats.std_dev,
            analytics.gc_stats.mean,
            analytics.pair_stats.total_pairs,
            analytics.pair_stats.mean_amplicon_size,
            analytics.pair_stats.mean_tm_diff
        )
    }

    fn generate_tm_histogram(analytics: &PrimerAnalytics) -> Result<String, Box<dyn std::error::Error>> {
        let mut data = vec![];
        
        // Create histogram data
        let hist = Histogram::new(vec![
            analytics.tm_stats.mean - 5.0,
            analytics.tm_stats.mean - 2.0,
            analytics.tm_stats.mean,
            analytics.tm_stats.mean + 2.0,
            analytics.tm_stats.mean + 5.0,
        ])
        .name("Temperature")
        .n_bins(20);
        
        data.push(hist);

        let layout = Layout::new()
            .title(Title::new("Melting Temperature Distribution"))
            .x_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("Temperature (°C)")))
            .y_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("Frequency")));

        let mut plot = Plot::new();
        plot.set_layout(layout);
        for trace in data {
            plot.add_trace(trace);
        }

        Ok(plot.to_html())
    }

    fn generate_gc_histogram(analytics: &PrimerAnalytics) -> Result<String, Box<dyn std::error::Error>> {
        let gc_values: Vec<f32> = analytics
            .gc_stats
            .distribution
            .iter()
            .flat_map(|(gc, count)| vec![*gc; *count])
            .collect();

        let hist = Histogram::new(gc_values)
            .name("GC Content")
            .n_bins(15);

        let layout = Layout::new()
            .title(Title::new("GC Content Distribution"))
            .x_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("GC Content (%)")))
            .y_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("Frequency")));

        let mut plot = Plot::new();
        plot.add_trace(hist);
        plot.set_layout(layout);

        Ok(plot.to_html())
    }

    fn generate_amplicon_histogram(analytics: &PrimerAnalytics) -> Result<String, Box<dyn std::error::Error>> {
        let amplicon_values: Vec<f32> = vec![
            analytics.pair_stats.min_amplicon_size as f32,
            analytics.pair_stats.mean_amplicon_size,
            analytics.pair_stats.max_amplicon_size as f32,
        ];

        let hist = Histogram::new(amplicon_values)
            .name("Amplicon Size")
            .n_bins(20);

        let layout = Layout::new()
            .title(Title::new("Amplicon Size Distribution"))
            .x_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("Size (bp)")))
            .y_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("Frequency")));

        let mut plot = Plot::new();
        plot.add_trace(hist);
        plot.set_layout(layout);

        Ok(plot.to_html())
    }

    fn generate_tm_diff_histogram(analytics: &PrimerAnalytics) -> Result<String, Box<dyn std::error::Error>> {
        let tm_diff_values: Vec<f32> = analytics
            .pair_stats
            .tm_diff_distribution
            .iter()
            .flat_map(|(tm, count)| vec![*tm; *count])
            .collect();

        let hist = Histogram::new(tm_diff_values)
            .name("Tm Difference")
            .n_bins(15);

        let layout = Layout::new()
            .title(Title::new("Melting Temperature Difference Distribution"))
            .x_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("Tm Difference (°C)")))
            .y_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("Frequency")));

        let mut plot = Plot::new();
        plot.add_trace(hist);
        plot.set_layout(layout);

        Ok(plot.to_html())
    }

    fn generate_scatter_plots(pairs: &[PrimerPair]) -> Result<String, Box<dyn std::error::Error>> {
        let indices: Vec<usize> = (0..pairs.len()).collect();
        let amplicon_sizes: Vec<i32> = pairs.iter().map(|p| p.amplicon_size() as i32).collect();
        let tm_diffs: Vec<f32> = pairs.iter().map(|p| p.temp_diff()).collect();

        let scatter = Scatter::new(
            indices.iter().map(|i| *i as i32).collect(),
            tm_diffs.clone(),
        )
        .mode(ScatterMode::Markers)
        .name("Tm Difference");

        let layout = Layout::new()
            .title(Title::new("Primer Pair Index vs Temperature Difference"))
            .x_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("Pair Index")))
            .y_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("Tm Difference (°C)")));

        let mut plot = Plot::new();
        plot.add_trace(scatter);
        plot.set_layout(layout);

        Ok(plot.to_html())
    }

    fn generate_box_plots(analytics: &PrimerAnalytics) -> Result<String, Box<dyn std::error::Error>> {
        let box_plot = Box::new(vec![
            analytics.tm_stats.min,
            analytics.tm_stats.q25,
            analytics.tm_stats.median,
            analytics.tm_stats.q75,
            analytics.tm_stats.max,
        ])
        .name("Temperature");

        let layout = Layout::new()
            .title(Title::new("Statistical Distribution Summary"))
            .y_axis(plotly::layout::Axis::new().title(plotly::layout::AxisTitle::new("Temperature (°C)")));

        let mut plot = Plot::new();
        plot.add_trace(box_plot);
        plot.set_layout(layout);

        Ok(plot.to_html())
    }
}
