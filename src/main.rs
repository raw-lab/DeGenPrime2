//! DeGenPrime command-line interface.
//!
//! Two principal modes:
//!   * design — read a FASTA/CLUSTAL alignment and emit primer pairs (CSV + JSON)
//!   * test   — evaluate a single user-provided primer sequence
//!
//! Output is always machine-readable so the Python GUI can pick it up. The
//! `--detailed` flag adds a human-friendly text report as well.

use clap::{Parser, Subcommand};
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use degenprime::{
    DataSequence, DesignReport, PrimerCalculator, ResultWriter, SequenceReader, Settings,
};

#[derive(Parser, Debug)]
#[command(
    name = "degenprime",
    version = "2.1.0",
    about = "Robust PCR primer design for phylogenetically similar DNA sequences",
    long_about = "Selects the top PCR primer pairs from an aligned set of DNA sequences. \
                  Filters: low degeneracy, few deletions, GC content, low repetition, \
                  non-complementary ends, melting temperature. Now with proper reverse-strand \
                  coordinates, settings-aware Tm filter, and JSON output for downstream tools."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Design primer pairs from an aligned FASTA or CLUSTAL file.
    Design(DesignArgs),
    /// Test/analyse a single primer sequence (sequence given on the CLI).
    Test(TestArgs),
}

#[derive(Parser, Debug)]
#[command(rename_all = "snake_case")]
struct DesignArgs {
    /// Input FASTA or CLUSTAL alignment file
    #[arg(short, long, value_name = "FILE")]
    input: PathBuf,

    /// Base name for output files; ".csv" / ".json" / ".txt" will be appended.
    /// If omitted, derived from the input file name.
    #[arg(short, long, value_name = "PATH")]
    output: Option<PathBuf>,

    /// Minimum primer length (18-25)
    #[arg(long, default_value_t = 20)]
    min_primer_len: usize,
    /// Maximum primer length (18-25)
    #[arg(long, default_value_t = 22)]
    max_primer_len: usize,

    /// Minimum amplicon length (bp)
    #[arg(long, default_value_t = 100)]
    amplicon: usize,
    /// Optional maximum amplicon length (bp)
    #[arg(long)]
    max_amplicon: Option<usize>,

    /// Begin nucleotide (0-indexed)
    #[arg(long)]
    begin: Option<usize>,
    /// End nucleotide (exclusive)
    #[arg(long)]
    end: Option<usize>,

    /// Min melting temperature (50-65 °C)
    #[arg(long, default_value_t = 58.0)]
    min_temp: f32,
    /// Max melting temperature (50-65 °C)
    #[arg(long, default_value_t = 62.0)]
    max_temp: f32,

    /// Min GC content percentage (0-100)
    #[arg(long, default_value_t = 40.0)]
    min_gc: f32,
    /// Max GC content percentage (0-100)
    #[arg(long, default_value_t = 60.0)]
    max_gc: f32,

    /// Primer concentration (nM)
    #[arg(long, default_value_t = 200.0)]
    primer_conc: f32,
    /// Salt (Na+) concentration (mM)
    #[arg(long, default_value_t = 50.0)]
    salt_conc: f32,

    /// Maximum allowed ΔTm between forward and reverse (°C)
    #[arg(long, default_value_t = 5.0)]
    max_tm_diff: f32,

    /// Maximum number of primer pairs to output
    #[arg(long, default_value_t = 5)]
    max_primers: usize,

    /// Write a human-readable detailed text report
    #[arg(long)]
    detailed: bool,

    /// Verbose progress output
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Parser, Debug)]
#[command(rename_all = "snake_case")]
struct TestArgs {
    /// Primer sequence (e.g. ATGCAATTGGCC)
    #[arg(short, long)]
    primer: String,

    /// Optional JSON output path
    #[arg(long)]
    json: Option<PathBuf>,

    /// Primer concentration (nM)
    #[arg(long, default_value_t = 200.0)]
    primer_conc: f32,
    /// Salt concentration (mM)
    #[arg(long, default_value_t = 50.0)]
    salt_conc: f32,
    /// Min Tm (°C)
    #[arg(long, default_value_t = 50.0)]
    min_temp: f32,
    /// Max Tm (°C)
    #[arg(long, default_value_t = 65.0)]
    max_temp: f32,
    /// Min GC (%)
    #[arg(long, default_value_t = 30.0)]
    min_gc: f32,
    /// Max GC (%)
    #[arg(long, default_value_t = 70.0)]
    max_gc: f32,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Design(args) => run_design(args),
        Command::Test(args) => run_test(args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run_design(args: DesignArgs) -> Result<(), String> {
    let start = Instant::now();

    let mut settings = Settings::new()
        .with_primer_len(args.min_primer_len, args.max_primer_len)
        .with_amplicon_len(args.amplicon)
        .with_tm_range(args.min_temp, args.max_temp)
        .with_gc_content(args.min_gc, args.max_gc)
        .with_primer_concentration(args.primer_conc)
        .with_salt_concentration(args.salt_conc)
        .with_max_primers(args.max_primers)
        .with_max_tm_diff(args.max_tm_diff)
        .with_input_file(args.input.to_string_lossy().into_owned());
    if let Some(m) = args.max_amplicon {
        settings = settings.with_max_amplicon_len(m);
    }
    if let (Some(b), Some(e)) = (args.begin, args.end) {
        settings = settings.with_region(b, e);
    }
    settings.validate()?;

    if args.verbose {
        eprintln!("DeGenPrime 2.1 — designing primers from {}", args.input.display());
    }

    let mut sequences = SequenceReader::read(&args.input)?;
    if sequences.is_empty() {
        return Err("input file contained no sequences".into());
    }

    if !SequenceReader::is_aligned(&sequences) {
        if args.verbose {
            eprintln!("warning: sequences not aligned; padding to longest with gaps. \
                       Use MAFFT for proper alignment: \
                       `mafft --localpair --maxiterate 1000 --clustalout {} > aligned.clust`",
                       args.input.display());
        }
        SequenceReader::pad_to_longest(&mut sequences);
    }

    let consensus = SequenceReader::alignment_to_consensus(&sequences)?;
    let rev_consensus = consensus.reverse_complement();
    let template_len = consensus.len();

    if args.verbose {
        eprintln!("Read {} sequences, consensus length {} bp", sequences.len(), template_len);
    }

    let pairs = PrimerCalculator::design_primers(&consensus, &rev_consensus, &settings)?;
    if args.verbose {
        eprintln!("Found {} primer pair{}", pairs.len(), if pairs.len() == 1 {""} else {"s"});
    }

    let rows = ResultWriter::build_rows(&pairs, &consensus, &rev_consensus, &settings);

    // Determine output paths.
    let base = args
        .output
        .clone()
        .unwrap_or_else(|| args.input.with_extension(""));
    let csv_path = base.with_extension("csv");
    let json_path = base.with_extension("json");

    settings.output_file = Some(csv_path.to_string_lossy().into_owned());

    let report = DesignReport {
        version: env!("CARGO_PKG_VERSION").into(),
        generated_at: chrono_now(),
        input_file: Some(args.input.to_string_lossy().into_owned()),
        template_length: template_len,
        settings: settings.clone(),
        primer_pairs: rows.clone(),
    };

    ResultWriter::write_csv(&csv_path, &rows)?;
    ResultWriter::write_json(&json_path, &report)?;
    if args.detailed {
        let txt_path = base.with_extension("txt");
        ResultWriter::write_detailed_report(&txt_path, &rows, &settings, template_len)?;
        if args.verbose {
            eprintln!("Detailed report: {}", txt_path.display());
        }
    }

    println!();
    println!("Designed {} primer pair{} in {:.2}s", pairs.len(),
        if pairs.len() == 1 {""} else {"s"}, start.elapsed().as_secs_f32());
    println!("CSV:  {}", csv_path.display());
    println!("JSON: {}", json_path.display());
    Ok(())
}

fn run_test(args: TestArgs) -> Result<(), String> {
    let settings = Settings::new()
        .with_primer_concentration(args.primer_conc)
        .with_salt_concentration(args.salt_conc)
        .with_tm_range(args.min_temp, args.max_temp)
        .with_gc_content(args.min_gc, args.max_gc);
    let seq = DataSequence::from_string(&args.primer);
    let analysis = PrimerCalculator::analyze_primer(&seq, &settings);

    println!("{}", analysis.report());
    if let Some(p) = &args.json {
        ResultWriter::write_primer_analysis_json(p, &analysis)?;
        println!("JSON written to {}", p.display());
    }
    Ok(())
}

/// Format current UTC time as ISO-8601 without pulling in a heavy dep.
fn chrono_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let now = SystemTime::now().duration_since(UNIX_EPOCH).unwrap_or_default();
    let secs = now.as_secs();
    // Days since epoch
    let days = secs / 86400;
    let h = ((secs % 86400) / 3600) as u32;
    let m = ((secs % 3600) / 60) as u32;
    let s = (secs % 60) as u32;
    // Civil-from-days algorithm (Hinnant)
    let z = days as i64 + 719468;
    let era = z.div_euclid(146097);
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mo <= 2 { y + 1 } else { y };
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{m:02}:{s:02}Z")
}
