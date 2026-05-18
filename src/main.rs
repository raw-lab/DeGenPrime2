use clap::Parser;
use std::path::PathBuf;
use std::time::Instant;

use degenprime_rs::{
    calculator::PrimerCalculator,
    io::{ResultWriter, SequenceReader},
    models::{DataNode, DataSequence},
    settings::Settings,
};

/// DeGenPrime - Robust PCR primer design for phylogenetically similar sequences
#[derive(Parser, Debug)]
#[command(
    name = "DeGenPrime",
    version = "2.0.0",
    about = "Design optimal PCR primers from aligned DNA sequences",
    long_about = "DeGenPrime selects the top PCR primer pairs for phylogenetically similar DNA sequences based on filter criteria including low degeneracy, few deletions, GC content, low repetition, non-complementary ends, and minimal risk of secondary structures."
)]
struct Args {
    /// Input FASTA or CLUSTAL alignment file
    #[arg(long, value_name = "FILE")]
    input_file: Option<PathBuf>,

    /// Output CSV file for primer results
    #[arg(long, value_name = "FILE")]
    output_file: Option<PathBuf>,

    /// Minimum primer length (18-25, default 20)
    #[arg(long, value_name = "INT", default_value = "20")]
    min_primer_len: usize,

    /// Maximum primer length (18-25, default 22)
    #[arg(long, value_name = "INT", default_value = "22")]
    max_primer_len: usize,

    /// Minimum amplicon length in base pairs
    #[arg(long, value_name = "INT", default_value = "100")]
    amplicon: usize,

    /// Beginning nucleotide position
    #[arg(long, value_name = "INT")]
    begin: Option<usize>,

    /// Ending nucleotide position
    #[arg(long, value_name = "INT")]
    end: Option<usize>,

    /// Minimum melting temperature (50-65°C, default 58)
    #[arg(long, value_name = "FLOAT", default_value = "58.0")]
    min_temp: f32,

    /// Maximum melting temperature (50-65°C, default 62)
    #[arg(long, value_name = "FLOAT", default_value = "62.0")]
    max_temp: f32,

    /// Minimum GC content percentage (0-100, default 40)
    #[arg(long, value_name = "FLOAT", default_value = "40")]
    min_gc: f32,

    /// Maximum GC content percentage (0-100, default 60)
    #[arg(long, value_name = "FLOAT", default_value = "60")]
    max_gc: f32,

    /// Primer concentration in nM (default 200)
    #[arg(long, value_name = "FLOAT", default_value = "200")]
    primer_conc: f32,

    /// Salt concentration in mM (default 50)
    #[arg(long, value_name = "FLOAT", default_value = "50")]
    salt_conc: f32,

    /// Maximum number of primer pairs to output (1-10, default 5)
    #[arg(long, value_name = "INT", default_value = "5")]
    max_primers: usize,

    /// Test a specific primer sequence
    #[arg(long, value_name = "STRING")]
    test: Option<String>,

    /// Analyze a primer sequence and show filter results
    #[arg(long, value_name = "STRING")]
    analyze: Option<String>,

    /// Show detailed report in addition to CSV
    #[arg(long)]
    detailed: bool,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let start = Instant::now();

    // Validate and build settings
    let mut settings = Settings::new();

    if let Some(input) = &args.input_file {
        settings = settings.with_input_file(input.to_string_lossy().to_string());
    } else {
        eprintln!("Error: --input_file is required");
        std::process::exit(1);
    }

    settings = settings
        .with_primer_len(args.min_primer_len, args.max_primer_len)
        .with_amplicon_len(args.amplicon)
        .with_tm_range(args.min_temp, args.max_temp)
        .with_gc_content(args.min_gc, args.max_gc)
        .with_primer_concentration(args.primer_conc)
        .with_salt_concentration(args.salt_conc)
        .with_max_primers(args.max_primers);

    if let (Some(begin), Some(end)) = (args.begin, args.end) {
        settings = settings.with_region(begin, end);
    }

    if let Some(output) = &args.output_file {
        settings = settings.with_output_file(output.to_string_lossy().to_string());
    }

    settings.validate()?;

    if args.verbose {
        println!("DeGenPrime 2.0 - Rust Implementation");
        println!("====================================");
        println!("Settings:");
        println!("  Primer length: {}-{}", settings.min_primer_len(), settings.max_primer_len());
        println!("  Tm range: {:.1}°C - {:.1}°C", settings.min_tm(), settings.max_tm());
        println!("  GC content: {:.0}% - {:.0}%", settings.min_gc_content() * 100.0, settings.max_gc_content() * 100.0);
        println!();
    }

    // Handle test mode
    if let Some(test_seq) = args.test {
        return handle_test_mode(&test_seq, &settings);
    }

    // Handle analyze mode
    if let Some(analyze_seq) = args.analyze {
        return handle_analyze_mode(&analyze_seq, &settings);
    }

    // Read input file
    let input_path = settings
        .input_file()
        .ok_or("Input file is required")?;

    if args.verbose {
        println!("Reading sequences from: {}", input_path);
    }

    let sequences = SequenceReader::read(std::path::Path::new(input_path))?;

    if args.verbose {
        println!("Found {} sequences", sequences.len());
        println!("Sequence lengths: {:?}", 
            sequences.iter().map(|s| s.len()).collect::<Vec<_>>());
    }

    // Check alignment
    if !SequenceReader::is_aligned(&sequences) {
        eprintln!("Error: Sequences are not aligned.");
        eprintln!("Please align your sequences using MAFFT or similar tool:");
        eprintln!("  mafft --localpair --maxiterate 1000 --clustalout input.fasta > output.clust");
        std::process::exit(1);
    }

    // Convert to DataSequence
    let consensus = SequenceReader::alignment_to_consensus(&sequences)?;

    if args.verbose {
        println!("Consensus sequence length: {} bp", consensus.len());
    }

    // Get reverse complement
    let rev_consensus = consensus.reverse_complement();

    // Design primers
    if args.verbose {
        println!("Generating candidate primers...");
    }

    let primer_pairs = PrimerCalculator::design_primers(&consensus, &rev_consensus, &settings)?;

    if args.verbose {
        println!("Found {} primer pairs", primer_pairs.len());
    }

    if primer_pairs.is_empty() {
        eprintln!("No suitable primer pairs found.");
        eprintln!("Try relaxing filter criteria or checking sequence quality.");
        std::process::exit(1);
    }

    // Write output
    let output_path = if let Some(output) = args.output_file {
        output
    } else {
        let input = std::path::Path::new(input_path);
        input
            .with_extension("csv")
            .to_path_buf()
    };

    if args.verbose {
        println!("Writing results to: {}", output_path.display());
    }

    ResultWriter::write_csv(&output_path, &primer_pairs, &consensus, &rev_consensus)?;

    if args.detailed {
        let detail_path = output_path.with_extension("detailed.txt");
        ResultWriter::write_detailed_report(&detail_path, &primer_pairs, &consensus, &rev_consensus)?;
        if args.verbose {
            println!("Detailed report: {}", detail_path.display());
        }
    }

    // Print summary
    let elapsed = start.elapsed();

    println!("\n✓ DeGenPrime Complete");
    println!("=====================");
    println!("Found {} primer pairs", primer_pairs.len());
    println!("Output file: {}", output_path.display());
    println!("Runtime: {:.2}s", elapsed.as_secs_f32());

    Ok(())
}

/// Handle test mode - evaluate a single primer
fn handle_test_mode(primer_seq: &str, settings: &Settings) -> Result<(), Box<dyn std::error::Error>> {
    println!("Testing primer: {}", primer_seq);
    println!("===============================\n");

    let normalized = primer_seq.to_uppercase();
    let chars: Vec<char> = normalized.chars().collect();
    let nodes: Vec<_> = chars.iter()
        .map(|c| DataNode::new(vec![*c]))
        .collect();

    let seq = DataSequence::new(nodes);
    let report = PrimerCalculator::analyze_primer(&seq, settings);
    println!("{}", report);

    Ok(())
}

/// Handle analyze mode - detailed analysis of a sequence
fn handle_analyze_mode(seq_str: &str, settings: &Settings) -> Result<(), Box<dyn std::error::Error>> {
    println!("Analyzing sequence: {}", seq_str);
    println!("===============================\n");

    let normalized = seq_str.to_uppercase();
    let chars: Vec<char> = normalized.chars().collect();
    let nodes: Vec<_> = chars.iter()
        .map(|c| DataNode::new(vec![*c]))
        .collect();

    let seq = DataSequence::new(nodes);
    let report = PrimerCalculator::analyze_primer(&seq, settings);
    println!("{}", report);

    Ok(())
}
