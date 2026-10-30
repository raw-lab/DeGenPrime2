//! `degenprime` — command-line interface.
//!
//! * `design` — design degenerate primer pairs from an alignment (CSV + JSON,
//!   optional text report and FASTA); unaligned input is aligned first with
//!   the built-in, pure-Rust MAFFT
//! * `align` — align sequences with the built-in MAFFT
//! * `test` (alias `check`) — analyse a primer you already have, optionally
//!   against a partner primer
//! * `gui` — open the desktop app (the web interface served on 127.0.0.1)
//! * `doctor` — check that everything works

use clap::{CommandFactory, Parser, Subcommand};
use degenprime::api::{self, AlignRequest, AlignWhen, CheckRequest, DesignRequest, SettingsInput};
use degenprime::io::{DesignReport, PrimerRow, ResultWriter};
use degenprime::msa::AlignMode;
use degenprime::settings::RankBy;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

const EXAMPLES: &str = "Examples:
  degenprime design -i data/example_anfD.clust -o anfd --detailed
  degenprime design -i genes.fasta --max_amplicon 600     (unaligned input is aligned first)
  degenprime align -i genes.fasta -o genes_aligned.fasta  (pure-Rust MAFFT)
  degenprime test -p GCTACCTGAACACCATCCCG --partner CCCTGGATGTTGTASTCGCC
  degenprime gui        open the desktop app in your browser
  degenprime doctor     check your installation";

#[derive(Parser, Debug)]
#[command(
    name = "degenprime",
    version,
    about = "Degenerate PCR primer design from DNA sequences, aligned or not",
    before_help = degenprime::BANNER,
    after_help = EXAMPLES
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Design degenerate primer pairs from FASTA or CLUSTAL (aligned or not)
    Design(DesignArgs),
    /// Align sequences with the built-in MAFFT (pure Rust, no install needed)
    Align(AlignArgs),
    /// Analyse a primer you already have (Tm, GC, degeneracy, dimers)
    #[command(alias = "check")]
    Test(TestArgs),
    /// Open the desktop app in your browser (runs locally; nothing is uploaded)
    #[cfg(feature = "gui")]
    Gui(GuiArgs),
    /// Check that DeGenPrime works on this computer
    Doctor,
}

fn parse_rank(s: &str) -> Result<RankBy, String> {
    RankBy::parse(s).ok_or_else(|| format!("unknown ranking '{s}' (use score or tm_diff)"))
}

fn parse_align_when(s: &str) -> Result<AlignWhen, String> {
    AlignWhen::parse(s).ok_or_else(|| format!("unknown value '{s}' (use auto, always or never)"))
}

fn parse_align_mode(s: &str) -> Result<AlignMode, String> {
    AlignMode::parse(s)
        .ok_or_else(|| format!("unknown MAFFT strategy '{s}' (use auto, fftns2, fftnsi, linsi, ginsi or einsi)"))
}

#[derive(Parser, Debug)]
#[command(rename_all = "snake_case")]
struct DesignArgs {
    /// Input: FASTA or CLUSTAL, aligned or not (unaligned input is aligned first)
    #[arg(short, long, value_name = "FILE")]
    input: PathBuf,
    /// Base name for the outputs (.csv and .json; .txt and .fasta on request)
    #[arg(short, long, value_name = "PATH")]
    output: Option<PathBuf>,
    /// Minimum primer length, nt (18-25)
    #[arg(long, default_value_t = 20)]
    min_primer_len: usize,
    /// Maximum primer length, nt (18-25)
    #[arg(long, default_value_t = 22)]
    max_primer_len: usize,
    /// Minimum product size, bp
    #[arg(long, default_value_t = 100)]
    amplicon: usize,
    /// Maximum product size, bp
    #[arg(long)]
    max_amplicon: Option<usize>,
    /// Design region start: first alignment column, 0-based
    #[arg(long)]
    begin: Option<usize>,
    /// Design region end: alignment column, 0-based, exclusive
    #[arg(long)]
    end: Option<usize>,
    /// Minimum Tm, °C (50-65)
    #[arg(long, default_value_t = 58.0)]
    min_temp: f32,
    /// Maximum Tm, °C (50-65)
    #[arg(long, default_value_t = 62.0)]
    max_temp: f32,
    /// Minimum GC content, %
    #[arg(long, default_value_t = 40.0)]
    min_gc: f32,
    /// Maximum GC content, %
    #[arg(long, default_value_t = 60.0)]
    max_gc: f32,
    /// Maximum GC in the last five 3' bases, %
    #[arg(long, default_value_t = 80.0)]
    max_gc_end: f32,
    /// Maximum share of degenerate positions per primer, %
    #[arg(long, default_value_t = 15.0)]
    max_degeneracy: f32,
    /// Maximum share of gap (deletion) columns per primer, %
    #[arg(long, default_value_t = 15.0)]
    max_deletions: f32,
    /// Longest allowed run of one base
    #[arg(long, default_value_t = 4)]
    max_poly_x: usize,
    /// Maximum self-complementarity, any position (Primer3-style score)
    #[arg(long, default_value_t = 8.0)]
    max_self_any: f32,
    /// Maximum self-complementarity at the 3' end
    #[arg(long, default_value_t = 3.0)]
    max_self_end: f32,
    /// Maximum forward/reverse complementarity, any position
    #[arg(long, default_value_t = 8.0)]
    max_pair_any: f32,
    /// Maximum forward/reverse complementarity at the 3' ends
    #[arg(long, default_value_t = 3.0)]
    max_pair_end: f32,
    /// Primer concentration, nM
    #[arg(long, default_value_t = 200.0)]
    primer_conc: f32,
    /// Monovalent salt (Na+), mM
    #[arg(long, default_value_t = 50.0)]
    salt_conc: f32,
    /// Maximum Tm difference between forward and reverse, °C
    #[arg(long, default_value_t = 5.0)]
    max_tm_diff: f32,
    /// Number of primer pairs to return (1-50)
    #[arg(long, default_value_t = 5)]
    max_primers: usize,
    /// Returned pairs must differ by this many bases at a 3' end (0 = off)
    #[arg(long = "min_3prime_distance", default_value_t = 3)]
    min_3prime_distance: usize,
    /// Ranking: score (conservation, then Tm) or tm_diff (2.1 behaviour)
    #[arg(long, default_value = "score", value_parser = parse_rank)]
    rank: RankBy,
    /// When to align with the built-in MAFFT: auto (when lengths differ), always, never
    #[arg(long, value_name = "WHEN", default_value = "auto", value_parser = parse_align_when)]
    align: AlignWhen,
    /// MAFFT strategy: auto (mafft --auto), fftns2, fftnsi, linsi, ginsi, einsi
    #[arg(long, value_name = "MODE", default_value = "auto", value_parser = parse_align_mode)]
    align_mode: AlignMode,
    /// Don't reverse-complement sequences that lie on the opposite strand
    #[arg(long)]
    keep_strand: bool,
    /// Also write the alignment DeGenPrime made (.aligned.fasta)
    #[arg(long)]
    save_alignment: bool,
    /// Version 2.2 option, now the default (unaligned input is always aligned)
    #[arg(long, hide = true)]
    mafft: bool,
    /// Also write a human-readable report (.txt)
    #[arg(long)]
    detailed: bool,
    /// Also write the primers as FASTA (.fasta)
    #[arg(long)]
    fasta: bool,
    /// Threads to use (0 = all cores)
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
    /// Don't print the summary table
    #[arg(short, long)]
    quiet: bool,
    /// Print progress details to stderr
    #[arg(short, long)]
    verbose: bool,
}

impl DesignArgs {
    fn settings(&self) -> SettingsInput {
        SettingsInput {
            min_primer_len: self.min_primer_len,
            max_primer_len: self.max_primer_len,
            min_amplicon: self.amplicon,
            max_amplicon: self.max_amplicon,
            begin: self.begin,
            end: self.end,
            min_temp: self.min_temp,
            max_temp: self.max_temp,
            max_tm_diff: self.max_tm_diff,
            min_gc: self.min_gc,
            max_gc: self.max_gc,
            max_gc_end: self.max_gc_end,
            max_degeneracy: self.max_degeneracy,
            max_deletions: self.max_deletions,
            max_poly_x: self.max_poly_x,
            max_self_any: self.max_self_any,
            max_self_end: self.max_self_end,
            max_pair_any: self.max_pair_any,
            max_pair_end: self.max_pair_end,
            primer_conc: self.primer_conc,
            salt_conc: self.salt_conc,
            max_primers: self.max_primers,
            min_3prime_distance: self.min_3prime_distance,
            rank_by: self.rank,
        }
    }
}

#[derive(Parser, Debug)]
#[command(rename_all = "snake_case")]
struct AlignArgs {
    /// Sequences to align: FASTA or CLUSTAL (gaps already present are removed)
    #[arg(short, long, value_name = "FILE")]
    input: PathBuf,
    /// Write the alignment here (FASTA); default: standard output
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,
    /// MAFFT strategy: auto (mafft --auto), fftns2, fftnsi, linsi, ginsi, einsi
    #[arg(long, value_name = "MODE", default_value = "auto", value_parser = parse_align_mode)]
    mode: AlignMode,
    /// Don't reverse-complement sequences that lie on the opposite strand
    #[arg(long)]
    keep_strand: bool,
    /// Threads to use (0 = all cores)
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

#[derive(Parser, Debug)]
#[command(rename_all = "snake_case")]
struct TestArgs {
    /// Primer sequence, 5'->3' (IUPAC codes allowed)
    #[arg(short, long)]
    primer: String,
    /// Partner primer, to check the pair (ΔTm and cross-dimers)
    #[arg(long)]
    partner: Option<String>,
    /// Write the analysis as JSON
    #[arg(long)]
    json: Option<PathBuf>,
    /// Primer concentration, nM
    #[arg(long, default_value_t = 200.0)]
    primer_conc: f32,
    /// Monovalent salt (Na+), mM
    #[arg(long, default_value_t = 50.0)]
    salt_conc: f32,
    /// Minimum Tm for the pass/fail check, °C
    #[arg(long, default_value_t = 50.0)]
    min_temp: f32,
    /// Maximum Tm for the pass/fail check, °C
    #[arg(long, default_value_t = 65.0)]
    max_temp: f32,
    /// Minimum GC for the pass/fail check, %
    #[arg(long, default_value_t = 30.0)]
    min_gc: f32,
    /// Maximum GC for the pass/fail check, %
    #[arg(long, default_value_t = 70.0)]
    max_gc: f32,
}

#[cfg(feature = "gui")]
#[derive(Parser, Debug)]
#[command(rename_all = "snake_case")]
struct GuiArgs {
    /// Port on 127.0.0.1 (another free port is used if this one is taken)
    #[arg(long, default_value_t = degenprime::gui::DEFAULT_PORT)]
    port: u16,
    /// Don't open a browser window
    #[arg(long)]
    no_browser: bool,
    /// Threads to use for designs and alignments (0 = all cores)
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        None => {
            let _ = Cli::command().print_help();
            println!();
            return ExitCode::SUCCESS;
        }
        Some(Command::Design(a)) => set_threads(a.threads).and_then(|_| run_design(a)),
        Some(Command::Align(a)) => set_threads(a.threads).and_then(|_| run_align(a)),
        Some(Command::Test(a)) => run_test(a),
        #[cfg(feature = "gui")]
        Some(Command::Gui(a)) => set_threads(a.threads).and_then(|_| {
            degenprime::gui::serve(degenprime::gui::GuiOptions { port: a.port, open_browser: !a.no_browser })
        }),
        Some(Command::Doctor) => run_doctor(),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Size the thread pool used for screening and alignment (0 = all cores).
fn set_threads(n: usize) -> Result<(), String> {
    if n > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|e| format!("cannot start {n} threads: {e}"))?;
    }
    Ok(())
}

fn plural(n: usize) -> &'static str {
    if n == 1 {
        ""
    } else {
        "s"
    }
}

/// Read a text file without copying it when it is valid UTF-8.
fn read_input(path: &Path) -> Result<String, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    Ok(String::from_utf8(bytes).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned()))
}

/// `base` + `.ext`, without doubling a known output extension.
fn out_path(base: &Path, ext: &str) -> PathBuf {
    let s = base.to_string_lossy();
    let trimmed = [".aligned.fasta", ".csv", ".json", ".txt", ".fasta"]
        .iter()
        .find_map(|e| s.strip_suffix(e))
        .unwrap_or(s.as_ref());
    PathBuf::from(format!("{trimmed}.{ext}"))
}

fn same_file(a: &Path, b: &Path) -> bool {
    match (a.canonicalize(), b.canonicalize()) {
        (Ok(x), Ok(y)) => x == y,
        _ => a == b,
    }
}

/// An output path that is never the input file: `genes.fasta` with
/// `--fasta` would otherwise overwrite the alignment it was designed from
/// (version 2.2 did exactly that).
fn safe_out(base: &Path, ext: &str, input: &Path) -> PathBuf {
    let p = out_path(base, ext);
    if same_file(&p, input) {
        let s = out_path(base, "").to_string_lossy().trim_end_matches('.').to_string();
        PathBuf::from(format!("{s}_primers.{ext}"))
    } else {
        p
    }
}

fn print_table(rows: &[PrimerRow]) {
    println!(
        "  {:>2}  {:<26}{:>11}{:>7}   {:<26}{:>11}{:>7}  {:>8}",
        "#", "Forward primer (5'->3')", "Location", "Tm", "Reverse primer (5'->3')", "Location", "Tm", "Product"
    );
    for r in rows {
        println!(
            "  {:>2}  {:<26}{:>11}{:>7.1}   {:<26}{:>11}{:>7.1}  {:>5} bp",
            r.rank,
            r.forward_primer,
            format!("{}-{}", r.forward_from, r.forward_to),
            r.forward_tm,
            r.reverse_primer,
            format!("{}-{}", r.reverse_from, r.reverse_to),
            r.reverse_tm,
            r.amplicon_size
        );
    }
}

fn run_design(args: DesignArgs) -> Result<(), String> {
    let start = Instant::now();
    let input = read_input(&args.input)?;
    if args.verbose {
        eprintln!("DeGenPrime {}: designing primers from {}", degenprime::VERSION, args.input.display());
    }
    let req = DesignRequest {
        input,
        settings: args.settings(),
        filename: args.input.file_name().map(|s| s.to_string_lossy().into_owned()),
        generated_at: Some(chrono_now()),
        align: args.align,
        align_mode: args.align_mode,
        adjust_direction: !args.keep_strand,
        mafft: args.mafft,
        max_site_sequences: Some(0),
    };
    let resp = api::design(&req);
    drop(req);
    for n in &resp.notes {
        eprintln!("note: {n}");
    }
    for w in &resp.warnings {
        eprintln!("warning: {w}");
    }
    if !resp.ok {
        return Err(resp.error.unwrap_or_else(|| "design failed".into()));
    }
    let input = resp.input.as_ref().ok_or("internal error: no input summary")?;
    if args.verbose {
        eprintln!(
            "Read {} sequences ({:?}), {} alignment columns{}",
            input.sequences,
            input.format,
            input.columns,
            input.aligned_with.as_ref().map(|a| format!(", aligned with {a}")).unwrap_or_default()
        );
    }

    let base = args.output.clone().unwrap_or_else(|| args.input.with_extension(""));
    let csv_path = safe_out(&base, "csv", &args.input);
    let json_path = safe_out(&base, "json", &args.input);
    let mut settings = resp.settings.clone().unwrap_or_default();
    settings.input_file = Some(args.input.to_string_lossy().into_owned());
    settings.output_file = Some(csv_path.to_string_lossy().into_owned());
    let report = DesignReport {
        version: degenprime::VERSION.into(),
        generated_at: chrono_now(),
        input_file: settings.input_file.clone(),
        template_length: input.columns,
        aligned_with: input.aligned_with.clone(),
        settings,
        primer_pairs: resp.pairs.clone(),
        warnings: resp.warnings.clone(),
        notes: resp.notes.clone(),
        explain: resp.explain.clone(),
    };
    ResultWriter::write_text(&csv_path, &resp.csv)?;
    ResultWriter::write_json(&json_path, &report)?;
    let mut written = vec![("CSV", csv_path), ("JSON", json_path)];
    if args.detailed {
        let p = safe_out(&base, "txt", &args.input);
        ResultWriter::write_text(&p, &resp.report)?;
        written.push(("Report", p));
    }
    if args.fasta {
        let p = safe_out(&base, "fasta", &args.input);
        ResultWriter::write_text(&p, &resp.fasta)?;
        written.push(("FASTA", p));
    }
    if args.save_alignment {
        match &resp.alignment {
            Some(aln) => {
                let p = safe_out(&base, "aligned.fasta", &args.input);
                ResultWriter::write_text(&p, aln)?;
                written.push(("Aligned", p));
            }
            None => eprintln!("note: the input was already aligned, so there is no new alignment to save"),
        }
    }

    let n = resp.pairs.len();
    if !args.quiet {
        println!();
        println!(
            "Designed {n} primer pair{} in {:.2}s ({} sequence{}, {} alignment columns{})",
            plural(n),
            start.elapsed().as_secs_f32(),
            input.sequences,
            plural(input.sequences),
            input.columns,
            input.aligned_with.as_ref().map(|a| format!(", aligned with {a}")).unwrap_or_default()
        );
        if n > 0 {
            println!();
            print_table(&resp.pairs);
        }
        println!();
    }
    for (label, p) in &written {
        println!("{:<9}{}", format!("{label}:"), p.display());
    }
    if n == 0 {
        return Err("no primer pairs passed the filters (see the warning above)".into());
    }
    Ok(())
}

fn run_align(args: AlignArgs) -> Result<(), String> {
    let start = Instant::now();
    let input = read_input(&args.input)?;
    let resp = api::align(&AlignRequest { input, mode: args.mode, adjust_direction: !args.keep_strand });
    if !resp.ok {
        return Err(resp.error.unwrap_or_else(|| "alignment failed".into()));
    }
    if !resp.reversed.is_empty() {
        eprintln!(
            "note: reverse-complemented {} sequence{} on the opposite strand: {}",
            resp.reversed.len(),
            plural(resp.reversed.len()),
            resp.reversed.join(", ")
        );
    }
    match &args.output {
        Some(p) => {
            if same_file(p, &args.input) {
                return Err(format!("refusing to overwrite the input file {}", p.display()));
            }
            ResultWriter::write_text(p, &resp.fasta)?;
        }
        None => {
            use std::io::Write;
            let mut out = std::io::stdout().lock();
            out.write_all(resp.fasta.as_bytes()).map_err(|e| e.to_string())?;
        }
    }
    eprintln!(
        "Aligned {} sequences with MAFFT {} (pure Rust): {} columns in {:.2}s{}",
        resp.sequences,
        resp.strategy.as_deref().unwrap_or("?"),
        resp.columns,
        start.elapsed().as_secs_f32(),
        args.output.as_ref().map(|p| format!(" -> {}", p.display())).unwrap_or_default()
    );
    Ok(())
}

fn run_test(args: TestArgs) -> Result<(), String> {
    let settings = SettingsInput {
        primer_conc: args.primer_conc,
        salt_conc: args.salt_conc,
        min_temp: args.min_temp,
        max_temp: args.max_temp,
        min_gc: args.min_gc,
        max_gc: args.max_gc,
        ..SettingsInput::check_defaults()
    };
    let resp = api::check(&CheckRequest { primer: args.primer.clone(), partner: args.partner.clone(), settings });
    for w in &resp.warnings {
        eprintln!("warning: {w}");
    }
    if !resp.ok {
        return Err(resp.error.unwrap_or_else(|| "check failed".into()));
    }
    let a = resp.primer.as_ref().ok_or("internal error: no analysis")?;
    println!("{}", a.report());
    if let (Some(b), Some(pair)) = (&resp.partner, &resp.pair) {
        println!("Partner primer\n--------------\n{}", b.report());
        println!("Pair\n----");
        println!("ΔTm:                  {:.2} °C", pair.tm_diff);
        println!("Pair complementarity: any {:.2}, 3' {:.2}", pair.pair_any, pair.pair_end);
        println!("3'-3' dimer:          {} bp", pair.dimer_score);
        if let Some(d) = &pair.diagram {
            println!("\nStrongest 3' cross-dimer:\n{d}");
        }
    }
    if let Some(p) = &args.json {
        ResultWriter::write_primer_analysis_json(p, a)?;
        println!("JSON written to {}", p.display());
    }
    Ok(())
}

fn on_path(program: &str) -> bool {
    std::env::var_os("PATH").is_some_and(|paths| {
        std::env::split_paths(&paths).any(|dir| {
            let p = dir.join(program);
            p.is_file() || p.with_extension("exe").is_file()
        })
    })
}

fn run_doctor() -> Result<(), String> {
    println!("{}", degenprime::banner());
    println!("DeGenPrime doctor\n");
    let ok = |label: &str, detail: String| println!("  [ OK ] {label}: {detail}");
    let warn = |label: &str, detail: String| println!("  [WARN] {label}: {detail}");
    let fail = |label: &str, detail: String| println!("  [FAIL] {label}: {detail}");
    let mut failures = 0;

    ok(
        "version",
        format!(
            "{} ({} {}, {} thread{})",
            degenprime::VERSION,
            std::env::consts::OS,
            std::env::consts::ARCH,
            rayon::current_num_threads(),
            plural(rayon::current_num_threads())
        ),
    );

    let req = DesignRequest {
        input: degenprime::EXAMPLE_ALIGNMENT.into(),
        filename: Some("example_anfD.clust".into()),
        max_site_sequences: Some(0),
        ..Default::default()
    };
    let t = Instant::now();
    let first = api::design(&req);
    let ms = t.elapsed().as_secs_f64() * 1000.0;
    if first.ok && !first.pairs.is_empty() {
        ok(
            "core self-test",
            format!("{} primer pairs from the bundled anfD alignment in {ms:.0} ms", first.pairs.len()),
        );
    } else {
        failures += 1;
        fail("core self-test", first.error.clone().unwrap_or_else(|| "no primer pairs".into()));
    }
    if api::design(&req).csv == first.csv {
        ok("deterministic", "two runs gave identical results".into());
    } else {
        failures += 1;
        fail("deterministic", "two runs of the same input differed".into());
    }

    // The bundled anfD alignment was made with C MAFFT (L-INS-i); the
    // built-in aligner must reproduce it from the ungapped sequences.
    let t = Instant::now();
    match degenprime::msa::self_test() {
        Ok(strategy) => ok(
            "aligner",
            format!(
                "built-in MAFFT ({strategy}) reproduced the bundled MAFFT alignment exactly in {:.0} ms",
                t.elapsed().as_secs_f64() * 1000.0
            ),
        ),
        Err(e) => {
            failures += 1;
            fail("aligner", e);
        }
    }

    let tmp = std::env::temp_dir().join(format!("degenprime-doctor-{}.tmp", std::process::id()));
    match std::fs::write(&tmp, b"ok") {
        Ok(()) => {
            let _ = std::fs::remove_file(&tmp);
            ok("temp directory", format!("{} is writable", std::env::temp_dir().display()));
        }
        Err(e) => {
            failures += 1;
            fail("temp directory", e.to_string());
        }
    }

    #[cfg(feature = "gui")]
    match std::net::TcpListener::bind(("127.0.0.1", degenprime::gui::DEFAULT_PORT)) {
        Ok(_) => ok("desktop app", format!("port {} is free for `degenprime gui`", degenprime::gui::DEFAULT_PORT)),
        Err(_) => warn(
            "desktop app",
            format!("port {} is busy; `degenprime gui` will pick another", degenprime::gui::DEFAULT_PORT),
        ),
    }
    let opener = if cfg!(target_os = "macos") {
        Some("/usr/bin/open")
    } else if cfg!(target_os = "windows") {
        None
    } else {
        Some("xdg-open")
    };
    match opener {
        Some(p) if p.starts_with('/') && Path::new(p).exists() => ok("browser launcher", p.to_string()),
        Some(p) if on_path(p) => ok("browser launcher", p.to_string()),
        Some(p) => {
            warn("browser launcher", format!("{p} not found; `degenprime gui` will print the address to open instead"))
        }
        None => ok("browser launcher", "cmd /C start".into()),
    }

    if failures == 0 {
        println!("\nAll essential checks passed.");
        Ok(())
    } else {
        Err(format!("{failures} essential check{} failed", plural(failures)))
    }
}

/// Current UTC time as ISO-8601, without a date-time dependency.
fn chrono_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now().duration_since(UNIX_EPOCH).unwrap_or_default().as_secs();
    let days = secs / 86400;
    let (h, m, s) = ((secs % 86400) / 3600, (secs % 3600) / 60, secs % 60);
    let z = days as i64 + 719_468;
    let era = z.div_euclid(146_097);
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe as i64 + era * 400 + i64::from(mo <= 2);
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{m:02}:{s:02}Z")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_defaults_match_api_defaults() {
        let cli = Cli::try_parse_from(["degenprime", "design", "-i", "x.fa"]).unwrap();
        let Some(Command::Design(a)) = cli.command else { panic!("expected design") };
        let ours = serde_json::to_value(a.settings()).unwrap();
        let api = serde_json::to_value(SettingsInput::default()).unwrap();
        assert_eq!(ours, api);
        let d = DesignRequest::default();
        assert_eq!((a.align, a.align_mode, !a.keep_strand), (d.align, d.align_mode, d.adjust_direction));
    }

    #[test]
    fn output_paths_do_not_clobber_dots() {
        assert_eq!(out_path(Path::new("res.v2"), "csv"), PathBuf::from("res.v2.csv"));
        assert_eq!(out_path(Path::new("res.csv"), "json"), PathBuf::from("res.json"));
        assert_eq!(out_path(Path::new("res"), "aligned.fasta"), PathBuf::from("res.aligned.fasta"));
    }

    #[test]
    fn outputs_never_overwrite_the_input() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("genes.fasta");
        std::fs::write(&input, ">a\nACGT\n").unwrap();
        let base = input.with_extension("");
        let fasta = safe_out(&base, "fasta", &input);
        assert_ne!(fasta, input);
        assert_eq!(fasta.file_name().unwrap(), "genes_primers.fasta");
        assert_eq!(safe_out(&base, "csv", &input).file_name().unwrap(), "genes.csv");
    }

    #[test]
    fn check_alias_parses() {
        assert!(Cli::try_parse_from(["degenprime", "check", "-p", "ACGTACGTAC"]).is_ok());
    }

    #[test]
    fn align_options_parse() {
        let cli =
            Cli::try_parse_from(["degenprime", "design", "-i", "x.fa", "--align", "never", "--align_mode", "linsi"])
                .unwrap();
        let Some(Command::Design(a)) = cli.command else { panic!("expected design") };
        assert_eq!(a.align, AlignWhen::Never);
        assert_eq!(a.align_mode, AlignMode::Linsi);
        assert!(Cli::try_parse_from(["degenprime", "design", "-i", "x.fa", "--mafft"]).is_ok());
        assert!(Cli::try_parse_from(["degenprime", "align", "-i", "x.fa", "--mode", "fftns2"]).is_ok());
    }
}
