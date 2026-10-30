//! JSON API shared by every front end: the CLI, the desktop GUI (`POST
//! /api/<endpoint>`) and the WebAssembly build (`dgp_call`). Requests and
//! responses are plain JSON, so the web app does not care which engine runs it.
//!
//! Endpoints: `design`, `inspect`, `align`, `check`, `defaults`, `version`.

use crate::calculator::{Explain, PrimerAnalysis, PrimerCalculator};
use crate::io::{FileFormat, PrimerRow, ResultWriter, SequenceReader};
use crate::models::{DataNode, DataSequence, Sequence};
use crate::msa::{self, AlignMode, AlignOptions};
use crate::settings::{limits, RankBy, Settings};
use crate::thermodynamics::Thermodynamics;
use crate::VERSION;
use serde::{Deserialize, Serialize};

/// Which engine answered: `native` (CLI / desktop) or `wasm` (browser).
pub const ENGINE: &str = if cfg!(target_arch = "wasm32") { "wasm" } else { "native" };

/// The aligner built into every engine.
pub const ALIGNER: &str = "MAFFT 7.526 algorithms (rust-MAFFT 0.2.0, pure Rust)";

const UNALIGNED_WARNING: &str = "The sequences have different lengths, so they are not aligned. They were padded with gaps at the 3' end, which is rarely correct: let DeGenPrime align them (the default), or align them yourself first.";

/// Design parameters in user units (percent, °C, nM, mM) with the same names
/// as the CLI flags. Every field is optional in JSON; missing fields take the
/// documented defaults.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct SettingsInput {
    pub min_primer_len: usize,
    pub max_primer_len: usize,
    #[serde(alias = "amplicon")]
    pub min_amplicon: usize,
    pub max_amplicon: Option<usize>,
    /// Design region, 0-based alignment columns `[begin, end)`.
    pub begin: Option<usize>,
    pub end: Option<usize>,
    #[serde(alias = "min_tm")]
    pub min_temp: f32,
    #[serde(alias = "max_tm")]
    pub max_temp: f32,
    pub max_tm_diff: f32,
    pub min_gc: f32,
    pub max_gc: f32,
    pub max_gc_end: f32,
    pub max_degeneracy: f32,
    pub max_deletions: f32,
    pub max_poly_x: usize,
    pub max_self_any: f32,
    pub max_self_end: f32,
    pub max_pair_any: f32,
    pub max_pair_end: f32,
    pub primer_conc: f32,
    pub salt_conc: f32,
    pub max_primers: usize,
    pub min_3prime_distance: usize,
    pub rank_by: RankBy,
}

fn pct(f: f32) -> f32 {
    (f * 1000.0).round() / 10.0
}

fn sorted<T: PartialOrd + Copy>(a: T, b: T) -> (T, T) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

impl Default for SettingsInput {
    fn default() -> Self {
        let s = Settings::default();
        SettingsInput {
            min_primer_len: s.min_primer_len,
            max_primer_len: s.max_primer_len,
            min_amplicon: s.min_amplicon_len,
            max_amplicon: s.max_amplicon_len,
            begin: None,
            end: None,
            min_temp: s.min_tm,
            max_temp: s.max_tm,
            max_tm_diff: s.max_tm_diff,
            min_gc: pct(s.min_gc_content),
            max_gc: pct(s.max_gc_content),
            max_gc_end: pct(s.max_gc_end),
            max_degeneracy: pct(s.max_degeneracy_ratio),
            max_deletions: pct(s.max_deletion_ratio),
            max_poly_x: s.max_poly_x,
            max_self_any: s.max_self_any,
            max_self_end: s.max_self_end,
            max_pair_any: s.max_pair_any,
            max_pair_end: s.max_pair_end,
            primer_conc: s.primer_concentration,
            salt_conc: s.salt_concentration,
            max_primers: s.max_primers_output,
            min_3prime_distance: s.min_three_prime_distance,
            rank_by: s.rank_by,
        }
    }
}

impl SettingsInput {
    /// The looser pass/fail windows used when checking a single primer.
    pub fn check_defaults() -> Self {
        SettingsInput { min_temp: 50.0, max_temp: 65.0, min_gc: 30.0, max_gc: 70.0, ..Self::default() }
    }

    /// Convert to [`Settings`], returning a warning for every value that had
    /// to be clamped or swapped (2.1 adjusted these silently).
    pub fn to_settings(&self) -> (Settings, Vec<String>) {
        let mut w = Vec::new();
        let mut s = Settings::new()
            .with_primer_len(self.min_primer_len, self.max_primer_len)
            .with_amplicon_len(self.min_amplicon)
            .with_tm_range(self.min_temp, self.max_temp)
            .with_max_tm_diff(self.max_tm_diff)
            .with_gc_content(self.min_gc, self.max_gc)
            .with_max_gc_end(self.max_gc_end)
            .with_max_degeneracy(self.max_degeneracy)
            .with_max_deletions(self.max_deletions)
            .with_max_poly_x(self.max_poly_x)
            .with_complementarity(self.max_self_any, self.max_self_end, self.max_pair_any, self.max_pair_end)
            .with_primer_concentration(self.primer_conc)
            .with_salt_concentration(self.salt_conc)
            .with_max_primers(self.max_primers)
            .with_min_three_prime_distance(self.min_3prime_distance)
            .with_rank_by(self.rank_by)
            .with_region_bounds(self.begin, self.end);
        if let Some(m) = self.max_amplicon {
            s = s.with_max_amplicon_len(m);
        }

        let (lo, hi) = limits::PRIMER_LEN;
        let (a, b) = sorted(self.min_primer_len, self.max_primer_len);
        if (a, b) != (s.min_primer_len, s.max_primer_len) {
            w.push(format!(
                "Primer length {a}–{b} nt is outside the supported {lo}–{hi} nt; using {}–{} nt.",
                s.min_primer_len, s.max_primer_len
            ));
        } else if self.min_primer_len > self.max_primer_len {
            w.push("The minimum and maximum primer length were swapped.".into());
        }
        let (tlo, thi) = limits::TM;
        let (ta, tb) = sorted(self.min_temp, self.max_temp);
        if (ta - s.min_tm).abs() > 1e-4 || (tb - s.max_tm).abs() > 1e-4 {
            w.push(format!(
                "Tm window {ta:.1}–{tb:.1} °C is outside the supported {tlo:.0}–{thi:.0} °C; using {:.1}–{:.1} °C.",
                s.min_tm, s.max_tm
            ));
        } else if self.min_temp > self.max_temp {
            w.push("The minimum and maximum Tm were swapped.".into());
        }
        if self.min_gc > self.max_gc {
            w.push("The minimum and maximum GC content were swapped.".into());
        }
        for (name, v) in [
            ("GC content", self.min_gc),
            ("GC content", self.max_gc),
            ("3' end GC", self.max_gc_end),
            ("degenerate positions", self.max_degeneracy),
            ("gap columns", self.max_deletions),
        ] {
            if v.is_finite() && !(0.0..=100.0).contains(&v) {
                w.push(format!("{name} {v}% is outside 0–100%; using {}%.", v.clamp(0.0, 100.0)));
            }
        }
        if let (Some(b), Some(e)) = (self.begin, self.end) {
            if b > e {
                w.push("The start and end of the design region were swapped.".into());
            }
        }
        if self.max_tm_diff < 0.0 {
            w.push("A negative maximum ΔTm is not possible; using 0 °C.".into());
        }
        if s.max_primers_output != self.max_primers {
            w.push(format!(
                "The number of pairs must be {}–{}; using {}.",
                limits::MAX_PRIMERS.0,
                limits::MAX_PRIMERS.1,
                s.max_primers_output
            ));
        }
        if self.primer_conc < limits::MIN_PRIMER_CONC_NM {
            w.push(format!(
                "Primer concentration below {} nM is not supported; using {} nM.",
                limits::MIN_PRIMER_CONC_NM,
                s.primer_concentration
            ));
        }
        if self.salt_conc < limits::MIN_SALT_MM {
            w.push(format!(
                "Na+ below {} mM is not supported; using {} mM.",
                limits::MIN_SALT_MM,
                s.salt_concentration
            ));
        }
        (s, w)
    }
}

/// Whether to align the input before designing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum AlignWhen {
    /// Align when the sequences are not aligned (their lengths differ).
    #[default]
    Auto,
    /// Always (re)align, removing any gaps first.
    Always,
    /// Never align; pad unaligned sequences with gaps (as version 2.2 did).
    Never,
}

impl AlignWhen {
    pub fn parse(s: &str) -> Option<AlignWhen> {
        match s.trim().to_ascii_lowercase().as_str() {
            "auto" | "unaligned" => Some(AlignWhen::Auto),
            "always" | "realign" | "yes" | "true" => Some(AlignWhen::Always),
            "never" | "no" | "false" | "off" => Some(AlignWhen::Never),
            _ => None,
        }
    }
}

// ----- design ------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct DesignRequest {
    /// Alignment text: FASTA, CLUSTAL or a bare sequence.
    pub input: String,
    pub settings: SettingsInput,
    pub filename: Option<String>,
    /// Timestamp supplied by the caller (WebAssembly has no clock).
    pub generated_at: Option<String>,
    /// When to align the sequences with the built-in MAFFT.
    pub align: AlignWhen,
    /// MAFFT strategy.
    pub align_mode: AlignMode,
    /// Reverse-complement sequences on the opposite strand before aligning.
    pub adjust_direction: bool,
    /// Version 2.2: `true` asked for unaligned input to be aligned with an
    /// external MAFFT. That is now the default; accepted and ignored.
    pub mafft: bool,
    /// How many sequences to show per binding site (default 60, 0 = none).
    pub max_site_sequences: Option<usize>,
}

impl Default for DesignRequest {
    fn default() -> Self {
        DesignRequest {
            input: String::new(),
            settings: SettingsInput::default(),
            filename: None,
            generated_at: None,
            align: AlignWhen::Auto,
            align_mode: AlignMode::Auto,
            adjust_direction: true,
            mafft: false,
            max_site_sequences: None,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct InputSummary {
    pub format: FileFormat,
    pub sequences: usize,
    pub columns: usize,
    pub ids: Vec<String>,
    pub filename: Option<String>,
    /// Whether the input was already aligned (all sequences the same length).
    pub aligned_input: bool,
    /// The aligner and strategy, when DeGenPrime aligned the input.
    pub aligned_with: Option<String>,
    /// Sequences that were reverse-complemented before aligning.
    pub reversed: Vec<String>,
}

/// Per-column consensus for drawing: IUPAC codes, consensus bases and
/// conservation (0–100; for gap columns, the share of sequences with a gap).
#[derive(Debug, Clone, Serialize)]
pub struct ConsensusTrack {
    pub codes: String,
    pub bases: String,
    pub conservation: Vec<u8>,
}

impl ConsensusTrack {
    pub fn of(c: &DataSequence) -> Self {
        ConsensusTrack {
            codes: c.codes(),
            bases: c.consensus_string(),
            conservation: c.nodes().iter().map(|n| (n.ratio() * 100.0).round().clamp(0.0, 100.0) as u8).collect(),
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct SiteRow {
    pub id: String,
    pub seq: String,
}

/// The aligned sequences under one primer, read 5'→3' along the primer.
#[derive(Debug, Clone, Serialize)]
pub struct SiteBlock {
    pub primer: String,
    pub consensus: String,
    pub from: usize,
    pub to: usize,
    pub rows: Vec<SiteRow>,
    pub hidden: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct PairSites {
    pub rank: usize,
    pub forward: SiteBlock,
    pub reverse: SiteBlock,
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct DesignResponse {
    pub ok: bool,
    pub error: Option<String>,
    pub version: String,
    pub engine: String,
    pub warnings: Vec<String>,
    /// Informational messages (e.g. that the input was aligned).
    pub notes: Vec<String>,
    pub input: Option<InputSummary>,
    pub consensus: Option<ConsensusTrack>,
    pub settings: Option<Settings>,
    pub pairs: Vec<PrimerRow>,
    pub sites: Vec<PairSites>,
    pub explain: Option<Explain>,
    pub csv: String,
    pub report: String,
    pub fasta: String,
    /// The alignment primers were designed from, as FASTA, when DeGenPrime
    /// aligned the input itself.
    pub alignment: Option<String>,
}

fn site_block(
    seqs: &[Sequence],
    window: &DataSequence,
    from: usize,
    to: usize,
    reverse: bool,
    limit: usize,
) -> SiteBlock {
    let rows = seqs
        .iter()
        .take(limit)
        .map(|s| {
            let bytes = s.sequence.as_bytes();
            let (a, b) = (from.min(bytes.len()), to.min(bytes.len()));
            let slice: String = bytes[a..b].iter().map(|&c| c as char).collect();
            let seq = if reverse { slice.chars().rev().map(DataNode::complement_base).collect() } else { slice };
            SiteRow { id: s.id.clone(), seq }
        })
        .collect();
    SiteBlock {
        primer: window.codes(),
        consensus: window.consensus_string(),
        from: from + 1,
        to,
        rows,
        hidden: seqs.len().saturating_sub(limit),
    }
}

fn stem(name: &str) -> &str {
    let base = name.rsplit(['/', '\\']).next().unwrap_or(name);
    match base.split('.').next() {
        Some(s) if !s.is_empty() => s,
        _ => "degenprime",
    }
}

fn no_pairs_hint(e: &Explain) -> String {
    let mut msg = String::from("No primer pair met every constraint.");
    if e.forward.ok == 0 || e.reverse.ok == 0 {
        let (side, st) = if e.forward.ok == 0 { ("forward", &e.forward) } else { ("reverse", &e.reverse) };
        if st.considered == 0 {
            msg.push_str(" The design region is shorter than the primer length.");
        } else if let Some((why, n)) = st.top_reasons().first() {
            msg.push_str(&format!(
                " No {side} primer passed the filters; the most common reason was {why} ({n} of {} candidates).",
                st.considered
            ));
            msg.push_str(match *why {
                "Tm too low" => " Try lowering the minimum Tm or allowing longer primers.",
                "Tm too high" => " Try raising the maximum Tm or allowing shorter primers.",
                "GC content out of range" => " Try widening the GC range.",
                "too degenerate" => " Try allowing more degenerate positions (Advanced).",
                "self-complementary" => " Try relaxing the self-complementarity limits (Advanced).",
                "too many gaps" => " The alignment is gappy here; try another region or allow more gaps (Advanced).",
                _ => " Try relaxing that filter.",
            });
        }
    } else {
        let p = &e.pairs;
        if p.combinations > 0 && p.overlap == p.combinations {
            msg.push_str(" No reverse primer lies downstream of a forward primer: the conserved stretches are too close together. Try a longer alignment or region.");
            return msg;
        }
        // Overlapping combinations are not a useful explanation (half of all
        // combinations put the reverse primer upstream), so they are left out.
        let reasons = [
            ("the product was too short", p.product_too_short, " Lower the minimum product size."),
            ("the product was too long", p.product_too_long, " Raise the maximum product size."),
            ("the Tm difference was too large", p.tm_diff, " Allow a larger ΔTm or narrow the Tm window."),
            (
                "the primers could form a 3' dimer",
                p.pair_complementarity,
                " Relax the pair complementarity limits (Advanced).",
            ),
        ];
        if let Some((why, n, tip)) = reasons.iter().max_by_key(|r| r.1).filter(|r| r.1 > 0) {
            msg.push_str(&format!(
                " Primers passed on both strands, but most pairs failed because {why} ({n} pairs).{tip}"
            ));
        }
    }
    msg
}

/// Run a full design.
pub fn design(req: &DesignRequest) -> DesignResponse {
    let mut resp = DesignResponse { version: VERSION.into(), engine: ENGINE.into(), ..Default::default() };
    match design_into(req, &mut resp) {
        Ok(()) => resp.ok = true,
        Err(e) => {
            resp.ok = false;
            resp.error = Some(e);
        }
    }
    resp
}

/// Align `seqs` when `req` asks for it; returns (aligner label, reversed IDs).
fn align_if_needed(
    req: &DesignRequest,
    seqs: &mut Vec<Sequence>,
    resp_warnings: &mut Vec<String>,
    notes: &mut Vec<String>,
) -> (Option<String>, Vec<String>) {
    let aligned = SequenceReader::is_aligned(seqs);
    let wanted = seqs.len() > 1
        && match req.align {
            AlignWhen::Never => false,
            AlignWhen::Always => true,
            AlignWhen::Auto => !aligned,
        };
    if !wanted {
        return (None, Vec::new());
    }
    let opts = AlignOptions { mode: req.align_mode, adjust_direction: req.adjust_direction };
    match msa::align(seqs, &opts) {
        Ok(out) => {
            let label = format!("MAFFT {}", out.strategy);
            notes.push(if aligned {
                format!("The sequences were realigned with {label} (pure-Rust MAFFT) before designing.")
            } else {
                format!(
                    "The sequences had different lengths, so DeGenPrime aligned them with {label} (pure-Rust MAFFT) before designing. Download the alignment to check it."
                )
            });
            if !out.reversed.is_empty() {
                let list: Vec<&str> = out.reversed.iter().take(10).map(String::as_str).collect();
                resp_warnings.push(format!(
                    "{} sequence{} lay on the opposite strand and {} reverse-complemented before aligning: {}{}.",
                    out.reversed.len(),
                    if out.reversed.len() == 1 { "" } else { "s" },
                    if out.reversed.len() == 1 { "was" } else { "were" },
                    list.join(", "),
                    if out.reversed.len() > 10 { ", …" } else { "" }
                ));
            }
            *seqs = out.sequences;
            (Some(label), out.reversed)
        }
        Err(e) => {
            resp_warnings.push(format!("Alignment failed: {e}"));
            (None, Vec::new())
        }
    }
}

fn design_into(req: &DesignRequest, resp: &mut DesignResponse) -> Result<(), String> {
    let (settings, warnings) = req.settings.to_settings();
    resp.warnings.extend(warnings);
    settings.validate()?;
    if req.input.trim().is_empty() {
        return Err("No sequences yet. Paste an alignment or choose a file.".into());
    }
    let (mut seqs, format) = SequenceReader::parse(&req.input)?;
    let aligned_input = SequenceReader::is_aligned(&seqs);
    let mut notes = Vec::new();
    let (aligned_with, reversed) = align_if_needed(req, &mut seqs, &mut resp.warnings, &mut notes);
    if aligned_with.is_some() {
        resp.alignment = Some(msa::to_fasta(&seqs));
    }
    if seqs.len() > 1 && !SequenceReader::is_aligned(&seqs) {
        resp.warnings.push(UNALIGNED_WARNING.into());
        SequenceReader::pad_to_longest(&mut seqs);
    }
    let consensus = SequenceReader::alignment_to_consensus(&seqs)?;
    let n = consensus.len();
    if n < settings.min_primer_len * 2 {
        return Err(format!(
            "The alignment is only {n} columns long, too short for two {}-nt primers.",
            settings.min_primer_len
        ));
    }
    if let Some(b) = settings.begin_nucleotide {
        if b >= n {
            return Err(format!("The region starts at column {} but the alignment has only {n} columns.", b + 1));
        }
    }
    if let Some(e) = settings.end_nucleotide {
        if e > n {
            resp.warnings.push(format!("The region end ({e}) is past the end of the alignment; using {n}."));
        }
    }
    if seqs.len() == 1 {
        resp.warnings.push("Only one sequence was given, so the primers cannot be degenerate. Paste an alignment of related sequences to design degenerate primers.".into());
    }
    let rc = consensus.reverse_complement();
    let outcome = PrimerCalculator::design(&consensus, &rc, &settings)?;
    let rows = ResultWriter::build_rows(&outcome.pairs, &consensus, &rc, &settings);
    let limit = req.max_site_sequences.unwrap_or(60);
    if limit > 0 {
        resp.sites = rows
            .iter()
            .zip(&outcome.pairs)
            .map(|(row, p)| PairSites {
                rank: row.rank,
                forward: site_block(
                    &seqs,
                    &consensus.sub_seq(p.forward.index, p.forward.length),
                    p.forward.index,
                    p.forward.index + p.forward.length,
                    false,
                    limit,
                ),
                reverse: site_block(
                    &seqs,
                    &rc.sub_seq(p.reverse.index, p.reverse.length),
                    row.reverse_position,
                    row.reverse_position + p.reverse.length,
                    true,
                    limit,
                ),
            })
            .collect();
    }
    if rows.is_empty() {
        resp.warnings.push(no_pairs_hint(&outcome.explain));
    }
    let name = req.filename.as_deref().map(stem).unwrap_or("degenprime");
    resp.csv = ResultWriter::to_csv(&rows);
    resp.report = ResultWriter::report_text(&rows, &settings, n, Some(&outcome.explain), &resp.warnings);
    resp.fasta = ResultWriter::to_fasta(&rows, name);
    resp.input = Some(InputSummary {
        format,
        sequences: seqs.len(),
        columns: n,
        ids: seqs.iter().take(500).map(|s| s.id.clone()).collect(),
        filename: req.filename.clone(),
        aligned_input,
        aligned_with,
        reversed,
    });
    resp.notes = notes;
    resp.consensus = Some(ConsensusTrack::of(&consensus));
    resp.settings = Some(settings);
    resp.pairs = rows;
    resp.explain = Some(outcome.explain);
    Ok(())
}

// ----- inspect (live preview while typing) ---------------------------------------

#[derive(Debug, Clone, Deserialize, Default)]
#[serde(default)]
pub struct InspectRequest {
    pub input: String,
    pub filename: Option<String>,
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct InspectResponse {
    pub ok: bool,
    pub error: Option<String>,
    pub warnings: Vec<String>,
    pub notes: Vec<String>,
    pub input: Option<InputSummary>,
    /// `None` for unaligned input: its consensus would be meaningless.
    pub consensus: Option<ConsensusTrack>,
    /// The strategy `mafft --auto` would pick for these sequences.
    pub auto_strategy: Option<String>,
}

pub fn inspect(req: &InspectRequest) -> InspectResponse {
    let mut resp = InspectResponse::default();
    if req.input.trim().is_empty() {
        resp.ok = true;
        return resp;
    }
    let (seqs, format) = match SequenceReader::parse(&req.input) {
        Ok(x) => x,
        Err(e) => {
            resp.error = Some(e);
            return resp;
        }
    };
    let aligned = SequenceReader::is_aligned(&seqs);
    let longest = seqs.iter().map(|s| s.sequence.bytes().filter(|&c| c != b'-').count()).max().unwrap_or(0);
    resp.auto_strategy = Some(msa::auto_strategy(seqs.len(), longest).to_string());
    let summary = |columns: usize| InputSummary {
        format,
        sequences: seqs.len(),
        columns,
        ids: seqs.iter().take(500).map(|s| s.id.clone()).collect(),
        filename: req.filename.clone(),
        aligned_input: aligned,
        aligned_with: None,
        reversed: Vec::new(),
    };
    if seqs.len() > 1 && !aligned {
        resp.notes.push(format!(
            "These {} sequences are not aligned (their lengths differ, up to {longest} nt). DeGenPrime will align them with MAFFT ({}) when you design primers.",
            seqs.len(),
            msa::auto_strategy(seqs.len(), longest)
        ));
        resp.input = Some(summary(longest));
        resp.ok = true;
        return resp;
    }
    match SequenceReader::alignment_to_consensus(&seqs) {
        Ok(c) => {
            resp.input = Some(summary(c.len()));
            resp.consensus = Some(ConsensusTrack::of(&c));
            resp.ok = true;
        }
        Err(e) => resp.error = Some(e),
    }
    resp
}

// ----- align -----------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct AlignRequest {
    pub input: String,
    pub mode: AlignMode,
    pub adjust_direction: bool,
}

impl Default for AlignRequest {
    fn default() -> Self {
        AlignRequest { input: String::new(), mode: AlignMode::Auto, adjust_direction: true }
    }
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct AlignResponse {
    pub ok: bool,
    pub error: Option<String>,
    pub version: String,
    pub strategy: Option<String>,
    pub sequences: usize,
    pub columns: usize,
    pub reversed: Vec<String>,
    pub fasta: String,
}

/// Align sequences with the built-in MAFFT (gaps in the input are removed).
pub fn align(req: &AlignRequest) -> AlignResponse {
    let mut resp = AlignResponse { version: VERSION.into(), ..Default::default() };
    let seqs = match SequenceReader::parse(&req.input) {
        Ok((s, _)) => s,
        Err(e) => {
            resp.error = Some(e);
            return resp;
        }
    };
    match msa::align(&seqs, &AlignOptions { mode: req.mode, adjust_direction: req.adjust_direction }) {
        Ok(out) => {
            resp.ok = true;
            resp.strategy = Some(out.strategy.to_string());
            resp.sequences = out.sequences.len();
            resp.columns = out.sequences.first().map_or(0, |s| s.sequence.len());
            resp.fasta = msa::to_fasta(&out.sequences);
            resp.reversed = out.reversed;
        }
        Err(e) => resp.error = Some(e),
    }
    resp
}

// ----- check a primer ---------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct CheckRequest {
    pub primer: String,
    pub partner: Option<String>,
    pub settings: SettingsInput,
}

impl Default for CheckRequest {
    fn default() -> Self {
        CheckRequest { primer: String::new(), partner: None, settings: SettingsInput::check_defaults() }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct PairCheck {
    pub tm_diff: f32,
    pub pair_any: f32,
    pub pair_end: f32,
    pub dimer_score: usize,
    pub diagram: Option<String>,
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct CheckResponse {
    pub ok: bool,
    pub error: Option<String>,
    pub version: String,
    pub warnings: Vec<String>,
    pub primer: Option<PrimerAnalysis>,
    pub partner: Option<PrimerAnalysis>,
    pub pair: Option<PairCheck>,
}

/// Normalise a typed primer: drop 5'/3' labels, spaces and hyphens.
pub fn clean_primer(raw: &str, label: &str) -> Result<String, String> {
    let mut s = raw.trim().to_string();
    for p in ["5'-", "5′-", "5'", "5′"] {
        if let Some(r) = s.strip_prefix(p) {
            s = r.to_string();
            break;
        }
    }
    for p in ["-3'", "-3′", "3'", "3′"] {
        if let Some(r) = s.strip_suffix(p) {
            s = r.to_string();
            break;
        }
    }
    let cleaned: String = s
        .chars()
        .filter(|c| !c.is_whitespace() && *c != '-')
        .map(|c| match c.to_ascii_uppercase() {
            'U' => 'T',
            other => other,
        })
        .collect();
    if cleaned.is_empty() {
        return Err(format!("{label}: enter a sequence, 5'→3'."));
    }
    let mut bad: Vec<char> = Vec::new();
    for c in cleaned.chars() {
        if crate::models::iupac::mask(c).is_none() && !bad.contains(&c) {
            bad.push(c);
        }
    }
    if !bad.is_empty() {
        let list: Vec<String> = bad.iter().take(8).map(|c| format!("'{c}'")).collect();
        return Err(format!("{label} contains characters that are not IUPAC nucleotide codes: {}.", list.join(", ")));
    }
    if cleaned.len() < 5 {
        return Err(format!("{label} is too short to analyse (at least 5 nt)."));
    }
    if cleaned.len() > 200 {
        return Err(format!("{label} is longer than 200 nt; this page checks primers, not templates."));
    }
    Ok(cleaned)
}

pub fn check(req: &CheckRequest) -> CheckResponse {
    let mut resp = CheckResponse { version: VERSION.into(), ..Default::default() };
    let (settings, w) = req.settings.to_settings();
    resp.warnings = w;
    let a = match clean_primer(&req.primer, "Primer") {
        Ok(s) => s,
        Err(e) => {
            resp.error = Some(e);
            return resp;
        }
    };
    let a_seq = DataSequence::from_string(&a);
    let a_an = PrimerCalculator::analyze_primer(&a_seq, &settings);
    if let Some(p) = req.partner.as_deref().filter(|p| !p.trim().is_empty()) {
        let b = match clean_primer(p, "Partner primer") {
            Ok(s) => s,
            Err(e) => {
                resp.error = Some(e);
                return resp;
            }
        };
        let b_seq = DataSequence::from_string(&b);
        let b_an = PrimerCalculator::analyze_primer(&b_seq, &settings);
        let c = Thermodynamics::complementarity(&a_seq, &b_seq);
        resp.pair = Some(PairCheck {
            tm_diff: (a_an.tm - b_an.tm).abs(),
            pair_any: c.any,
            pair_end: c.end,
            dimer_score: Thermodynamics::dimer_score(&a_seq, &b_seq),
            diagram: if c.end > 0.0 { Thermodynamics::duplex_diagram(&a_seq, &b_seq, &c) } else { None },
        });
        resp.partner = Some(b_an);
    }
    resp.primer = Some(a_an);
    resp.ok = true;
    resp
}

// ----- dispatcher ---------------------------------------------------------------------

#[derive(Serialize)]
struct Defaults {
    design: SettingsInput,
    check: SettingsInput,
}

#[derive(Serialize)]
pub struct VersionInfo {
    pub version: String,
    pub engine: String,
    /// The built-in aligner (version 2.2 reported an external MAFFT here).
    pub mafft: Option<String>,
    pub aligner: String,
}

#[derive(Serialize)]
struct ErrorResponse {
    ok: bool,
    error: String,
    version: String,
}

fn to_json<T: Serialize>(v: &T) -> String {
    serde_json::to_string(v).unwrap_or_else(|e| format!("{{\"ok\":false,\"error\":\"serialization failed: {e}\"}}"))
}

fn parse<T: serde::de::DeserializeOwned>(body: &str) -> Result<T, String> {
    serde_json::from_str(body).map_err(|e| format!("invalid request: {e}"))
}

/// Answer one JSON request. Never panics on bad input.
pub fn handle(endpoint: &str, body: &str) -> String {
    let body = if body.trim().is_empty() { "{}" } else { body };
    let result: Result<String, String> = match endpoint.trim_matches('/') {
        "design" => parse::<DesignRequest>(body).map(|r| to_json(&design(&r))),
        "inspect" => parse::<InspectRequest>(body).map(|r| to_json(&inspect(&r))),
        "align" => parse::<AlignRequest>(body).map(|r| to_json(&align(&r))),
        "check" => parse::<CheckRequest>(body).map(|r| to_json(&check(&r))),
        "defaults" => {
            Ok(to_json(&Defaults { design: SettingsInput::default(), check: SettingsInput::check_defaults() }))
        }
        "version" => Ok(to_json(&VersionInfo {
            version: VERSION.into(),
            engine: ENGINE.into(),
            mafft: Some(ALIGNER.into()),
            aligner: ALIGNER.into(),
        })),
        other => Err(format!("unknown endpoint '{other}'")),
    };
    result.unwrap_or_else(|e| to_json(&ErrorResponse { ok: false, error: e, version: VERSION.into() }))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EXAMPLE: &str = include_str!("../data/example_anfD.clust");

    #[test]
    fn design_example_end_to_end() {
        let req = DesignRequest { input: EXAMPLE.into(), ..Default::default() };
        let r = design(&req);
        assert!(r.ok, "{:?}", r.error);
        assert!(!r.pairs.is_empty() && r.pairs.len() <= 5);
        assert_eq!(r.csv.lines().count(), r.pairs.len() + 1);
        assert_eq!(r.sites.len(), r.pairs.len());
        assert_eq!(r.sites[0].forward.rows.len(), 5);
        assert_eq!(r.sites[0].forward.rows[0].seq.len(), r.sites[0].forward.primer.len());
        let c = r.consensus.unwrap();
        assert_eq!(c.codes.len(), c.conservation.len());
    }

    #[test]
    fn clamps_are_reported() {
        let s = SettingsInput { min_temp: 45.0, max_primers: 500, ..Default::default() };
        let (settings, w) = s.to_settings();
        assert_eq!(settings.min_tm, 50.0);
        assert_eq!(w.len(), 2, "{w:?}");
    }

    #[test]
    fn check_handles_labels_and_partner() {
        let req = CheckRequest {
            primer: "5'-GCTACCTGAACACCATCCCG-3'".into(),
            partner: Some("cccTGGATGTTGTASTCGCC".into()),
            ..Default::default()
        };
        let r = check(&req);
        assert!(r.ok, "{:?}", r.error);
        assert_eq!(r.primer.as_ref().unwrap().sequence, "GCTACCTGAACACCATCCCG");
        assert_eq!(r.partner.as_ref().unwrap().fold, 2);
        assert!(r.pair.is_some());
    }

    #[test]
    fn dispatcher_rejects_bad_requests_gracefully() {
        assert!(handle("nope", "{}").contains("unknown endpoint"));
        assert!(handle("design", "{not json").contains("invalid request"));
        assert!(handle("design", "{}").contains("No sequences yet"));
        assert!(handle("defaults", "").contains("min_primer_len"));
    }

    fn ungapped_example() -> String {
        let (seqs, _) = SequenceReader::parse(EXAMPLE).unwrap();
        seqs.iter().map(|s| format!(">{}\n{}\n", s.id, s.sequence.replace('-', ""))).collect()
    }

    #[test]
    fn unaligned_input_is_aligned_with_the_built_in_mafft() {
        let r = design(&DesignRequest { input: ungapped_example(), ..Default::default() });
        assert!(r.ok, "{:?}", r.error);
        let inp = r.input.as_ref().unwrap();
        assert!(!inp.aligned_input);
        assert_eq!(inp.aligned_with.as_deref(), Some("MAFFT L-INS-i"));
        assert!(r.notes.iter().any(|n| n.contains("aligned them")));
        // L-INS-i reproduces the bundled alignment, so the primers are the same.
        let aligned = design(&DesignRequest { input: EXAMPLE.into(), ..Default::default() });
        assert_eq!(r.csv, aligned.csv);
        assert!(r.alignment.as_ref().unwrap().starts_with('>'));
        assert!(aligned.alignment.is_none());
    }

    #[test]
    fn alignment_can_be_turned_off() {
        let r = design(&DesignRequest { input: ungapped_example(), align: AlignWhen::Never, ..Default::default() });
        assert!(r.warnings.iter().any(|w| w.contains("not aligned")));
        assert!(r.input.map_or(true, |i| i.aligned_with.is_none()));
    }

    #[test]
    fn inspect_does_not_draw_unaligned_input() {
        let r = inspect(&InspectRequest { input: ungapped_example(), filename: None });
        assert!(r.ok);
        assert!(r.consensus.is_none());
        assert_eq!(r.auto_strategy.as_deref(), Some("L-INS-i"));
        assert!(r.notes[0].contains("not aligned"));
    }

    #[test]
    fn align_endpoint() {
        let out = handle("align", &serde_json::json!({ "input": ungapped_example(), "mode": "fftns2" }).to_string());
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["ok"], true);
        assert_eq!(v["strategy"], "FFT-NS-2");
        assert_eq!(v["sequences"], 5);
    }

    #[test]
    fn non_finite_and_out_of_range_settings_are_reported() {
        let s = SettingsInput { min_gc: -5.0, begin: Some(500), end: Some(100), ..Default::default() };
        let (settings, w) = s.to_settings();
        assert!(w.iter().any(|x| x.contains("outside 0–100")), "{w:?}");
        assert!(w.iter().any(|x| x.contains("swapped")), "{w:?}");
        assert_eq!((settings.begin_nucleotide, settings.end_nucleotide), (Some(100), Some(500)));
        let bad = SettingsInput { min_temp: f32::NAN, ..Default::default() };
        let r = design(&DesignRequest { input: EXAMPLE.into(), settings: bad, ..Default::default() });
        assert!(!r.ok);
    }

    #[test]
    fn region_past_the_end_is_an_error_not_a_panic() {
        let s = SettingsInput { begin: Some(99_999), ..Default::default() };
        let r = design(&DesignRequest { input: EXAMPLE.into(), settings: s, ..Default::default() });
        assert!(!r.ok);
    }
}
