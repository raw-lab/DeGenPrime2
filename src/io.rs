//! Reading alignments and writing results.
//!
//! Parsing works on strings so the CLI, the desktop GUI and the WebAssembly
//! build share one code path. Accepted input: FASTA (with or without headers),
//! CLUSTAL blocks written by Clustal, MAFFT or MUSCLE, and a bare sequence.
//! Whitespace and digits are ignored, `.` reads as a gap, `U` as `T` and `X`
//! as `N`; any other non-IUPAC character is reported instead of guessed at.

use crate::calculator::{Explain, PrimerAnalysis};
use crate::models::{iupac, ColumnCounts, DataSequence, PrimerPair, Sequence};
use crate::settings::Settings;
use crate::thermodynamics::Thermodynamics;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum FileFormat {
    Fasta,
    Clustal,
    Raw,
}

pub struct SequenceReader;

fn read_text(path: &Path) -> Result<String, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let text = String::from_utf8_lossy(&bytes).into_owned();
    Ok(text.trim_start_matches('\u{feff}').to_string())
}

fn is_clustal_header(line: &str) -> bool {
    let u = line.trim().to_ascii_uppercase();
    u.starts_with("CLUSTAL")
        || u.starts_with("MUSCLE")
        || u.starts_with("PROBCONS")
        || u.contains("MULTIPLE SEQUENCE ALIGNMENT")
}

fn looks_like_clustal_row(line: &str) -> bool {
    let t: Vec<&str> = line.split_whitespace().collect();
    (2..=3).contains(&t.len())
        && !t[0].chars().all(iupac::is_sequence_char)
        && t[1].chars().all(iupac::is_sequence_char)
}

/// Uppercase, drop whitespace/digits/stops, map `.`→`-`, `U`→`T`, `X`→`N`,
/// and reject anything that is not an IUPAC nucleotide or a gap.
pub fn clean_sequence(id: &str, raw: &str) -> Result<String, String> {
    let mut out = String::with_capacity(raw.len());
    let mut bad: Vec<char> = Vec::new();
    for c in raw.chars() {
        if c.is_whitespace() || c.is_ascii_digit() || c == '*' {
            continue;
        }
        let u = match c.to_ascii_uppercase() {
            '.' | '~' => '-',
            'U' => 'T',
            'X' => 'N',
            other => other,
        };
        if iupac::is_sequence_char(u) {
            out.push(u);
        } else if !bad.contains(&u) && bad.len() < 8 {
            bad.push(u);
        }
    }
    if !bad.is_empty() {
        let list: Vec<String> = bad.iter().map(|c| format!("'{c}'")).collect();
        return Err(format!(
            "Sequence '{id}' contains characters that are not DNA/IUPAC nucleotide codes ({}). DeGenPrime needs nucleotide sequences; protein alignments are not supported.",
            list.join(", ")
        ));
    }
    Ok(out)
}

fn push_record(out: &mut Vec<Sequence>, id: Option<String>, raw: &str) -> Result<(), String> {
    let id = match id {
        Some(s) if !s.is_empty() => s,
        _ => format!("sequence_{}", out.len() + 1),
    };
    let seq = clean_sequence(&id, raw)?;
    if !seq.is_empty() {
        out.push(Sequence::new(id, seq));
    }
    Ok(())
}

impl SequenceReader {
    /// Sniff the format of a file from its content.
    pub fn detect_format(path: &Path) -> Result<FileFormat, String> {
        let text = read_text(path)?;
        if text.trim().is_empty() {
            return Err("empty file".into());
        }
        Ok(Self::sniff(&text))
    }

    /// Sniff the format of text.
    pub fn sniff(text: &str) -> FileFormat {
        let first = match text.lines().map(str::trim).find(|l| !l.is_empty()) {
            Some(l) => l,
            None => return FileFormat::Raw,
        };
        if first.starts_with('>') {
            return FileFormat::Fasta;
        }
        if is_clustal_header(first) {
            return FileFormat::Clustal;
        }
        if text.lines().any(|l| l.trim_start().starts_with('>')) {
            return FileFormat::Fasta;
        }
        let rows: Vec<&str> =
            text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with(char::is_whitespace)).collect();
        if rows.len() >= 2 && rows.iter().filter(|l| looks_like_clustal_row(l)).count() * 2 > rows.len() {
            return FileFormat::Clustal;
        }
        FileFormat::Raw
    }

    /// Parse any supported text input.
    pub fn parse(text: &str) -> Result<(Vec<Sequence>, FileFormat), String> {
        let text = text.trim_start_matches('\u{feff}');
        // Classic Mac line endings (a lone CR) would make the whole file one
        // line; Windows CRLF is handled by `lines()` already.
        let owned;
        let text = if text.contains('\r') && !text.contains('\n') {
            owned = text.replace('\r', "\n");
            owned.as_str()
        } else {
            text
        };
        let format = Self::sniff(text);
        let seqs = match format {
            FileFormat::Fasta => Self::parse_fasta(text)?,
            FileFormat::Clustal => Self::parse_clustal(text)?,
            FileFormat::Raw => {
                let s = clean_sequence("sequence_1", text)?;
                if s.is_empty() {
                    return Err("No sequence found in the input.".into());
                }
                vec![Sequence::new("sequence_1".into(), s)]
            }
        };
        Ok((seqs, format))
    }

    /// FASTA. Records without a header (a bare sequence before the first `>`)
    /// or with an empty header get generated IDs instead of being dropped.
    pub fn parse_fasta(text: &str) -> Result<Vec<Sequence>, String> {
        let mut out = Vec::new();
        let mut id: Option<String> = None;
        let mut buf = String::new();
        let mut started = false;
        for line in text.lines() {
            let t = line.trim();
            if t.is_empty() || t.starts_with(';') {
                continue;
            }
            if let Some(rest) = t.strip_prefix('>') {
                if started {
                    push_record(&mut out, id.take(), &buf)?;
                }
                started = true;
                buf.clear();
                id = Some(rest.split_whitespace().next().unwrap_or("").to_string());
            } else {
                started = true;
                buf.push_str(t);
            }
        }
        if started {
            push_record(&mut out, id.take(), &buf)?;
        }
        if out.is_empty() {
            return Err("The FASTA input contained no sequences.".into());
        }
        Ok(out)
    }

    /// CLUSTAL-style blocks. Skips the header, conservation tracks (they start
    /// with whitespace) and blank lines. Rows are joined by their position in
    /// each block, so sequences whose (often truncated) names are the same
    /// stay separate; when blocks differ in size, rows are joined by name.
    pub fn parse_clustal(text: &str) -> Result<Vec<Sequence>, String> {
        let mut blocks: Vec<Vec<(String, String)>> = Vec::new();
        let mut current: Vec<(String, String)> = Vec::new();
        for raw in text.lines() {
            let line = raw.trim_end();
            let separator = line.trim().is_empty()
                || raw.starts_with(char::is_whitespace)
                || line.chars().all(|c| matches!(c, '*' | ':' | '.' | ' ' | '\t'));
            if separator {
                if !current.is_empty() {
                    blocks.push(std::mem::take(&mut current));
                }
                continue;
            }
            if is_clustal_header(line) {
                continue;
            }
            let mut parts = line.split_whitespace();
            let (id, seg) = match (parts.next(), parts.next()) {
                (Some(a), Some(b)) => (a, b),
                _ => continue,
            };
            if !seg.chars().all(|c| iupac::is_sequence_char(c) || c == '*' || c == '~') {
                continue;
            }
            current.push((id.to_string(), seg.to_string()));
        }
        if !current.is_empty() {
            blocks.push(current);
        }
        let Some(first) = blocks.first() else {
            return Err("The CLUSTAL input contained no sequence rows.".into());
        };
        let rows = first.len();
        let positional = blocks.iter().all(|b| b.len() == rows && b.iter().zip(first).all(|(x, y)| x.0 == y.0));
        let joined: Vec<(String, String)> = if positional {
            let mut out: Vec<(String, String)> = first.iter().map(|(id, _)| (id.clone(), String::new())).collect();
            for b in &blocks {
                for (k, (_, seg)) in b.iter().enumerate() {
                    out[k].1.push_str(seg);
                }
            }
            out
        } else {
            let mut order: Vec<String> = Vec::new();
            let mut map: HashMap<String, String> = HashMap::new();
            for (id, seg) in blocks.iter().flatten() {
                map.entry(id.clone())
                    .or_insert_with(|| {
                        order.push(id.clone());
                        String::new()
                    })
                    .push_str(seg);
            }
            order
                .into_iter()
                .map(|id| {
                    let seq = map.remove(&id).unwrap_or_default();
                    (id, seq)
                })
                .collect()
        };
        joined
            .into_iter()
            .map(|(id, raw)| {
                let s = clean_sequence(&id, &raw)?;
                Ok(Sequence::new(id, s))
            })
            .collect()
    }

    pub fn read_fasta(path: &Path) -> Result<Vec<Sequence>, String> {
        Self::parse_fasta(&read_text(path)?)
    }

    pub fn read_clustal(path: &Path) -> Result<Vec<Sequence>, String> {
        Self::parse_clustal(&read_text(path)?)
    }

    /// Read and parse a file, sniffing the format from its content.
    pub fn read(path: &Path) -> Result<Vec<Sequence>, String> {
        Ok(Self::parse(&read_text(path)?)?.0)
    }

    /// True when every sequence has the same length.
    pub fn is_aligned(seqs: &[Sequence]) -> bool {
        match seqs.first() {
            None => true,
            Some(s0) => {
                let n = residue_len(s0);
                seqs.iter().all(|s| residue_len(s) == n)
            }
        }
    }

    /// Pad shorter sequences with trailing gaps to the longest length.
    pub fn pad_to_longest(seqs: &mut [Sequence]) {
        let max_len = seqs.iter().map(residue_len).max().unwrap_or(0);
        for s in seqs {
            if !is_clean(&s.sequence) {
                s.sequence = s.normalize();
            }
            let n = s.sequence.len();
            if n < max_len {
                s.sequence.extend(std::iter::repeat('-').take(max_len - n));
            }
        }
    }

    /// Consensus of an alignment, one [`DataNode`] per column.
    ///
    /// The alignment is read row by row into one tally per column, in
    /// parallel over stretches of columns: no copy of the sequences, and
    /// memory-friendly for thousands of sequences. Each column sees the
    /// sequences in input order, so the result is the same as building every
    /// column on its own ([`crate::models::DataNode::from_column`]).
    pub fn alignment_to_consensus(seqs: &[Sequence]) -> Result<DataSequence, String> {
        if seqs.is_empty() {
            return Err("no sequences provided".into());
        }
        // Sequences built by hand may carry spaces or lower case; parsed ones never do.
        let cleaned: Vec<std::borrow::Cow<'_, str>> = seqs
            .iter()
            .map(|s| {
                if is_clean(&s.sequence) {
                    std::borrow::Cow::Borrowed(s.sequence.as_str())
                } else {
                    std::borrow::Cow::Owned(s.normalize())
                }
            })
            .collect();
        let len = cleaned[0].len();
        if cleaned.iter().any(|s| s.len() != len) {
            return Err("sequences are not aligned (different lengths)".into());
        }
        const CHUNK: usize = 4096;
        let mut cols = vec![ColumnCounts::default(); len];
        cols.par_chunks_mut(CHUNK).enumerate().for_each(|(ci, chunk)| {
            let start = ci * CHUNK;
            for s in &cleaned {
                let bytes = &s.as_bytes()[start..start + chunk.len()];
                for (acc, &b) in chunk.iter_mut().zip(bytes) {
                    acc.add(b);
                }
            }
        });
        Ok(DataSequence::new(cols.iter().map(ColumnCounts::node).collect()))
    }
}

/// Upper-case ASCII with no whitespace: what the parsers produce.
fn is_clean(s: &str) -> bool {
    s.bytes().all(|b| !b.is_ascii_whitespace() && !b.is_ascii_lowercase())
}

/// Length of a sequence as [`Sequence::normalize`] would leave it.
fn residue_len(s: &Sequence) -> usize {
    if s.sequence.bytes().any(|b| b.is_ascii_whitespace()) {
        s.normalize().len()
    } else {
        s.sequence.len()
    }
}

// ----- Output ----------------------------------------------------------------

/// One primer pair, fully self-describing (CSV and JSON share it).
///
/// Coordinates: `*_from` / `*_to` are 1-based, inclusive alignment columns of
/// the primer footprint on the forward strand (as in Primer3web and GenBank).
/// `forward_position` / `reverse_position` keep the 2.x meaning: 0-based
/// left-most column of each footprint.
#[derive(Debug, Clone, Serialize)]
pub struct PrimerRow {
    pub rank: usize,
    pub score: f32,
    pub forward_primer: String,
    pub forward_consensus: String,
    pub forward_position: usize,
    pub forward_from: usize,
    pub forward_to: usize,
    pub forward_length: usize,
    pub forward_tm: f32,
    pub forward_gc: f32,
    pub forward_degeneracy: usize,
    pub forward_fold: u64,
    pub forward_conservation: f32,
    pub forward_self_any: f32,
    pub forward_self_end: f32,
    pub reverse_primer: String,
    pub reverse_consensus: String,
    pub reverse_position: usize,
    pub reverse_from: usize,
    pub reverse_to: usize,
    pub reverse_length: usize,
    pub reverse_tm: f32,
    pub reverse_gc: f32,
    pub reverse_degeneracy: usize,
    pub reverse_fold: u64,
    pub reverse_conservation: f32,
    pub reverse_self_any: f32,
    pub reverse_self_end: f32,
    pub amplicon_size: usize,
    pub tm_difference: f32,
    pub pair_any: f32,
    pub pair_end: f32,
    pub hairpin_score_fwd: usize,
    pub hairpin_score_rev: usize,
    pub dimer_score: usize,
}

/// The JSON file written by `degenprime design`.
#[derive(Debug, Clone, Serialize)]
pub struct DesignReport {
    pub version: String,
    pub generated_at: String,
    pub input_file: Option<String>,
    pub template_length: usize,
    /// The built-in MAFFT strategy, when DeGenPrime aligned the input.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aligned_with: Option<String>,
    pub settings: Settings,
    pub primer_pairs: Vec<PrimerRow>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub explain: Option<Explain>,
}

pub const CSV_HEADER: &str = "rank,score,forward_primer,forward_consensus,forward_position,forward_from,forward_to,forward_length,forward_tm,forward_gc,forward_degeneracy,forward_fold,forward_conservation,forward_self_any,forward_self_end,reverse_primer,reverse_consensus,reverse_position,reverse_from,reverse_to,reverse_length,reverse_tm,reverse_gc,reverse_degeneracy,reverse_fold,reverse_conservation,reverse_self_any,reverse_self_end,amplicon_size,tm_difference,pair_any,pair_end,hairpin_score_fwd,hairpin_score_rev,dimer_score";

pub struct ResultWriter;

impl ResultWriter {
    pub fn build_rows(
        pairs: &[PrimerPair],
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
        _settings: &Settings,
    ) -> Vec<PrimerRow> {
        let n = fwd_seq.len();
        pairs
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let fw = fwd_seq.sub_seq(p.forward.index, p.forward.length);
                let rw = rev_seq.sub_seq(p.reverse.index, p.reverse.length);
                let r_left = p.reverse_fwd_strand_position(n);
                let fs = Thermodynamics::self_complementarity(&fw);
                let rs = Thermodynamics::self_complementarity(&rw);
                let pc = Thermodynamics::complementarity(&fw, &rw);
                PrimerRow {
                    rank: i + 1,
                    score: p.score,
                    forward_primer: fw.oligo(),
                    forward_consensus: fw.consensus_oligo(),
                    forward_position: p.forward.index,
                    forward_from: p.forward.index + 1,
                    forward_to: p.forward.index + p.forward.length,
                    forward_length: fw.actual_size(),
                    forward_tm: p.fwd_tm,
                    forward_gc: fw.gc_ratio() * 100.0,
                    forward_degeneracy: fw.degeneracy_count(),
                    forward_fold: fw.fold(),
                    forward_conservation: fw.mean_conservation() * 100.0,
                    forward_self_any: fs.any,
                    forward_self_end: fs.end,
                    reverse_primer: rw.oligo(),
                    reverse_consensus: rw.consensus_oligo(),
                    reverse_position: r_left,
                    reverse_from: r_left + 1,
                    reverse_to: r_left + p.reverse.length,
                    reverse_length: rw.actual_size(),
                    reverse_tm: p.rev_tm,
                    reverse_gc: rw.gc_ratio() * 100.0,
                    reverse_degeneracy: rw.degeneracy_count(),
                    reverse_fold: rw.fold(),
                    reverse_conservation: rw.mean_conservation() * 100.0,
                    reverse_self_any: rs.any,
                    reverse_self_end: rs.end,
                    amplicon_size: p.amplicon_size,
                    tm_difference: p.tm_diff,
                    pair_any: pc.any,
                    pair_end: pc.end,
                    hairpin_score_fwd: Thermodynamics::hairpin_score(&fw),
                    hairpin_score_rev: Thermodynamics::hairpin_score(&rw),
                    dimer_score: Thermodynamics::dimer_score(&fw, &rw),
                }
            })
            .collect()
    }

    pub fn to_csv(rows: &[PrimerRow]) -> String {
        let mut s = String::from(CSV_HEADER);
        s.push('\n');
        for r in rows {
            s.push_str(&format!(
                concat!(
                    "{},{:.2},{},{},{},{},{},{},{:.2},{:.1},{},{},{:.1},{:.2},{:.2},",
                    "{},{},{},{},{},{},{:.2},{:.1},{},{},{:.1},{:.2},{:.2},",
                    "{},{:.2},{:.2},{:.2},{},{},{}\n"
                ),
                r.rank,
                r.score,
                r.forward_primer,
                r.forward_consensus,
                r.forward_position,
                r.forward_from,
                r.forward_to,
                r.forward_length,
                r.forward_tm,
                r.forward_gc,
                r.forward_degeneracy,
                r.forward_fold,
                r.forward_conservation,
                r.forward_self_any,
                r.forward_self_end,
                r.reverse_primer,
                r.reverse_consensus,
                r.reverse_position,
                r.reverse_from,
                r.reverse_to,
                r.reverse_length,
                r.reverse_tm,
                r.reverse_gc,
                r.reverse_degeneracy,
                r.reverse_fold,
                r.reverse_conservation,
                r.reverse_self_any,
                r.reverse_self_end,
                r.amplicon_size,
                r.tm_difference,
                r.pair_any,
                r.pair_end,
                r.hairpin_score_fwd,
                r.hairpin_score_rev,
                r.dimer_score,
            ));
        }
        s
    }

    /// Primers as FASTA, ready to paste into an oligo order form.
    pub fn to_fasta(rows: &[PrimerRow], name: &str) -> String {
        let mut s = String::new();
        for r in rows {
            s.push_str(&format!(
                ">{name}_pair{}_F from={} to={} tm={:.1} product={}bp\n{}\n",
                r.rank, r.forward_from, r.forward_to, r.forward_tm, r.amplicon_size, r.forward_primer
            ));
            s.push_str(&format!(
                ">{name}_pair{}_R from={} to={} tm={:.1} product={}bp\n{}\n",
                r.rank, r.reverse_from, r.reverse_to, r.reverse_tm, r.amplicon_size, r.reverse_primer
            ));
        }
        s
    }

    /// Primer3-style plain-text report.
    pub fn report_text(
        rows: &[PrimerRow],
        settings: &Settings,
        template_len: usize,
        explain: Option<&Explain>,
        warnings: &[String],
    ) -> String {
        let mut s = String::new();
        let line = |s: &mut String, t: String| {
            s.push_str(&t);
            s.push('\n');
        };
        line(&mut s, format!("DeGenPrime {} primer design report", crate::VERSION));
        line(&mut s, "=".repeat(44));
        line(&mut s, format!("Alignment length:  {template_len} columns"));
        line(&mut s, format!("Primer length:     {}-{} nt", settings.min_primer_len, settings.max_primer_len));
        line(
            &mut s,
            format!(
                "Tm window:         {:.1}-{:.1} °C (max ΔTm {:.1} °C)",
                settings.min_tm, settings.max_tm, settings.max_tm_diff
            ),
        );
        line(
            &mut s,
            format!(
                "GC content:        {:.0}-{:.0}% (3' end ≤ {:.0}%)",
                settings.min_gc_content * 100.0,
                settings.max_gc_content * 100.0,
                settings.max_gc_end * 100.0
            ),
        );
        line(
            &mut s,
            format!(
                "Product size:      ≥ {} bp{}",
                settings.min_amplicon_len,
                settings.max_amplicon_len.map(|m| format!(", ≤ {m} bp")).unwrap_or_default()
            ),
        );
        if settings.begin_nucleotide.is_some() || settings.end_nucleotide.is_some() {
            line(
                &mut s,
                format!(
                    "Region:            columns {}-{}",
                    settings.begin_nucleotide.map(|b| b + 1).unwrap_or(1),
                    settings.end_nucleotide.map(|e| e.to_string()).unwrap_or_else(|| template_len.to_string())
                ),
            );
        }
        line(
            &mut s,
            format!(
                "Reaction:          {} nM primer, {} mM Na+",
                settings.primer_concentration, settings.salt_concentration
            ),
        );
        for w in warnings {
            line(&mut s, format!("WARNING: {w}"));
        }
        s.push('\n');
        for r in rows {
            line(
                &mut s,
                format!(
                    "PAIR {}   product {} bp   ΔTm {:.2} °C   score {:.2}",
                    r.rank, r.amplicon_size, r.tm_difference, r.score
                ),
            );
            line(
                &mut s,
                "  OLIGO      from    to  len     tm    gc%  degen   fold   any    3'   sequence (5'->3')".into(),
            );
            for (label, from, to, len, tm, gc, dg, fold, any, end, sq) in [
                (
                    "FORWARD",
                    r.forward_from,
                    r.forward_to,
                    r.forward_length,
                    r.forward_tm,
                    r.forward_gc,
                    r.forward_degeneracy,
                    r.forward_fold,
                    r.forward_self_any,
                    r.forward_self_end,
                    &r.forward_primer,
                ),
                (
                    "REVERSE",
                    r.reverse_from,
                    r.reverse_to,
                    r.reverse_length,
                    r.reverse_tm,
                    r.reverse_gc,
                    r.reverse_degeneracy,
                    r.reverse_fold,
                    r.reverse_self_any,
                    r.reverse_self_end,
                    &r.reverse_primer,
                ),
            ] {
                line(&mut s, format!(
                    "  {label:<8} {from:>6} {to:>5} {len:>4} {tm:>6.2} {gc:>6.1} {dg:>6} {fold:>6} {any:>5.2} {end:>5.2}   {sq}"
                ));
            }
            line(&mut s, format!("  Pair complementarity: any {:.2}, 3' {:.2}", r.pair_any, r.pair_end));
            s.push('\n');
        }
        if rows.is_empty() {
            line(&mut s, "No primer pairs passed the filters.".into());
            s.push('\n');
        }
        if let Some(e) = explain {
            line(&mut s, "STATISTICS".into());
            for (label, st) in [("Forward", &e.forward), ("Reverse", &e.reverse)] {
                line(&mut s, format!(
                    "  {label}: considered {}, too degenerate {}, gaps {}, GC {}, 3' GC {}, repeats {}, comp. ends {}, Tm low {}, Tm high {}, self-compl. {}, ok {}",
                    st.considered, st.degeneracy, st.deletions, st.gc_content, st.gc_3prime, st.repeats,
                    st.complementary_ends, st.tm_low, st.tm_high, st.self_complementarity, st.ok
                ));
            }
            let p = &e.pairs;
            line(&mut s, format!(
                "  Pairs: examined {}, overlapping {}, product too short {}, too long {}, ΔTm too high {}, pair compl. {}, near-duplicate {}, returned {}",
                p.examined, p.overlap, p.product_too_short, p.product_too_long, p.tm_diff, p.pair_complementarity, p.too_similar, p.ok
            ));
        }
        s
    }

    pub fn write_text(path: &Path, text: &str) -> Result<(), String> {
        std::fs::write(path, text).map_err(|e| format!("write {}: {e}", path.display()))
    }

    pub fn write_csv(path: &Path, rows: &[PrimerRow]) -> Result<(), String> {
        Self::write_text(path, &Self::to_csv(rows))
    }

    pub fn write_json(path: &Path, report: &DesignReport) -> Result<(), String> {
        let s = serde_json::to_string_pretty(report).map_err(|e| e.to_string())?;
        Self::write_text(path, &s)
    }

    pub fn write_detailed_report(
        path: &Path,
        rows: &[PrimerRow],
        settings: &Settings,
        template_len: usize,
    ) -> Result<(), String> {
        Self::write_text(path, &Self::report_text(rows, settings, template_len, None, &[]))
    }

    pub fn write_primer_analysis_json(path: &Path, analysis: &PrimerAnalysis) -> Result<(), String> {
        let s = serde_json::to_string_pretty(analysis).map_err(|e| e.to_string())?;
        Self::write_text(path, &s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp(content: &str, ext: &str) -> NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(&format!(".{ext}")).tempfile().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    #[test]
    fn read_fasta_basic() {
        let f = write_temp(">a\nACGT\n>b\nTGCA\n", "fasta");
        let seqs = SequenceReader::read_fasta(f.path()).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].id, "a");
        assert_eq!(seqs[0].sequence, "ACGT");
    }

    #[test]
    fn read_clustal_skips_conservation_line() {
        let content =
            "CLUSTAL W (1.83) multiple sequence alignment\n\nseqA   ACGTACGT\nseqB   ACGTACGT\n       ********\n\n";
        let f = write_temp(content, "clust");
        let seqs = SequenceReader::read_clustal(f.path()).unwrap();
        assert_eq!(seqs.len(), 2);
        assert!(seqs.iter().all(|s| s.sequence == "ACGTACGT"));
    }

    #[test]
    fn muscle_header_is_not_a_sequence() {
        let text = "MUSCLE (3.8) multiple sequence alignment\n\nseqA   ACGT-CGT\nseqB   ACGTACGT\n";
        let (seqs, fmt) = SequenceReader::parse(text).unwrap();
        assert_eq!(fmt, FileFormat::Clustal);
        assert_eq!(seqs.len(), 2);
    }

    #[test]
    fn consensus_uses_iupac_for_mixed_columns() {
        let seqs = vec![Sequence::new("a".into(), "ACGT".into()), Sequence::new("b".into(), "ATGT".into())];
        let cons = SequenceReader::alignment_to_consensus(&seqs).unwrap();
        assert_eq!(cons.codes(), "AYGT");
    }

    #[test]
    fn bare_sequence_is_accepted() {
        let (seqs, fmt) = SequenceReader::parse("  1 acgtacgtac gtacgtacgt\n 21 ACGUN\n").unwrap();
        assert_eq!(fmt, FileFormat::Raw);
        assert_eq!(seqs[0].sequence, "ACGTACGTACGTACGTACGTACGTN");
    }

    #[test]
    fn headerless_and_empty_header_fasta_get_ids() {
        let seqs = SequenceReader::parse_fasta("ACGTACGT\n>\nTTTTAAAA\n>named desc\nGGGG\n").unwrap();
        let ids: Vec<&str> = seqs.iter().map(|s| s.id.as_str()).collect();
        assert_eq!(ids, ["sequence_1", "sequence_2", "named"]);
    }

    #[test]
    fn protein_is_rejected_with_a_clear_message() {
        let err = SequenceReader::parse(">p\nMKLVEFQ\n").unwrap_err();
        assert!(err.contains("not DNA"), "{err}");
    }

    #[test]
    fn dots_are_gaps() {
        let (seqs, _) = SequenceReader::parse(">a\nAC..GT\n>b\nACGTGT\n").unwrap();
        assert_eq!(seqs[0].sequence, "AC--GT");
    }

    #[test]
    fn clustal_rows_with_the_same_name_stay_separate() {
        // MAFFT --clustalout truncates names; two sequences can end up with one.
        let text = "CLUSTAL\n\nNC_0455 ACGTAC\nNC_0455 ACGTTC\nMN9089 ACG-AC\n       ** * *\n\nNC_0455 GG\nNC_0455 GA\nMN9089 GG\n";
        let (seqs, fmt) = SequenceReader::parse(text).unwrap();
        assert_eq!(fmt, FileFormat::Clustal);
        let got: Vec<(&str, &str)> = seqs.iter().map(|s| (s.id.as_str(), s.sequence.as_str())).collect();
        assert_eq!(got, [("NC_0455", "ACGTACGG"), ("NC_0455", "ACGTTCGA"), ("MN9089", "ACG-ACGG")]);
    }

    #[test]
    fn classic_mac_line_endings() {
        let (seqs, _) = SequenceReader::parse(">a\rACGT\r>b\rACGA\r").unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[1].sequence, "ACGA");
    }

    #[test]
    fn row_wise_consensus_matches_column_by_column() {
        let mut x = 12345u64;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let alphabet = b"ACGTACGTACGT--RYSWKMN.acgt";
        let len = 9000; // spans several parallel chunks
        let seqs: Vec<Sequence> = (0..13)
            .map(|i| {
                let s: String = (0..len).map(|_| alphabet[(next() % alphabet.len() as u64) as usize] as char).collect();
                Sequence::new(format!("s{i}"), s)
            })
            .collect();
        let fast = SequenceReader::alignment_to_consensus(&seqs).unwrap();
        for c in (0..len).step_by(7) {
            let col: Vec<char> = seqs.iter().map(|s| s.sequence.as_bytes()[c] as char).collect();
            assert_eq!(fast.nodes()[c], crate::models::DataNode::from_column(col), "column {c}");
        }
    }

    #[test]
    fn csv_has_one_value_per_header_column() {
        let text = include_str!("../data/example_anfD.clust");
        let (seqs, _) = SequenceReader::parse(text).unwrap();
        let cons = SequenceReader::alignment_to_consensus(&seqs).unwrap();
        let rc = cons.reverse_complement();
        let settings = Settings::default();
        let out = crate::calculator::PrimerCalculator::design(&cons, &rc, &settings).unwrap();
        let rows = ResultWriter::build_rows(&out.pairs, &cons, &rc, &settings);
        let csv = ResultWriter::to_csv(&rows);
        let cols = CSV_HEADER.split(',').count();
        for l in csv.lines() {
            assert_eq!(l.split(',').count(), cols, "{l}");
        }
    }
}
