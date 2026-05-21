//! File I/O: read FASTA / CLUSTAL alignments, write CSV / JSON results.
//!
//! Bug fixes vs v2:
//!   - `read_clustal`: the v2 version checked `trimmed.starts_with(' ')` AFTER
//!     trimming the line, which is always false. Conservation lines were
//!     therefore parsed as sequence rows. Fixed by recognising conservation
//!     characters on the *raw* line.
//!   - Detection of CLUSTAL header is now lenient: any first line containing
//!     "CLUSTAL" (case-insensitive) qualifies.
//!   - We no longer reject otherwise-valid alignments with mixed line lengths;
//!     instead we pad short sequences with gaps to the longest length and
//!     warn the caller.

use crate::calculator::PrimerAnalysis;
use crate::models::{DataNode, DataSequence, PrimerPair, Sequence};
use crate::settings::Settings;
use crate::thermodynamics::Thermodynamics;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    Fasta,
    Clustal,
}

pub struct SequenceReader;

impl SequenceReader {
    pub fn detect_format(path: &Path) -> Result<FileFormat, String> {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_lowercase();
        if matches!(ext.as_str(), "clust" | "aln" | "clustal") {
            return Ok(FileFormat::Clustal);
        }
        if matches!(ext.as_str(), "fasta" | "fa" | "faa" | "fna" | "ffn" | "frn" | "fas") {
            return Ok(FileFormat::Fasta);
        }
        // Sniff first non-empty line.
        let f = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
        let mut reader = BufReader::new(f);
        let mut line = String::new();
        loop {
            line.clear();
            let n = reader.read_line(&mut line).map_err(|e| e.to_string())?;
            if n == 0 { return Err("empty file".into()); }
            let t = line.trim();
            if t.is_empty() { continue; }
            if t.starts_with('>') { return Ok(FileFormat::Fasta); }
            if t.to_ascii_uppercase().contains("CLUSTAL") { return Ok(FileFormat::Clustal); }
            // Otherwise assume FASTA (bare sequence files happen).
            return Ok(FileFormat::Fasta);
        }
    }

    pub fn read_fasta(path: &Path) -> Result<Vec<Sequence>, String> {
        let f = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
        let reader = BufReader::new(f);
        let mut out = Vec::new();
        let mut id = String::new();
        let mut seq = String::new();

        for line in reader.lines() {
            let line = line.map_err(|e| e.to_string())?;
            let t = line.trim();
            if t.is_empty() { continue; }
            if let Some(rest) = t.strip_prefix('>') {
                if !id.is_empty() && !seq.is_empty() {
                    out.push(Sequence::new(id.clone(), seq.clone()));
                }
                id = rest.trim().to_string();
                seq.clear();
            } else {
                seq.push_str(t);
            }
        }
        if !id.is_empty() && !seq.is_empty() {
            out.push(Sequence::new(id, seq));
        }
        if out.is_empty() {
            return Err("FASTA file contained no sequences".into());
        }
        Ok(out)
    }

    /// CLUSTAL parser. Skips header line, the conservation row, and blank
    /// lines. Concatenates each ID's sequence segments in file order.
    pub fn read_clustal(path: &Path) -> Result<Vec<Sequence>, String> {
        let f = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
        let reader = BufReader::new(f);

        // Preserve first-seen order while accumulating segments.
        let mut order: Vec<String> = Vec::new();
        let mut seqs: BTreeMap<String, String> = BTreeMap::new();
        let mut saw_header = false;

        for line in reader.lines() {
            let raw = line.map_err(|e| e.to_string())?;
            let trimmed = raw.trim_end();

            // Header line, only once.
            if !saw_header && trimmed.to_ascii_uppercase().contains("CLUSTAL") {
                saw_header = true;
                continue;
            }
            if !saw_header {
                // Some tools omit the header; accept anyway.
                saw_header = true;
            }

            // Blank line between blocks.
            if trimmed.trim().is_empty() { continue; }

            // Conservation lines: start with whitespace OR contain only
            // conservation symbols ('*', ':', '.', ' ').
            if raw.starts_with(' ') || raw.starts_with('\t') {
                continue;
            }
            let only_conservation = trimmed
                .chars()
                .all(|c| matches!(c, '*' | ':' | '.' | ' ' | '\t'));
            if only_conservation { continue; }

            // Sequence row: "ID  ACGT...  [count]"
            let mut parts = trimmed.split_whitespace();
            let id = match parts.next() {
                Some(x) => x.to_string(),
                None => continue,
            };
            let seg = match parts.next() {
                Some(x) => x.to_string(),
                None => continue,
            };

            if !seqs.contains_key(&id) {
                order.push(id.clone());
                seqs.insert(id.clone(), String::new());
            }
            if let Some(buf) = seqs.get_mut(&id) {
                buf.push_str(&seg);
            }
        }

        if order.is_empty() {
            return Err("CLUSTAL file contained no sequence rows".into());
        }
        Ok(order
            .into_iter()
            .map(|id| {
                let s = seqs.remove(&id).unwrap_or_default();
                Sequence::new(id, s)
            })
            .collect())
    }

    pub fn read(path: &Path) -> Result<Vec<Sequence>, String> {
        match Self::detect_format(path)? {
            FileFormat::Fasta => Self::read_fasta(path),
            FileFormat::Clustal => Self::read_clustal(path),
        }
    }

    /// True when every sequence is the same length.
    pub fn is_aligned(seqs: &[Sequence]) -> bool {
        if seqs.is_empty() { return true; }
        let n = seqs[0].normalize().len();
        seqs.iter().all(|s| s.normalize().len() == n)
    }

    /// Pad shorter sequences with gaps to the longest length (for tolerant
    /// handling of slightly ragged CLUSTAL output).
    pub fn pad_to_longest(seqs: &mut [Sequence]) {
        let max_len = seqs.iter().map(|s| s.normalize().len()).max().unwrap_or(0);
        for s in seqs {
            let norm = s.normalize();
            if norm.len() < max_len {
                let mut padded = norm;
                padded.extend(std::iter::repeat('-').take(max_len - padded.len()));
                s.sequence = padded;
            } else {
                s.sequence = norm;
            }
        }
    }

    /// Build a consensus DataSequence from an MSA. Each column becomes one
    /// DataNode (IUPAC code derived from observed bases).
    pub fn alignment_to_consensus(seqs: &[Sequence]) -> Result<DataSequence, String> {
        if seqs.is_empty() {
            return Err("no sequences provided".into());
        }
        let normed: Vec<String> = seqs.iter().map(|s| s.normalize()).collect();
        let len = normed[0].len();
        if !normed.iter().all(|s| s.len() == len) {
            return Err("sequences are not aligned (different lengths)".into());
        }

        let mut nodes = Vec::with_capacity(len);
        let cols: Vec<Vec<char>> = (0..len)
            .map(|c| normed.iter().filter_map(|s| s.chars().nth(c)).collect())
            .collect();
        for col in cols {
            nodes.push(DataNode::new(col));
        }
        Ok(DataSequence::new(nodes))
    }
}

// ----- Output ----------------------------------------------------------------

/// One row in the CSV / JSON output — fully self-describing, so the Python
/// GUI can read it without any cross-references.
#[derive(Debug, Clone, Serialize)]
pub struct PrimerRow {
    pub rank: usize,
    pub forward_primer: String,
    pub forward_position: usize,
    pub forward_length: usize,
    pub forward_tm: f32,
    pub forward_gc: f32,
    pub forward_degeneracy: usize,
    pub reverse_primer: String,
    pub reverse_position: usize, // on the forward strand
    pub reverse_length: usize,
    pub reverse_tm: f32,
    pub reverse_gc: f32,
    pub reverse_degeneracy: usize,
    pub amplicon_size: usize,
    pub tm_difference: f32,
    pub hairpin_score_fwd: usize,
    pub hairpin_score_rev: usize,
    pub dimer_score: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct DesignReport {
    pub version: String,
    pub generated_at: String,
    pub input_file: Option<String>,
    pub template_length: usize,
    pub settings: Settings,
    pub primer_pairs: Vec<PrimerRow>,
}

pub struct ResultWriter;

impl ResultWriter {
    /// Build the full rows from primer pairs (needed by both CSV and JSON output).
    pub fn build_rows(
        pairs: &[PrimerPair],
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
        _settings: &Settings,
    ) -> Vec<PrimerRow> {
        let template_len = fwd_seq.len();
        pairs
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let fwd_sub = fwd_seq.sub_seq(p.forward.index, p.forward.length);
                let rev_sub = rev_seq.sub_seq(p.reverse.index, p.reverse.length);
                let rev_5p = template_len
                    .saturating_sub(p.reverse.index)
                    .saturating_sub(p.reverse.length);
                PrimerRow {
                    rank: i + 1,
                    forward_primer: fwd_sub.codes(),
                    forward_position: p.forward.index,
                    forward_length: p.forward.length,
                    forward_tm: p.fwd_tm,
                    forward_gc: fwd_sub.gc_ratio() * 100.0,
                    forward_degeneracy: fwd_sub.degeneracy_count(),
                    reverse_primer: rev_sub.codes(),
                    reverse_position: rev_5p,
                    reverse_length: p.reverse.length,
                    reverse_tm: p.rev_tm,
                    reverse_gc: rev_sub.gc_ratio() * 100.0,
                    reverse_degeneracy: rev_sub.degeneracy_count(),
                    amplicon_size: p.amplicon_size,
                    tm_difference: p.tm_diff,
                    hairpin_score_fwd: Thermodynamics::hairpin_score(&fwd_sub),
                    hairpin_score_rev: Thermodynamics::hairpin_score(&rev_sub),
                    dimer_score: Thermodynamics::dimer_score(&fwd_sub, &rev_sub),
                }
            })
            .collect()
    }

    pub fn write_csv(path: &Path, rows: &[PrimerRow]) -> Result<(), String> {
        let mut f = File::create(path).map_err(|e| format!("create {}: {e}", path.display()))?;
        // Header
        writeln!(
            f,
            "rank,forward_primer,forward_position,forward_length,forward_tm,forward_gc,forward_degeneracy,reverse_primer,reverse_position,reverse_length,reverse_tm,reverse_gc,reverse_degeneracy,amplicon_size,tm_difference,hairpin_score_fwd,hairpin_score_rev,dimer_score"
        ).map_err(|e| e.to_string())?;
        for r in rows {
            writeln!(
                f,
                "{},{},{},{},{:.2},{:.1},{},{},{},{},{:.2},{:.1},{},{},{:.2},{},{},{}",
                r.rank,
                r.forward_primer,
                r.forward_position,
                r.forward_length,
                r.forward_tm,
                r.forward_gc,
                r.forward_degeneracy,
                r.reverse_primer,
                r.reverse_position,
                r.reverse_length,
                r.reverse_tm,
                r.reverse_gc,
                r.reverse_degeneracy,
                r.amplicon_size,
                r.tm_difference,
                r.hairpin_score_fwd,
                r.hairpin_score_rev,
                r.dimer_score,
            ).map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    pub fn write_json(path: &Path, report: &DesignReport) -> Result<(), String> {
        let s = serde_json::to_string_pretty(report).map_err(|e| e.to_string())?;
        let mut f = File::create(path).map_err(|e| format!("create {}: {e}", path.display()))?;
        f.write_all(s.as_bytes()).map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn write_detailed_report(
        path: &Path,
        rows: &[PrimerRow],
        settings: &Settings,
        template_len: usize,
    ) -> Result<(), String> {
        let mut f = File::create(path).map_err(|e| e.to_string())?;
        writeln!(f, "DeGenPrime detailed report").map_err(|e| e.to_string())?;
        writeln!(f, "===========================").map_err(|e| e.to_string())?;
        writeln!(f, "Template length: {template_len} bp").map_err(|e| e.to_string())?;
        writeln!(f, "Primer length:   {}-{} bp", settings.min_primer_len, settings.max_primer_len).map_err(|e| e.to_string())?;
        writeln!(f, "Tm range:        {:.1}-{:.1} °C", settings.min_tm, settings.max_tm).map_err(|e| e.to_string())?;
        writeln!(f, "GC range:        {:.0}-{:.0}%", settings.min_gc_content * 100.0, settings.max_gc_content * 100.0).map_err(|e| e.to_string())?;
        writeln!(f, "Primer conc:     {} nM", settings.primer_concentration).map_err(|e| e.to_string())?;
        writeln!(f, "Salt:            {} mM", settings.salt_concentration).map_err(|e| e.to_string())?;
        writeln!(f, "Amplicon ≥:      {} bp\n", settings.min_amplicon_len).map_err(|e| e.to_string())?;

        for r in rows {
            writeln!(f, "Pair {}", r.rank).map_err(|e| e.to_string())?;
            writeln!(f, "--------").map_err(|e| e.to_string())?;
            writeln!(f, "  Forward: {} ({} bp) @ pos {}  Tm={:.2}°C GC={:.1}% degen={}",
                r.forward_primer, r.forward_length, r.forward_position,
                r.forward_tm, r.forward_gc, r.forward_degeneracy).map_err(|e| e.to_string())?;
            writeln!(f, "  Reverse: {} ({} bp) @ pos {}  Tm={:.2}°C GC={:.1}% degen={}",
                r.reverse_primer, r.reverse_length, r.reverse_position,
                r.reverse_tm, r.reverse_gc, r.reverse_degeneracy).map_err(|e| e.to_string())?;
            writeln!(f, "  Amplicon: {} bp, ΔTm: {:.2}°C", r.amplicon_size, r.tm_difference).map_err(|e| e.to_string())?;
            writeln!(f, "  Secondary structure scores: fwd hairpin={} rev hairpin={} dimer={}",
                r.hairpin_score_fwd, r.hairpin_score_rev, r.dimer_score).map_err(|e| e.to_string())?;
            writeln!(f).map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    pub fn write_primer_analysis_json(
        path: &Path,
        analysis: &PrimerAnalysis,
    ) -> Result<(), String> {
        let s = serde_json::to_string_pretty(analysis).map_err(|e| e.to_string())?;
        let mut f = File::create(path).map_err(|e| e.to_string())?;
        f.write_all(s.as_bytes()).map_err(|e| e.to_string())?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp(content: &str, ext: &str) -> NamedTempFile {
        let mut f = tempfile::Builder::new()
            .suffix(&format!(".{ext}"))
            .tempfile()
            .unwrap();
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
        let content = "CLUSTAL W (1.83) multiple sequence alignment\n\
            \n\
            seqA   ACGTACGT\n\
            seqB   ACGTACGT\n\
                   ********\n\
            \n";
        let f = write_temp(content, "clust");
        let seqs = SequenceReader::read_clustal(f.path()).unwrap();
        assert_eq!(seqs.len(), 2, "got {:?}", seqs.iter().map(|s| &s.id).collect::<Vec<_>>());
        for s in &seqs {
            assert_eq!(s.sequence, "ACGTACGT");
        }
    }

    #[test]
    fn consensus_uses_iupac_for_mixed_columns() {
        let seqs = vec![
            Sequence::new("a".into(), "ACGT".into()),
            Sequence::new("b".into(), "ATGT".into()),
        ];
        let cons = SequenceReader::alignment_to_consensus(&seqs).unwrap();
        // Column 2: C and T → Y
        assert_eq!(cons.codes(), "AYGT");
    }
}
