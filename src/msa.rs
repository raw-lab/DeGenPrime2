//! Aligning unaligned input with MAFFT, in pure Rust.
//!
//! DeGenPrime designs primers from a multiple alignment. When the sequences
//! it is given have different lengths they are not aligned, and this module
//! aligns them with the vendored rust-MAFFT port ([`crate::mafft`]), which
//! reproduces C MAFFT's alignments. It needs no external program, so it works
//! the same in the command-line tool, the desktop app and the website
//! (WebAssembly).
//!
//! [`AlignMode::Auto`] follows `mafft --auto`: L-INS-i for small sets,
//! FFT-NS-i for medium ones and FFT-NS-2 for large ones. Sequences that lie
//! on the opposite strand can be reverse-complemented first, like
//! `mafft --adjustdirection`.

use crate::mafft::pipeline::adjust_direction::{adjust_direction_mode, AdjustMode};
use crate::mafft::pipeline::{AlignmentMode, MafftEngine};
use crate::mafft::types::{SeqType, Sequence as MafftSeq, SequenceSet};
use crate::models::Sequence;
use serde::{Deserialize, Serialize};

/// The MAFFT strategy to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum AlignMode {
    /// `mafft --auto`: picks a strategy from the number and length of the
    /// sequences.
    #[default]
    Auto,
    /// FFT-NS-2: fast progressive alignment (`mafft --retree 2`).
    Fftns2,
    /// FFT-NS-i: progressive alignment plus iterative refinement.
    Fftnsi,
    /// L-INS-i: most accurate for up to ~200 sequences (`mafft --localpair`).
    Linsi,
    /// G-INS-i: for sequences that align over their whole length.
    Ginsi,
    /// E-INS-i: for sequences with long unalignable insertions.
    Einsi,
}

impl AlignMode {
    pub fn parse(s: &str) -> Option<AlignMode> {
        let k: String = s.trim().to_ascii_lowercase().chars().filter(|c| c.is_ascii_alphanumeric()).collect();
        Some(match k.as_str() {
            "auto" | "mafft" => AlignMode::Auto,
            "fftns2" | "fftns" | "fast" | "fft" => AlignMode::Fftns2,
            "fftnsi" => AlignMode::Fftnsi,
            "linsi" | "localpair" | "accurate" => AlignMode::Linsi,
            "ginsi" | "globalpair" => AlignMode::Ginsi,
            "einsi" | "genafpair" => AlignMode::Einsi,
            _ => return None,
        })
    }

    /// MAFFT's name for the strategy.
    pub fn label(self) -> &'static str {
        match self {
            AlignMode::Auto => "auto",
            AlignMode::Fftns2 => "FFT-NS-2",
            AlignMode::Fftnsi => "FFT-NS-i",
            AlignMode::Linsi => "L-INS-i",
            AlignMode::Ginsi => "G-INS-i",
            AlignMode::Einsi => "E-INS-i",
        }
    }
}

/// How to align.
#[derive(Debug, Clone, Copy)]
pub struct AlignOptions {
    pub mode: AlignMode,
    /// Reverse-complement sequences that lie on the opposite strand to the
    /// first one before aligning (`mafft --adjustdirection`).
    pub adjust_direction: bool,
}

impl Default for AlignOptions {
    fn default() -> Self {
        AlignOptions { mode: AlignMode::Auto, adjust_direction: true }
    }
}

/// An alignment and what was done to get it.
#[derive(Debug, Clone)]
pub struct AlignOutcome {
    /// The aligned sequences, in input order, upper case, gaps as `-`.
    pub sequences: Vec<Sequence>,
    /// The strategy MAFFT ran (`FFT-NS-2`, `L-INS-i`, ...).
    pub strategy: &'static str,
    /// IDs of sequences that were reverse-complemented first.
    pub reversed: Vec<String>,
}

/// The strategy `mafft --auto` picks (MAFFT 7.526 `scripts/mafft`), as
/// (mode, guide-tree cycles, memory-saving tree).
fn decide_auto(nseq: usize, nlen: usize) -> (AlignmentMode, usize, bool) {
    if nlen < 3000 && nseq < 100 {
        (AlignmentMode::LInsi { iterations: 1000 }, 1, false)
    } else if nlen < 1000 && nseq < 200 {
        (AlignmentMode::LInsi { iterations: 2 }, 1, false)
    } else if nlen < 10000 && nseq < 500 {
        (AlignmentMode::FftNsi { iterations: 2 }, 2, false)
    } else if nseq < 20000 {
        (AlignmentMode::FftNs2, 2, false)
    } else if nseq < 100_000 {
        (AlignmentMode::FftNs2, 2, true)
    } else {
        (AlignmentMode::FftNs2, 1, true)
    }
}

/// The strategy [`AlignMode::Auto`] picks for `nseq` sequences whose
/// longest is `nlen` nt.
pub fn auto_strategy(nseq: usize, nlen: usize) -> &'static str {
    strategy_label(&decide_auto(nseq, nlen).0)
}

fn strategy_label(m: &AlignmentMode) -> &'static str {
    match m {
        AlignmentMode::FftNs2 => "FFT-NS-2",
        AlignmentMode::FftNsi { .. } => "FFT-NS-i",
        AlignmentMode::LInsi { .. } => "L-INS-i",
        AlignmentMode::GInsi { .. } => "G-INS-i",
        AlignmentMode::EInsi { .. } => "E-INS-i",
        _ => "MAFFT",
    }
}

/// The MAFFT engine settings for `mode` and the size of the input.
fn engine_for(mode: AlignMode, nseq: usize, nlen: usize) -> (MafftEngine, &'static str) {
    let (m, retree, memsavetree) = match mode {
        AlignMode::Auto => decide_auto(nseq, nlen),
        AlignMode::Fftns2 => (AlignmentMode::FftNs2, 2, false),
        AlignMode::Fftnsi => (AlignmentMode::FftNsi { iterations: 1000 }, 2, false),
        AlignMode::Linsi => (AlignmentMode::LInsi { iterations: 1000 }, 1, false),
        AlignMode::Ginsi => (AlignmentMode::GInsi { iterations: 1000 }, 1, false),
        AlignMode::Einsi => (AlignmentMode::EInsi { iterations: 1000 }, 1, false),
    };
    let label = strategy_label(&m);
    let mut engine = MafftEngine::new(m).with_retree(retree);
    engine.memsavetree = memsavetree;
    (engine, label)
}

/// Align `seqs` (gaps already in them are removed first).
pub fn align(seqs: &[Sequence], opts: &AlignOptions) -> Result<AlignOutcome, String> {
    if seqs.len() < 2 {
        return Err("Alignment needs at least two sequences.".into());
    }
    let mut set = SequenceSet::new(SeqType::Dna);
    for (i, s) in seqs.iter().enumerate() {
        // MAFFT's DNA convention: lower case, no gaps.
        let data: Vec<u8> =
            s.sequence.bytes().filter(|c| !matches!(c, b'-' | b'.')).map(|c| c.to_ascii_lowercase()).collect();
        if data.is_empty() {
            return Err(format!("Sequence '{}' is empty once its gaps are removed.", s.id));
        }
        set.sequences.push(MafftSeq { name: format!("s{i}"), data });
    }

    let mut reversed = Vec::new();
    if opts.adjust_direction {
        let adjusted = adjust_direction_mode(&set, AdjustMode::Kmer);
        if adjusted.sequences.len() == set.sequences.len() {
            for (i, s) in adjusted.sequences.iter().enumerate() {
                if s.name.starts_with("_R_") {
                    reversed.push(seqs[i].id.clone());
                }
            }
            set = adjusted;
        }
    }

    let nlen = set.max_len();
    let (engine, strategy) = engine_for(opts.mode, set.nseq(), nlen);
    let msa = run_engine(&engine, &set)?;
    if msa.sequences.len() != seqs.len() {
        return Err(format!("MAFFT returned {} sequences for {} input sequences.", msa.sequences.len(), seqs.len()));
    }
    let mut out = Vec::with_capacity(seqs.len());
    for (i, row) in msa.sequences.iter().enumerate() {
        let aligned: String =
            row.iter().map(|&c| if c == b'.' { '-' } else { (c as char).to_ascii_uppercase() }).collect();
        // The alignment must only add gaps to what went in.
        let before = set.sequences[i].data.iter().map(|c| (*c as char).to_ascii_uppercase());
        if !aligned.chars().filter(|&c| c != '-').eq(before) {
            return Err(format!("MAFFT changed the residues of sequence '{}'.", seqs[i].id));
        }
        out.push(Sequence::new(seqs[i].id.clone(), aligned));
    }
    Ok(AlignOutcome { sequences: out, strategy, reversed })
}

#[cfg(not(target_arch = "wasm32"))]
fn run_engine(engine: &MafftEngine, set: &SequenceSet) -> Result<crate::mafft::pipeline::MultipleAlignment, String> {
    // A panic in the aligner becomes an error message, not a crash of the
    // desktop app's server thread.
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| engine.align(set))).map_err(|e| {
        let why = e
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| e.downcast_ref::<&str>().copied())
            .unwrap_or("unknown error");
        format!("The aligner failed ({why}). Align the sequences with another tool and paste the alignment.")
    })
}

#[cfg(target_arch = "wasm32")]
fn run_engine(engine: &MafftEngine, set: &SequenceSet) -> Result<crate::mafft::pipeline::MultipleAlignment, String> {
    Ok(engine.align(set))
}

/// Check the aligner against C MAFFT: the bundled anfD alignment was made
/// with C MAFFT (L-INS-i), and aligning its ungapped sequences must give it
/// back exactly. Returns the strategy used.
#[cfg(not(target_arch = "wasm32"))]
pub fn self_test() -> Result<&'static str, String> {
    let (seqs, _) = crate::io::SequenceReader::parse(crate::EXAMPLE_ALIGNMENT)?;
    let out = align(&seqs, &AlignOptions::default())?;
    for (got, want) in out.sequences.iter().zip(&seqs) {
        if got.sequence != want.sequence {
            return Err(format!("the alignment of {} differs from C MAFFT's", want.id));
        }
    }
    Ok(out.strategy)
}

/// FASTA text of an alignment (60 columns per line).
pub fn to_fasta(seqs: &[Sequence]) -> String {
    let mut s = String::new();
    for q in seqs {
        s.push('>');
        s.push_str(&q.id);
        s.push('\n');
        let b = q.sequence.as_bytes();
        for chunk in b.chunks(60) {
            s.push_str(std::str::from_utf8(chunk).unwrap_or(""));
            s.push('\n');
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(id: &str, s: &str) -> Sequence {
        Sequence::new(id.into(), s.into())
    }

    #[test]
    fn matches_c_mafft_on_a_small_case() {
        let input = [seq("a", "ACGTACGTAAGG"), seq("b", "ACGTACTAAGG")];
        // C MAFFT 7.525: `mafft --auto` (L-INS-i here) and plain `mafft` (FFT-NS-2).
        let auto = align(&input, &AlignOptions { mode: AlignMode::Auto, adjust_direction: false }).unwrap();
        assert_eq!(auto.strategy, "L-INS-i");
        assert_eq!(auto.sequences[0].sequence, "ACGTACGTAAGG");
        assert_eq!(auto.sequences[1].sequence, "ACGTAC-TAAGG");
        let fast = align(&input, &AlignOptions { mode: AlignMode::Fftns2, adjust_direction: false }).unwrap();
        assert_eq!(fast.strategy, "FFT-NS-2");
        assert_eq!(fast.sequences[1].sequence, "ACGTACTAAGG-");
    }

    #[test]
    fn auto_follows_mafft_size_rules() {
        assert!(matches!(decide_auto(5, 1500).0, AlignmentMode::LInsi { iterations: 1000 }));
        assert!(matches!(decide_auto(150, 800).0, AlignmentMode::LInsi { iterations: 2 }));
        assert!(matches!(decide_auto(300, 5000).0, AlignmentMode::FftNsi { iterations: 2 }));
        assert!(matches!(decide_auto(50, 30000).0, AlignmentMode::FftNs2));
    }

    #[test]
    fn reverse_strand_sequences_are_flipped() {
        let fwd = "ATGGCTAGCTTGGACCAGTCAAGGCTTACGATCGGATCCAGTAGCTAGGCATCGATCGGATTACAGG";
        let rc: String = fwd
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                _ => 'C',
            })
            .collect();
        let out = align(
            &[seq("a", fwd), seq("b", fwd), seq("c", &rc)],
            &AlignOptions { mode: AlignMode::Fftns2, adjust_direction: true },
        )
        .unwrap();
        assert_eq!(out.reversed, vec!["c".to_string()]);
        assert_eq!(out.sequences[2].sequence, fwd);
    }

    #[test]
    fn reproduces_the_bundled_c_mafft_alignment() {
        assert_eq!(self_test().unwrap(), "L-INS-i");
    }

    #[test]
    fn mode_names_parse() {
        assert_eq!(AlignMode::parse("L-INS-i"), Some(AlignMode::Linsi));
        assert_eq!(AlignMode::parse("fft-ns-2"), Some(AlignMode::Fftns2));
        assert_eq!(AlignMode::parse("auto"), Some(AlignMode::Auto));
        assert_eq!(AlignMode::parse("bogus"), None);
    }
}
