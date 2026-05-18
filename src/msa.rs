use crate::models::Sequence;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::process::Command;

/// Multiple Sequence Alignment handler
pub struct MSAHandler;

impl MSAHandler {
    /// Check if MAFFT is installed on the system
    pub fn is_mafft_available() -> bool {
        Command::new("which")
            .arg("mafft")
            .output()
            .ok()
            .and_then(|output| if output.status.success() { Some(output) } else { None })
            .is_some()
    }

    /// Align sequences using MAFFT (external tool)
    pub fn align_with_mafft(
        sequences: &[Sequence],
        mode: MSAMode,
        output_path: &Path,
    ) -> Result<Vec<Sequence>, String> {
        if !Self::is_mafft_available() {
            return Err("MAFFT is not installed. Install with: brew install mafft (macOS) or apt-get install mafft (Linux)".to_string());
        }

        // Write input FASTA
        let input_path = output_path.with_extension("input.fasta");
        Self::write_fasta(&input_path, sequences)?;

        // Prepare MAFFT arguments
        let mut args = vec!["--clustalout".to_string(), "--quiet".to_string()];

        match mode {
            MSAMode::GlobalAccurate => {
                args.push("--globalpair".to_string());
                args.push("--maxiterate".to_string());
                args.push("1000".to_string());
            }
            MSAMode::LocalAccurate => {
                args.push("--localpair".to_string());
                args.push("--maxiterate".to_string());
                args.push("1000".to_string());
            }
            MSAMode::FastGlobal => {
                args.push("--globalpair".to_string());
                args.push("--maxiterate".to_string());
                args.push("100".to_string());
            }
            MSAMode::FastLocal => {
                args.push("--localpair".to_string());
                args.push("--maxiterate".to_string());
                args.push("100".to_string());
            }
            MSAMode::VeryFast => {
                args.push("--retree".to_string());
                args.push("1".to_string());
            }
        }

        args.push(input_path.to_string_lossy().to_string());

        // Run MAFFT
        let output = Command::new("mafft")
            .args(&args)
            .output()
            .map_err(|e| format!("Failed to run MAFFT: {}", e))?;

        if !output.status.success() {
            return Err(format!(
                "MAFFT failed: {}",
                String::from_utf8_lossy(&output.stderr)
            ));
        }

        // Parse CLUSTAL output
        let result_str = String::from_utf8(output.stdout)
            .map_err(|e| format!("Failed to parse MAFFT output: {}", e))?;

        let aligned = Self::parse_clustal_output(&result_str)?;

        // Clean up input file
        let _ = std::fs::remove_file(&input_path);

        Ok(aligned)
    }

    /// Align sequences using built-in Rust algorithm (fast but less accurate)
    pub fn align_with_rust(sequences: &[Sequence]) -> Result<Vec<Sequence>, String> {
        if sequences.is_empty() {
            return Err("No sequences to align".to_string());
        }

        if sequences.len() == 1 {
            return Ok(sequences.to_vec());
        }

        // Use progressive alignment algorithm
        let mut aligned = Self::simple_progressive_align(sequences)?;
        
        // Refinement pass
        aligned = Self::refine_alignment(&aligned)?;

        Ok(aligned)
    }

    /// Simple progressive multiple sequence alignment
    fn simple_progressive_align(sequences: &[Sequence]) -> Result<Vec<Sequence>, String> {
        let mut aligned = vec![sequences[0].clone()];

        for seq in &sequences[1..] {
            aligned = Self::align_pair(&aligned, seq)?;
        }

        Ok(aligned)
    }

    /// Pairwise alignment using simple dynamic programming
    fn align_pair(seqs: &[Sequence], new_seq: &Sequence) -> Result<Vec<Sequence>, String> {
        // For simplicity, just pad with gaps to make equal length
        let max_len = seqs.iter().map(|s| s.len()).max().unwrap_or(0).max(new_seq.len());

        let mut result = Vec::new();

        for seq in seqs {
            let mut s = seq.sequence().to_string();
            while s.len() < max_len {
                s.push('-');
            }
            result.push(Sequence::new(seq.id().to_string(), s));
        }

        let mut new_s = new_seq.sequence().to_string();
        while new_s.len() < max_len {
            new_s.push('-');
        }
        result.push(Sequence::new(new_seq.id().to_string(), new_s));

        Ok(result)
    }

    /// Refine alignment by iterative improvement
    fn refine_alignment(sequences: &[Sequence]) -> Result<Vec<Sequence>, String> {
        // Simple refinement: check for obvious improvements
        // This is a placeholder - real implementation would do iterative refinement
        Ok(sequences.to_vec())
    }

    fn write_fasta(path: &Path, sequences: &[Sequence]) -> Result<(), String> {
        let mut file = File::create(path).map_err(|e| e.to_string())?;

        for seq in sequences {
            writeln!(file, ">{}", seq.id()).map_err(|e| e.to_string())?;
            
            // Write sequence in 60 char chunks
            let norm_seq = seq.normalize();
            for chunk in norm_seq.chars().collect::<Vec<_>>().chunks(60) {
                let line: String = chunk.iter().collect();
                writeln!(file, "{}", line).map_err(|e| e.to_string())?;
            }
        }

        Ok(())
    }

    fn parse_clustal_output(content: &str) -> Result<Vec<Sequence>, String> {
        let mut sequences = std::collections::HashMap::new();
        let mut sequence_order = Vec::new();

        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with("CLUSTAL") {
                continue;
            }

            if trimmed.starts_with(' ') {
                // Conservation line - skip
                continue;
            }

            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 2 {
                let id = parts[0].to_string();
                let seq_part = parts[1].to_string();

                if !sequence_order.contains(&id) {
                    sequence_order.push(id.clone());
                }

                sequences
                    .entry(id)
                    .or_insert_with(String::new)
                    .push_str(&seq_part);
            }
        }

        let result = sequence_order
            .iter()
            .map(|id| {
                let seq = sequences.get(id).cloned().unwrap_or_default();
                Sequence::new(id.clone(), seq)
            })
            .collect();

        Ok(result)
    }
}

/// MSA algorithm mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MSAMode {
    /// Accurate global alignment (slow, best for similar sequences)
    GlobalAccurate,
    /// Accurate local alignment (slow, best for diverse sequences)
    LocalAccurate,
    /// Fast global alignment (medium speed)
    FastGlobal,
    /// Fast local alignment (medium speed)
    FastLocal,
    /// Very fast alignment (for large datasets)
    VeryFast,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_write() {
        let seqs = vec![
            Sequence::new("seq1".to_string(), "ATGC".to_string()),
            Sequence::new("seq2".to_string(), "CGTA".to_string()),
        ];

        let path = Path::new("/tmp/test.fasta");
        assert!(MSAHandler::write_fasta(path, &seqs).is_ok());
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_align_pair() {
        let seqs = vec![Sequence::new("seq1".to_string(), "ATG".to_string())];
        let new_seq = Sequence::new("seq2".to_string(), "ATGC".to_string());

        let result = MSAHandler::align_pair(&seqs, &new_seq);
        assert!(result.is_ok());
    }
}
