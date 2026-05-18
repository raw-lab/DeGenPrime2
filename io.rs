use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use crate::models::{Sequence, DataSequence, DataNode, PrimerPair, Primer};

/// Sequence file format types
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    Fasta,
    Clustal,
}

/// Read sequences from a file
pub struct SequenceReader;

impl SequenceReader {
    /// Detect file format from content
    pub fn detect_format(path: &Path) -> Result<FileFormat, String> {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_lowercase();

        match ext.as_str() {
            "clust" | "aln" => Ok(FileFormat::Clustal),
            "fasta" | "fa" | "faa" | "fna" | "ffn" | "frn" => Ok(FileFormat::Fasta),
            _ => {
                // Try to detect from first line
                let file = File::open(path).map_err(|e| e.to_string())?;
                let reader = BufReader::new(file);
                let first_line = reader
                    .lines()
                    .next()
                    .ok_or("Empty file")??;

                if first_line.starts_with('>') {
                    Ok(FileFormat::Fasta)
                } else if first_line.contains("CLUSTAL") || first_line.starts_with(' ') {
                    Ok(FileFormat::Clustal)
                } else {
                    Err("Unknown format".to_string())
                }
            }
        }
    }

    /// Read FASTA format file
    pub fn read_fasta(path: &Path) -> Result<Vec<Sequence>, String> {
        let file = File::open(path).map_err(|e| e.to_string())?;
        let reader = BufReader::new(file);

        let mut sequences = Vec::new();
        let mut current_id = String::new();
        let mut current_seq = String::new();

        for line_result in reader.lines() {
            let line = line_result.map_err(|e| e.to_string())?;
            let trimmed = line.trim();

            if trimmed.is_empty() {
                continue;
            }

            if trimmed.starts_with('>') {
                // Save previous sequence if any
                if !current_id.is_empty() && !current_seq.is_empty() {
                    sequences.push(Sequence::new(current_id.clone(), current_seq.clone()));
                }

                // Start new sequence
                current_id = trimmed[1..].trim().to_string();
                current_seq.clear();
            } else {
                current_seq.push_str(trimmed);
            }
        }

        // Save last sequence
        if !current_id.is_empty() && !current_seq.is_empty() {
            sequences.push(Sequence::new(current_id, current_seq));
        }

        Ok(sequences)
    }

    /// Read CLUSTAL format file
    pub fn read_clustal(path: &Path) -> Result<Vec<Sequence>, String> {
        let file = File::open(path).map_err(|e| e.to_string())?;
        let reader = BufReader::new(file);

        let mut sequences = std::collections::HashMap::new();
        let mut skip_header = true;

        for line_result in reader.lines() {
            let line = line_result.map_err(|e| e.to_string())?;
            let trimmed = line.trim();

            if trimmed.is_empty() {
                continue;
            }

            if trimmed.contains("CLUSTAL") {
                skip_header = false;
                continue;
            }

            if skip_header || trimmed.starts_with(' ') {
                continue;
            }

            // CLUSTAL format: ID sequence [score]
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 2 {
                let id = parts[0].to_string();
                let seq = parts[1].to_string();
                sequences
                    .entry(id)
                    .or_insert_with(String::new)
                    .push_str(&seq);
            }
        }

        Ok(sequences
            .into_iter()
            .map(|(id, seq)| Sequence::new(id, seq))
            .collect())
    }

    /// Read sequences from file (auto-detect format)
    pub fn read(path: &Path) -> Result<Vec<Sequence>, String> {
        let format = Self::detect_format(path)?;

        match format {
            FileFormat::Fasta => Self::read_fasta(path),
            FileFormat::Clustal => Self::read_clustal(path),
        }
    }

    /// Check if sequences are aligned (all same length)
    pub fn is_aligned(sequences: &[Sequence]) -> bool {
        if sequences.is_empty() {
            return true;
        }

        let first_len = sequences[0].len();
        sequences.iter().all(|s| s.len() == first_len)
    }

    /// Convert sequences to DataSequence for processing
    pub fn to_data_sequences(sequences: &[Sequence]) -> Vec<DataSequence> {
        sequences
            .iter()
            .map(|seq| {
                let normalized = seq.normalize();
                let chars: Vec<char> = normalized.chars().collect();

                let nodes = (0..chars.len())
                    .scan(0, |col, _| {
                        // Build a column across all sequences
                        let col_idx = *col;
                        *col += 1;

                        // For a single sequence, create a DataNode from each position
                        if col_idx < normalized.len() {
                            Some(DataNode::new(vec![chars[col_idx]]))
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();

                // This is for single sequence - for alignment we'd need different logic
                if nodes.is_empty() {
                    DataSequence::new(chars.iter().map(|c| DataNode::new(vec![*c])).collect())
                } else {
                    DataSequence::new(nodes)
                }
            })
            .collect()
    }

    /// Convert aligned sequences to a single consensus DataSequence
    pub fn alignment_to_consensus(sequences: &[Sequence]) -> Result<DataSequence, String> {
        if sequences.is_empty() {
            return Err("No sequences provided".to_string());
        }

        if !Self::is_aligned(sequences) {
            return Err("Sequences are not aligned".to_string());
        }

        let seq_len = sequences[0].len();
        let mut nodes = Vec::new();

        // For each alignment column
        for col in 0..seq_len {
            let column: Vec<char> = sequences
                .iter()
                .filter_map(|seq| seq.normalize().chars().nth(col))
                .collect();

            if !column.is_empty() {
                nodes.push(DataNode::new(column));
            }
        }

        Ok(DataSequence::new(nodes))
    }
}

/// Output results to file
pub struct ResultWriter;

impl ResultWriter {
    /// Write primer pairs to CSV file
    pub fn write_csv(
        path: &Path,
        pairs: &[PrimerPair],
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
    ) -> Result<(), String> {
        let mut file = File::create(path).map_err(|e| e.to_string())?;

        // Write header
        writeln!(
            file,
            "Rank,Forward Primer,Forward Position,Forward Tm,Reverse Primer,Reverse Position,Reverse Tm,Amplicon Size (bp),Tm Difference"
        )
        .map_err(|e| e.to_string())?;

        // Write primers
        for (rank, pair) in pairs.iter().enumerate() {
            let fwd_primer = fwd_seq.sub_seq(pair.forward().index(), pair.forward().length());
            let rev_primer = rev_seq.sub_seq(pair.reverse().index(), pair.reverse().length());

            let fwd_tm = crate::thermodynamics::Thermodynamics::calculate_tm(
                &fwd_primer,
                &crate::settings::Settings::default(),
            );
            let rev_tm = crate::thermodynamics::Thermodynamics::calculate_tm(
                &rev_primer,
                &crate::settings::Settings::default(),
            );

            writeln!(
                file,
                "{},{},{},{:.2},{},{},{:.2},{},{}",
                rank + 1,
                fwd_primer.codes(),
                pair.forward().index(),
                fwd_tm,
                rev_primer.codes(),
                pair.reverse().index(),
                rev_tm,
                pair.amplicon_size(),
                pair.temp_diff()
            )
            .map_err(|e| e.to_string())?;
        }

        Ok(())
    }

    /// Write detailed analysis
    pub fn write_detailed_report(
        path: &Path,
        pairs: &[PrimerPair],
        fwd_seq: &DataSequence,
        rev_seq: &DataSequence,
    ) -> Result<(), String> {
        let mut file = File::create(path).map_err(|e| e.to_string())?;

        writeln!(file, "DeGenPrime Detailed Report").map_err(|e| e.to_string())?;
        writeln!(file, "=========================\n").map_err(|e| e.to_string())?;

        for (rank, pair) in pairs.iter().enumerate() {
            writeln!(file, "Primer Pair #{}", rank + 1).map_err(|e| e.to_string())?;
            writeln!(file, "-------------------").map_err(|e| e.to_string())?;

            let fwd_primer = fwd_seq.sub_seq(pair.forward().index(), pair.forward().length());
            let rev_primer = rev_seq.sub_seq(pair.reverse().index(), pair.reverse().length());

            writeln!(
                file,
                "Forward:  {} (position {})",
                fwd_primer.codes(),
                pair.forward().index()
            )
            .map_err(|e| e.to_string())?;
            writeln!(
                file,
                "Reverse:  {} (position {})",
                rev_primer.codes(),
                pair.reverse().index()
            )
            .map_err(|e| e.to_string())?;
            writeln!(
                file,
                "Amplicon: {} bp",
                pair.amplicon_size()
            )
            .map_err(|e| e.to_string())?;
            writeln!(
                file,
                "Tm Difference: {:.2}°C\n",
                pair.temp_diff()
            )
            .map_err(|e| e.to_string())?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_aligned() {
        let seqs = vec![
            Sequence::new("seq1".to_string(), "ATGC".to_string()),
            Sequence::new("seq2".to_string(), "CGTA".to_string()),
        ];
        assert!(SequenceReader::is_aligned(&seqs));

        let unaligned = vec![
            Sequence::new("seq1".to_string(), "ATGC".to_string()),
            Sequence::new("seq2".to_string(), "CG".to_string()),
        ];
        assert!(!SequenceReader::is_aligned(&unaligned));
    }
}
