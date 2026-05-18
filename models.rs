use std::cmp::Ordering;
use std::collections::HashMap;

/// Represents a single nucleotide with degeneracy information
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DataNode {
    /// IUPAC degenerate base code (A, T, G, C, W, S, R, Y, K, M, B, D, H, V, N)
    code: char,
    /// Most common nucleotide at this position
    most_common: char,
    /// Frequency ratio of most common to total [0.0, 1.0]
    ratio: f32,
}

impl DataNode {
    /// Create a DataNode from a vector of nucleotides at this alignment position
    pub fn new(chars: Vec<char>) -> Self {
        let filtered: Vec<char> = chars
            .into_iter()
            .filter(|&c| c != '-' && c != 'N' && c != 'n')
            .map(|c| c.to_ascii_uppercase())
            .collect();

        if filtered.is_empty() {
            return DataNode {
                code: 'N',
                most_common: 'N',
                ratio: 0.0,
            };
        }

        let mut counts = HashMap::new();
        for &c in &filtered {
            *counts.entry(c).or_insert(0) += 1;
        }

        let most_common = *counts
            .iter()
            .max_by_key(|&(_, count)| count)
            .map(|(c, _)| c)
            .unwrap_or(&'N');

        let ratio = counts[&most_common] as f32 / filtered.len() as f32;
        let code = Self::choose_code(&counts, most_common, ratio);

        DataNode {
            code,
            most_common,
            ratio,
        }
    }

    /// Choose the appropriate IUPAC code based on nucleotide frequencies
    fn choose_code(counts: &HashMap<char, usize>, most_common: char, ratio: f32) -> char {
        if counts.len() == 1 {
            most_common
        } else if counts.len() == 2 {
            let has_a = counts.contains_key(&'A');
            let has_t = counts.contains_key(&'T');
            let has_g = counts.contains_key(&'G');
            let has_c = counts.contains_key(&'C');

            match (has_a, has_t, has_g, has_c) {
                (true, true, false, false) => 'W',
                (false, false, true, true) => 'S',
                (true, false, true, false) => 'R',
                (false, true, false, true) => 'Y',
                (false, true, true, false) => 'K',
                (true, false, false, true) => 'M',
                (false, true, true, true) => 'B',
                (true, false, true, true) => 'D',
                (true, true, false, true) => 'H',
                (true, true, true, false) => 'V',
                _ => 'N',
            }
        } else {
            'N'
        }
    }

    /// Get the nucleotide code
    pub fn code(&self) -> char {
        self.code
    }

    /// Get the most common nucleotide
    pub fn most_common(&self) -> char {
        self.most_common
    }

    /// Get the frequency ratio
    pub fn ratio(&self) -> f32 {
        self.ratio
    }

    /// Inverse complement: A↔T, C↔G
    pub fn inverse(&self) -> Self {
        let inv_code = Self::inverse_complement(self.code);
        let inv_common = Self::inverse_complement(self.most_common);
        DataNode {
            code: inv_code,
            most_common: inv_common,
            ratio: self.ratio,
        }
    }

    /// Get the inverse complement of a nucleotide
    pub fn inverse_complement(base: char) -> char {
        match base.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'W' => 'W',
            'S' => 'S',
            'R' => 'Y',
            'Y' => 'R',
            'K' => 'M',
            'M' => 'K',
            'B' => 'V',
            'V' => 'B',
            'D' => 'H',
            'H' => 'D',
            '-' => '-',
            _ => 'N',
        }
    }

    /// Get enthalpy value for nearest-neighbor thermodynamics
    pub fn enthalpy_to(&self, other: &DataNode) -> f32 {
        Self::nn_enthalpy(self.most_common, other.most_common)
    }

    /// Get entropy value for nearest-neighbor thermodynamics
    pub fn entropy_to(&self, other: &DataNode) -> f32 {
        Self::nn_entropy(self.most_common, other.most_common)
    }

    /// Get Gibbs free energy for nearest-neighbor thermodynamics
    pub fn gibbs_to(&self, other: &DataNode) -> f32 {
        Self::nn_gibbs(self.most_common, other.most_common)
    }

    // Nearest-neighbor thermodynamic values (from PNAS 1998)
    fn nn_enthalpy(base1: char, base2: char) -> f32 {
        let pair = format!("{}{}", base1, base2);
        match pair.as_str() {
            "AA" | "TT" => -7.9,
            "AT" => -7.2,
            "TA" => -7.2,
            "CA" | "TG" => -8.5,
            "GT" | "AC" => -8.4,
            "CT" | "AG" => -7.8,
            "GA" | "TC" => -8.2,
            "CG" => -10.6,
            "GC" => -9.8,
            "GG" | "CC" => -8.0,
            _ => 0.0,
        }
    }

    fn nn_entropy(base1: char, base2: char) -> f32 {
        let pair = format!("{}{}", base1, base2);
        match pair.as_str() {
            "AA" | "TT" => -24.4,
            "AT" => -20.4,
            "TA" => -21.3,
            "CA" | "TG" => -22.7,
            "GT" | "AC" => -22.4,
            "CT" | "AG" => -21.0,
            "GA" | "TC" => -22.2,
            "CG" => -27.8,
            "GC" => -24.4,
            "GG" | "CC" => -19.9,
            _ => 0.0,
        }
    }

    fn nn_gibbs(base1: char, base2: char) -> f32 {
        let pair = format!("{}{}", base1, base2);
        match pair.as_str() {
            "AA" | "TT" => -1.03,
            "AT" => -0.88,
            "TA" => -0.58,
            "CA" | "TG" => -1.45,
            "GT" | "AC" => -1.44,
            "CT" | "AG" => -1.28,
            "GA" | "TC" => -1.30,
            "CG" => -2.17,
            "GC" => -2.24,
            "GG" | "CC" => -1.84,
            _ => 0.0,
        }
    }
}

/// A data sequence - a vector of DataNodes representing aligned positions
#[derive(Debug, Clone)]
pub struct DataSequence {
    nodes: Vec<DataNode>,
    /// Starting position in the original sequence
    start_pos: usize,
}

impl DataSequence {
    /// Create a DataSequence from a vector of DataNodes
    pub fn new(nodes: Vec<DataNode>) -> Self {
        DataSequence {
            nodes,
            start_pos: 0,
        }
    }

    /// Create a DataSequence with a starting position
    pub fn with_start(nodes: Vec<DataNode>, start_pos: usize) -> Self {
        DataSequence { nodes, start_pos }
    }

    /// Get the length of the sequence
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Check if sequence is empty
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Get a reference to the nodes
    pub fn nodes(&self) -> &[DataNode] {
        &self.nodes
    }

    /// Get a subsequence
    pub fn sub_seq(&self, start: usize, length: usize) -> DataSequence {
        let end = (start + length).min(self.nodes.len());
        let sub_nodes = self.nodes[start..end].to_vec();
        DataSequence::with_start(sub_nodes, self.start_pos + start)
    }

    /// Get the actual size (without gap characters)
    pub fn actual_size(&self) -> usize {
        self.nodes.iter().filter(|n| n.code() != '-').count()
    }

    /// Get the reverse sequence
    pub fn reverse(&self) -> DataSequence {
        let mut rev = self.nodes.clone();
        rev.reverse();
        DataSequence::new(rev)
    }

    /// Get the inverse complement
    pub fn inverse_complement(&self) -> DataSequence {
        let inv = self
            .nodes
            .iter()
            .map(|n| n.inverse())
            .collect::<Vec<_>>();
        DataSequence::new(inv)
    }

    /// Get the reverse complement (inverse + reverse)
    pub fn reverse_complement(&self) -> DataSequence {
        self.inverse_complement().reverse()
    }

    /// Get sequence as string of codes
    pub fn codes(&self) -> String {
        self.nodes.iter().map(|n| n.code()).collect()
    }

    /// Get GC content ratio
    pub fn gc_ratio(&self) -> f32 {
        let gc_count = self
            .nodes
            .iter()
            .filter(|n| matches!(n.most_common(), 'G' | 'C'))
            .count();
        gc_count as f32 / self.actual_size().max(1) as f32
    }

    /// Calculate penalty based on degeneracy
    pub fn penalty(&self) -> f32 {
        self.nodes
            .iter()
            .map(|n| {
                // Penalty increases with lower frequency of most common base
                // and with degeneracy (more ambiguous codes)
                (1.0 - n.ratio()) * 100.0
            })
            .sum()
    }

    /// Count degenerate positions
    pub fn degeneracy_count(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| !matches!(n.code(), 'A' | 'T' | 'G' | 'C' | '-'))
            .count()
    }

    /// Count gaps (deletions)
    pub fn gap_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.code() == '-').count()
    }

    /// Check for complementary ends (potential secondary structures)
    pub fn has_complementary_ends(&self) -> bool {
        let (first3, last3) = (self.len().min(3), self.len().min(3));
        let start = self.nodes.iter().take(first3);
        let end = self
            .nodes
            .iter()
            .rev()
            .take(last3)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>();

        start.zip(end.iter()).all(|(a, b)| {
            DataNode::inverse_complement(a.code()) == b.code()
        })
    }
}

/// A single primer with position and length
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Primer {
    index: usize,
    length: usize,
    penalty: f32,
}

impl Primer {
    /// Create a new primer at a position with a length
    pub fn new(index: usize, length: usize) -> Self {
        Primer {
            index,
            length,
            penalty: 0.0,
        }
    }

    /// Set the penalty value
    pub fn with_penalty(mut self, penalty: f32) -> Self {
        self.penalty = penalty;
        self
    }

    /// Get the starting index
    pub fn index(&self) -> usize {
        self.index
    }

    /// Get the length
    pub fn length(&self) -> usize {
        self.length
    }

    /// Get the penalty
    pub fn penalty(&self) -> f32 {
        self.penalty
    }

    /// Get the ending index (inclusive)
    pub fn end_index(&self) -> usize {
        self.index + self.length - 1
    }
}

impl Ord for Primer {
    fn cmp(&self, other: &Self) -> Ordering {
        self.penalty
            .partial_cmp(&other.penalty)
            .unwrap_or(Ordering::Equal)
    }
}

impl PartialOrd for Primer {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for Primer {}

/// A pair of primers (forward and reverse)
#[derive(Debug, Clone, Copy)]
pub struct PrimerPair {
    forward: Primer,
    reverse: Primer,
    /// Amplicon length in base pairs
    amplicon_size: usize,
    /// Melting temperature difference
    temp_diff: f32,
}

impl PrimerPair {
    /// Create a new primer pair
    pub fn new(forward: Primer, reverse: Primer, amplicon_size: usize, temp_diff: f32) -> Self {
        PrimerPair {
            forward,
            reverse,
            amplicon_size,
            temp_diff,
        }
    }

    /// Get the forward primer
    pub fn forward(&self) -> &Primer {
        &self.forward
    }

    /// Get the reverse primer
    pub fn reverse(&self) -> &Primer {
        &self.reverse
    }

    /// Get amplicon size
    pub fn amplicon_size(&self) -> usize {
        self.amplicon_size
    }

    /// Get temperature difference
    pub fn temp_diff(&self) -> f32 {
        self.temp_diff
    }
}

impl Ord for PrimerPair {
    fn cmp(&self, other: &Self) -> Ordering {
        self.temp_diff
            .partial_cmp(&other.temp_diff)
            .unwrap_or(Ordering::Equal)
    }
}

impl PartialOrd for PrimerPair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for PrimerPair {}

impl PartialEq for PrimerPair {
    fn eq(&self, other: &Self) -> bool {
        self.temp_diff == other.temp_diff
    }
}

/// A biological sequence with optional identifier
#[derive(Debug, Clone)]
pub struct Sequence {
    id: String,
    sequence: String,
}

impl Sequence {
    /// Create a new sequence
    pub fn new(id: String, sequence: String) -> Self {
        Sequence { id, sequence }
    }

    /// Get the identifier
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Get the sequence string
    pub fn sequence(&self) -> &str {
        &self.sequence
    }

    /// Get length
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Normalize sequence (uppercase, remove whitespace)
    pub fn normalize(&self) -> String {
        self.sequence
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(|c| c.to_ascii_uppercase())
            .collect()
    }

    /// Get the reverse complement
    pub fn reverse_complement(&self) -> String {
        self.sequence
            .chars()
            .rev()
            .map(|c| DataNode::inverse_complement(c.to_ascii_uppercase()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_datanode_inverse() {
        let node = DataNode::new(vec!['A', 'A', 'A']);
        let inv = node.inverse();
        assert_eq!(inv.code(), 'T');
    }

    #[test]
    fn test_datasequence_gc_ratio() {
        let nodes = vec![
            DataNode::new(vec!['G']),
            DataNode::new(vec!['C']),
            DataNode::new(vec!['A']),
            DataNode::new(vec!['T']),
        ];
        let seq = DataSequence::new(nodes);
        assert!((seq.gc_ratio() - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_primer_pair_ordering() {
        let p1 = PrimerPair::new(
            Primer::new(0, 20),
            Primer::new(100, 20),
            100,
            0.5,
        );
        let p2 = PrimerPair::new(
            Primer::new(0, 20),
            Primer::new(100, 20),
            100,
            1.5,
        );
        assert!(p1 < p2);
    }
}
