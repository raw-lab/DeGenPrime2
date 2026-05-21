//! Core data models for DeGenPrime.
//!
//! `DataNode` represents a single aligned column with IUPAC degeneracy info.
//! `DataSequence` is a vector of DataNodes (either a single sequence or a consensus
//! built from an MSA).
//!
//! Thermodynamic nearest-neighbour parameters follow SantaLucia (1998) PNAS 95(4):1460.

use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::HashMap;

/// A single aligned column. Tracks the IUPAC code, the most common base, and
/// the frequency ratio of that most-common base (1.0 means fully conserved).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DataNode {
    code: char,
    most_common: char,
    ratio: f32,
}

impl DataNode {
    /// Build a DataNode from the bases observed at this alignment column.
    /// Gaps (`-`) and Ns are filtered out before counting.
    pub fn new(chars: Vec<char>) -> Self {
        let filtered: Vec<char> = chars
            .into_iter()
            .filter(|c| *c != '-' && *c != 'N' && *c != 'n')
            .map(|c| c.to_ascii_uppercase())
            .collect();

        if filtered.is_empty() {
            return DataNode {
                code: 'N',
                most_common: 'N',
                ratio: 0.0,
            };
        }

        let mut counts: HashMap<char, usize> = HashMap::new();
        for &c in &filtered {
            *counts.entry(c).or_insert(0) += 1;
        }

        let most_common = *counts
            .iter()
            .max_by_key(|&(_, v)| v)
            .map(|(c, _)| c)
            .unwrap_or(&'N');

        let ratio = counts[&most_common] as f32 / filtered.len() as f32;
        let code = Self::choose_iupac(&counts);

        DataNode { code, most_common, ratio }
    }

    /// Convenience: a fully-conserved column with the given base.
    pub fn from_base(base: char) -> Self {
        DataNode::new(vec![base.to_ascii_uppercase()])
    }

    fn choose_iupac(counts: &HashMap<char, usize>) -> char {
        // Bases actually observed in this column.
        let has = |b: char| counts.contains_key(&b);
        let (a, t, g, c) = (has('A'), has('T'), has('G'), has('C'));
        match (a, t, g, c) {
            (true, false, false, false) => 'A',
            (false, true, false, false) => 'T',
            (false, false, true, false) => 'G',
            (false, false, false, true) => 'C',
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
    }

    pub fn code(&self) -> char { self.code }
    pub fn most_common(&self) -> char { self.most_common }
    pub fn ratio(&self) -> f32 { self.ratio }

    /// Returns the IUPAC complement (single-letter mapping).
    /// Note: for ambiguous bases this returns the complementary ambiguity code.
    pub fn complement_base(base: char) -> char {
        match base.to_ascii_uppercase() {
            'A' => 'T', 'T' => 'A', 'U' => 'A',
            'G' => 'C', 'C' => 'G',
            'W' => 'W', 'S' => 'S',
            'R' => 'Y', 'Y' => 'R',
            'K' => 'M', 'M' => 'K',
            'B' => 'V', 'V' => 'B',
            'D' => 'H', 'H' => 'D',
            'N' => 'N', '-' => '-',
            other => other,
        }
    }

    /// Return a node with the complementary base (same ratio).
    pub fn complement(&self) -> Self {
        DataNode {
            code: Self::complement_base(self.code),
            most_common: Self::complement_base(self.most_common),
            ratio: self.ratio,
        }
    }

    // ---- nearest-neighbour thermodynamics (SantaLucia 1998) -------------
    //
    // Values are per-stack contributions in kcal/mol (ΔH) and cal/(mol·K) (ΔS).
    // When the column's most_common base is ambiguous (N or a missing IUPAC
    // code) we return 0 — these stacks contribute nothing to the Tm, which is
    // a deliberate conservative behaviour.

    /// Enthalpy contribution for the stack (self -> other).
    pub fn enthalpy_to(&self, other: &DataNode) -> f32 {
        Self::nn_enthalpy(self.most_common, other.most_common)
    }

    /// Entropy contribution for the stack (self -> other).
    pub fn entropy_to(&self, other: &DataNode) -> f32 {
        Self::nn_entropy(self.most_common, other.most_common)
    }

    fn nn_enthalpy(b1: char, b2: char) -> f32 {
        match (b1, b2) {
            ('A', 'A') | ('T', 'T') => -7.9,
            ('A', 'T') => -7.2,
            ('T', 'A') => -7.2,
            ('C', 'A') | ('T', 'G') => -8.5,
            ('G', 'T') | ('A', 'C') => -8.4,
            ('C', 'T') | ('A', 'G') => -7.8,
            ('G', 'A') | ('T', 'C') => -8.2,
            ('C', 'G') => -10.6,
            ('G', 'C') => -9.8,
            ('G', 'G') | ('C', 'C') => -8.0,
            _ => 0.0,
        }
    }

    fn nn_entropy(b1: char, b2: char) -> f32 {
        match (b1, b2) {
            ('A', 'A') | ('T', 'T') => -22.2,
            ('A', 'T') => -20.4,
            ('T', 'A') => -21.3,
            ('C', 'A') | ('T', 'G') => -22.7,
            ('G', 'T') | ('A', 'C') => -22.4,
            ('C', 'T') | ('A', 'G') => -21.0,
            ('G', 'A') | ('T', 'C') => -22.2,
            ('C', 'G') => -27.2,
            ('G', 'C') => -24.4,
            ('G', 'G') | ('C', 'C') => -19.9,
            _ => 0.0,
        }
    }
}

/// A sequence of DataNodes representing either a real sequence or a consensus.
#[derive(Debug, Clone)]
pub struct DataSequence {
    nodes: Vec<DataNode>,
    start_pos: usize,
}

impl DataSequence {
    pub fn new(nodes: Vec<DataNode>) -> Self {
        DataSequence { nodes, start_pos: 0 }
    }

    pub fn with_start(nodes: Vec<DataNode>, start_pos: usize) -> Self {
        DataSequence { nodes, start_pos }
    }

    /// Build directly from a plain DNA string (e.g. for primer testing).
    pub fn from_string(seq: &str) -> Self {
        let nodes = seq
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(|c| DataNode::from_base(c))
            .collect();
        DataSequence::new(nodes)
    }

    pub fn len(&self) -> usize { self.nodes.len() }
    pub fn is_empty(&self) -> bool { self.nodes.is_empty() }
    pub fn nodes(&self) -> &[DataNode] { &self.nodes }
    pub fn start_pos(&self) -> usize { self.start_pos }

    /// Slice a subsequence by start position and length.
    pub fn sub_seq(&self, start: usize, length: usize) -> DataSequence {
        let end = (start + length).min(self.nodes.len());
        let sub_nodes = self.nodes[start..end].to_vec();
        DataSequence::with_start(sub_nodes, self.start_pos + start)
    }

    /// Number of non-gap positions.
    pub fn actual_size(&self) -> usize {
        self.nodes.iter().filter(|n| n.code() != '-').count()
    }

    /// Returns a new sequence with positions reversed.
    pub fn reverse(&self) -> DataSequence {
        let mut rev = self.nodes.clone();
        rev.reverse();
        DataSequence::new(rev)
    }

    /// Base-by-base IUPAC complement (no reversal).
    pub fn complement(&self) -> DataSequence {
        let comp = self.nodes.iter().map(|n| n.complement()).collect();
        DataSequence::new(comp)
    }

    /// Reverse complement: the standard DNA operation.
    pub fn reverse_complement(&self) -> DataSequence {
        self.complement().reverse()
    }

    /// Sequence as a string of IUPAC codes.
    pub fn codes(&self) -> String {
        self.nodes.iter().map(|n| n.code()).collect()
    }

    /// Sequence using only the most-common base at each column (no degeneracy).
    pub fn consensus_string(&self) -> String {
        self.nodes.iter().map(|n| n.most_common()).collect()
    }

    /// GC ratio computed against non-gap positions; degenerate codes that
    /// always pair with G or C (`S`, `B`, `V`) count toward GC. Returns 0
    /// for an all-gap sequence.
    pub fn gc_ratio(&self) -> f32 {
        let total = self.actual_size().max(1);
        let gc_count = self
            .nodes
            .iter()
            .filter(|n| matches!(n.most_common(), 'G' | 'C'))
            .count();
        gc_count as f32 / total as f32
    }

    /// Penalty score for ranking: low conservation (low ratio) hurts the score.
    /// Lower penalty = better.
    pub fn penalty(&self) -> f32 {
        self.nodes
            .iter()
            .map(|n| (1.0 - n.ratio()) * 100.0)
            .sum::<f32>()
    }

    /// Number of degenerate (ambiguous) positions in the consensus.
    pub fn degeneracy_count(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| !matches!(n.code(), 'A' | 'T' | 'G' | 'C' | '-'))
            .count()
    }

    /// Number of gap positions.
    pub fn gap_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.code() == '-').count()
    }

    /// Check the last 3 bases against the first 3 (reversed) — a quick
    /// secondary-structure heuristic for 3' end hairpin potential.
    pub fn has_complementary_ends(&self) -> bool {
        if self.len() < 6 {
            return false;
        }
        let n = self.len();
        let head: Vec<char> = self.nodes[..3].iter().map(|x| x.most_common()).collect();
        let tail: Vec<char> = self.nodes[n - 3..].iter().map(|x| x.most_common()).collect();
        // Are head and reversed tail base-pair complements at all 3 positions?
        head.iter()
            .zip(tail.iter().rev())
            .all(|(a, b)| DataNode::complement_base(*a) == *b)
    }
}

/// A candidate primer — a slice of a parent sequence.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Primer {
    pub index: usize,
    pub length: usize,
    pub penalty: f32,
}

impl Primer {
    pub fn new(index: usize, length: usize) -> Self {
        Primer { index, length, penalty: 0.0 }
    }

    pub fn with_penalty(mut self, penalty: f32) -> Self {
        self.penalty = penalty;
        self
    }

    pub fn end_index(&self) -> usize {
        self.index + self.length - 1
    }
}

impl Eq for Primer {}

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

/// A forward + reverse primer pair.
///
/// IMPORTANT: `reverse.index` is the offset on the **reverse-complement**
/// sequence. To translate that to a position on the original (forward) strand
/// use `PrimerPair::reverse_fwd_strand_position(template_len)`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PrimerPair {
    pub forward: Primer,
    pub reverse: Primer,
    pub amplicon_size: usize,
    pub tm_diff: f32,
    pub fwd_tm: f32,
    pub rev_tm: f32,
}

impl PrimerPair {
    pub fn new(
        forward: Primer,
        reverse: Primer,
        amplicon_size: usize,
        tm_diff: f32,
        fwd_tm: f32,
        rev_tm: f32,
    ) -> Self {
        PrimerPair { forward, reverse, amplicon_size, tm_diff, fwd_tm, rev_tm }
    }

    /// Convert the reverse primer's index (which is in reverse-complement
    /// coordinates) to its 5' position on the original forward strand.
    pub fn reverse_fwd_strand_position(&self, template_len: usize) -> usize {
        // The reverse primer covers positions [rev.index, rev.index+rev.length)
        // on the RC strand. On the forward strand this corresponds to
        // [template_len - rev.index - rev.length, template_len - rev.index).
        // We report the 5'-most position on the forward strand of the
        // reverse-primer footprint (i.e. where on the template the primer
        // anneals in 3'->5' direction).
        template_len
            .saturating_sub(self.reverse.index)
            .saturating_sub(self.reverse.length)
    }
}

impl PartialEq for PrimerPair {
    fn eq(&self, other: &Self) -> bool {
        self.forward == other.forward && self.reverse == other.reverse
    }
}
impl Eq for PrimerPair {}

impl Ord for PrimerPair {
    fn cmp(&self, other: &Self) -> Ordering {
        // Lower tm_diff is better; tie-break on combined penalty.
        self.tm_diff
            .partial_cmp(&other.tm_diff)
            .unwrap_or(Ordering::Equal)
            .then_with(|| {
                (self.forward.penalty + self.reverse.penalty)
                    .partial_cmp(&(other.forward.penalty + other.reverse.penalty))
                    .unwrap_or(Ordering::Equal)
            })
    }
}

impl PartialOrd for PrimerPair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A raw biological sequence read from a file (FASTA / CLUSTAL).
#[derive(Debug, Clone)]
pub struct Sequence {
    pub id: String,
    pub sequence: String,
}

impl Sequence {
    pub fn new(id: String, sequence: String) -> Self {
        Sequence { id, sequence }
    }

    pub fn len(&self) -> usize { self.sequence.len() }
    pub fn is_empty(&self) -> bool { self.sequence.is_empty() }

    /// Uppercase, strip whitespace.
    pub fn normalize(&self) -> String {
        self.sequence
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(|c| c.to_ascii_uppercase())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_base_works_for_ambig() {
        assert_eq!(DataNode::complement_base('A'), 'T');
        assert_eq!(DataNode::complement_base('R'), 'Y');
        assert_eq!(DataNode::complement_base('N'), 'N');
        assert_eq!(DataNode::complement_base('-'), '-');
    }

    #[test]
    fn datanode_picks_iupac_pair() {
        let node = DataNode::new(vec!['A', 'A', 'A', 'G']);
        assert_eq!(node.code(), 'R'); // A/G is R
        assert_eq!(node.most_common(), 'A');
        assert!((node.ratio() - 0.75).abs() < 1e-6);
    }

    #[test]
    fn reverse_complement_is_involutive() {
        let s = DataSequence::from_string("ATGCAATTGC");
        let rc = s.reverse_complement();
        let back = rc.reverse_complement();
        assert_eq!(s.codes(), back.codes());
    }

    #[test]
    fn gc_ratio_correct() {
        let s = DataSequence::from_string("GGCC");
        assert!((s.gc_ratio() - 1.0).abs() < 1e-6);
        let s2 = DataSequence::from_string("ATAT");
        assert!((s2.gc_ratio() - 0.0).abs() < 1e-6);
        let s3 = DataSequence::from_string("ATGC");
        assert!((s3.gc_ratio() - 0.5).abs() < 1e-6);
    }
}
