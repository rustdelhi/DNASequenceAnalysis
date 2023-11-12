//! This module is used to align two or more DNA/RNA sequences
//! to "align" them, see: https://en.wikipedia.org/wiki/Sequence_alignment

use std::fmt::Display;

use bio::alignment::{
    distance::{hamming, levenshtein},
    pairwise::{MatchFunc, Scoring},
    Alignment,
};

type PairwiseAlignment = bio::alignment::Alignment;
type PartialorderAlignment = bio::alignment::poa::Alignment;

///  Scoring rule for [Substitution matrix](https://en.wikipedia.org/wiki/Smith_Waterman_algorithm#Substitution_matrix)
#[derive(Debug, Clone)]
pub struct Score {
    r#match: i32,
    miss_match: i32,
}

impl Display for Score {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Score(match={},miss-match={})",
            self.r#match, self.miss_match
        )
    }
}

impl Score {
    pub fn new(r#match: i32, miss_match: i32) -> Self {
        tracing::info!(
            "Generating Score match={} miss-match={}",
            r#match,
            miss_match
        );
        Self {
            r#match,
            miss_match,
        }
    }
}

impl From<(i32, i32)> for Score {
    fn from(value: (i32, i32)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl MatchFunc for Score {
    fn score(&self, a: u8, b: u8) -> i32 {
        if a == b {
            self.r#match
        } else {
            self.miss_match
        }
    }
}

/// Specifying gap penalty for Smith Waterman algorithm
/// See: https://en.wikipedia.org/wiki/Smith_Waterman_algorithm#Gap_penalty
#[derive(Debug)]
pub struct GapPanelty {
    pub open: i32,
    pub extend: i32,
}

impl Display for GapPanelty {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Gap(open={},extend={})", self.open, self.extend)
    }
}

impl GapPanelty {
    pub fn new(open: i32, extend: i32) -> Self {
        assert!(open < 0, "Gap open penalty cant be positive");
        assert!(extend < 0, "Gap extend penalty cant be positive");
        tracing::info!("Generating GapPanelty open={} extend={}", open, extend);
        Self { open, extend }
    }
}

impl From<(i32, i32)> for GapPanelty {
    fn from(value: (i32, i32)) -> Self {
        Self::new(value.0, value.1)
    }
}

/// Compare two sequences and align them
#[derive(Debug)]
pub struct DiffStat<'seq, F>
where
    F: MatchFunc + Clone + Display,
{
    /// Master sequence
    reference: &'seq [u8],
    /// The one that will be aligned
    query: &'seq [u8],
    /// Gap penalty used while alignment
    gap_penalty: GapPanelty,
    /// Match and miss-match score used while alignment
    score: F,
    /// Alignment of query sequence wrt reference
    alignment: Option<PairwiseAlignment>,
}

impl<'seq, F> AsRef<Self> for DiffStat<'seq, F>
where
    F: MatchFunc + Clone + Display,
{
    fn as_ref(&self) -> &Self {
        self
    }
}

impl<'seq, F> DiffStat<'seq, F>
where
    F: MatchFunc + Clone + Display,
{
    /// Generate a new [DiffStat] that aligns query with reference as master sequence
    pub fn new<T, G>(reference: &'seq T, query: &'seq T, gap_penalty: G, score: F) -> Self
    where
        T: AsRef<[u8]> + ?Sized,
        G: Into<GapPanelty>,
    {
        Self {
            reference: reference.as_ref(),
            query: query.as_ref(),
            alignment: None,
            gap_penalty: gap_penalty.into(),
            score,
        }
    }

    /// Calculate [Levenshtein](https://en.wikipedia.org/wiki/Levenshtein_distance) distance
    pub fn levenshtein(&self) -> u32 {
        tracing::info!("Calculating Lavenshtein distance");
        levenshtein(self.reference, self.query)
    }

    pub fn levenshtein_simd(&self) -> u32 {
        tracing::info!("Calculating Lavenshtein distance(simd)");
        bio::alignment::distance::simd::levenshtein(self.reference, self.query)
    }

    /// Calculate [Hamming](https://en.wikipedia.org/wiki/Hamming_distance) distance
    pub fn hamming_distance(&self) -> u64 {
        tracing::info!("Calculating Hamming distance");
        hamming(self.reference, self.query)
    }

    pub fn hamming_distance_simd(&self) -> u64 {
        tracing::info!("Calculating Hamming distance(simd)");
        bio::alignment::distance::simd::hamming(self.reference, self.query)
    }

    fn aligner(&self) -> bio::alignment::pairwise::Aligner<F>
    where
        F: MatchFunc,
    {
        bio::alignment::pairwise::Aligner::with_capacity(
            self.reference.len(),
            self.query.len(),
            self.gap_penalty.open,
            self.gap_penalty.extend,
            self.score.clone(),
        )
    }

    /// Pairwise alignment using Smith Waterman algorithm (Semiglobal)
    pub fn pairwise_aligner_semiglobal(&mut self) {
        tracing::info!(
            "Performing pairwise alignment (semiglobal) using {} and {}",
            self.gap_penalty,
            self.score
        );
        self.alignment = Some(self.aligner().semiglobal(self.reference, self.query));
    }

    /// Pairwise alignment using Smith Waterman algorithm (Global)
    pub fn pairwise_aligner_global(&mut self) {
        tracing::info!(
            "Performing pairwise alignment (global) using {} and {}",
            self.gap_penalty,
            self.score
        );
        self.alignment = Some(self.aligner().global(self.reference, self.query));
    }

    /// Pairwise alignment using Smith Waterman algorithm (Local)
    pub fn pairwise_aligner_local(&mut self) {
        tracing::info!(
            "Performing pairwise alignment (local) using {} and {}",
            self.gap_penalty,
            self.score
        );
        self.alignment = Some(self.aligner().local(self.reference, self.query));
    }

    /// CAUTION: Use for small sequence only, its running time complexity is
    /// `O(N^2 * L^2)`, where `N` is the number of sequences and `L` is the length of each sequence.
    ///
    /// Partial order alignment is not sustained in [DiffStat] struct.
    pub fn partial_order_alignment<S>(
        &self,
        scoring: Scoring<F>,
        references: Option<Vec<S>>,
    ) -> PartialorderAlignment
    where
        S: AsRef<[u8]>,
    {
        tracing::info!("Performing pairwise alignment (semiglobal)");
        let mut aligner = bio::alignment::poa::Aligner::new(scoring, self.reference);
        references.into_iter().flatten().for_each(|reference| {
            aligner.global(reference.as_ref()).add_to_graph();
        });
        aligner.global(self.query.as_ref()).alignment()
    }

    /// Pretty print the alignment, see [bio::alignment::Alignment::pretty]
    pub fn pretty_print(&self, coloumn: usize) {
        tracing::info!("Pretty print with {} coloumns", coloumn);
        if let Some(pretty) = self.alignment
            .as_ref()
            .map(|alignment| alignment.pretty(self.reference, self.query, coloumn)) { println!("{pretty}") }
    }

    pub fn alignment(&self) -> Option<&Alignment> {
        self.alignment.as_ref()
    }
}
