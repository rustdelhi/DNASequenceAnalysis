//! This module is used to align two or more DNA/RNA sequences
//! to "align" them, see: https://en.wikipedia.org/wiki/Sequence_alignment

use bio::alignment::{
    distance::{hamming, levenshtein},
    pairwise::{MatchFunc, Scoring},
    Alignment,
};

type PairwiseAlignment = bio::alignment::Alignment;
type PartialorderAlignment = bio::alignment::poa::Alignment;

///  Scoring rule for [Substitution matrix](https://en.wikipedia.org/wiki/Smith_Waterman_algorithm#Substitution_matrix)
pub struct Score {
    r#match: i32,
    miss_match: i32,
}

impl Score {
    pub fn new(r#match: i32, miss_match: i32) -> Self {
        Self {
            r#match,
            miss_match,
        }
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
pub struct GapPanelty {
    pub open: i32,
    pub extend: i32,
}

impl GapPanelty {
    pub fn new(open: i32, extend: i32) -> Self {
        assert!(open < 0, "Gap open penalty cant be positive");
        assert!(extend < 0, "Gap extend penalty cant be positive");
        Self { open, extend }
    }
}

impl From<(i32, i32)> for GapPanelty {
    fn from(value: (i32, i32)) -> Self {
        Self::new(value.0, value.1)
    }
}

/// Compare two sequences and align them
pub struct DiffStat<'seq> {
    /// Master sequence
    reference: &'seq [u8],
    /// The one that will be aligned
    query: &'seq [u8],
    /// Alignment of query sequence wrt reference
    alignment: Option<PairwiseAlignment>,
}

impl<'seq> DiffStat<'seq> {
    /// Generate a new [DiffStat] that aligns query with reference as master sequence
    pub fn new<T>(reference: &'seq T, query: &'seq T) -> Self
    where
        T: AsRef<[u8]> + ?Sized,
    {
        Self {
            reference: reference.as_ref(),
            query: query.as_ref(),
            alignment: None,
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

    fn aligner<F>(&self, gap: GapPanelty, score: F) -> bio::alignment::pairwise::Aligner<F>
    where
        F: MatchFunc,
    {
        bio::alignment::pairwise::Aligner::with_capacity(
            self.reference.len(),
            self.query.len(),
            gap.open,
            gap.extend,
            score,
        )
    }

    /// Pairwise alignment using Smith Waterman algorithm (Semiglobal)
    pub fn pairwise_aligner_semiglobal(&mut self, gap: GapPanelty, score: Score) {
        tracing::info!("Performing pairwise alignment (semiglobal)");
        self.alignment = Some(
            self.aligner(gap, score)
                .semiglobal(self.reference, self.query),
        );
    }

    /// Pairwise alignment using Smith Waterman algorithm (Global)
    pub fn pairwise_aligner_global(&mut self, gap: GapPanelty, score: Score) {
        tracing::info!("Performing pairwise alignment (global)");
        self.alignment = Some(self.aligner(gap, score).global(self.reference, self.query));
    }

    /// Pairwise alignment using Smith Waterman algorithm (Local)
    pub fn pairwise_aligner_local(&mut self, gap: GapPanelty, score: Score) {
        tracing::info!("Performing pairwise alignment (local)");
        self.alignment = Some(self.aligner(gap, score).local(self.reference, self.query));
    }

    /// CAUTION: Use for small sequence only, its running time complexity is
    /// `O(N^2 * L^2)`, where `N` is the number of sequences and `L` is the length of each sequence.
    ///
    /// Partial order alignment is not sustained in [DiffStat] struct.
    pub fn partial_order_alignment<F, S>(
        &self,
        scoring: Scoring<F>,
        references: Option<Vec<S>>,
    ) -> PartialorderAlignment
    where
        F: MatchFunc,
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
    pub fn pretty(&self, coloumn: usize) -> Option<String> {
        tracing::info!("Pretty print with {} coloumns", coloumn);
        self.alignment
            .as_ref()
            .map(|alignment| alignment.pretty(self.reference, self.query, coloumn))
    }

    pub fn alignment(&self) -> Option<&Alignment> {
        self.alignment.as_ref()
    }
}
