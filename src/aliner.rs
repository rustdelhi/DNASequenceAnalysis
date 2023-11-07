//! This module is used to align two or more DNA/RNA sequences
//! to "align" them, see: https://en.wikipedia.org/wiki/Sequence_alignment

use bio::alignment::{
    distance::{hamming, levenshtein},
    pairwise::{MatchFunc, Scoring},
};

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
}

impl<'seq> DiffStat<'seq> {
    pub fn new<T>(reference: &'seq T, query: &'seq T) -> Self
    where
        T: AsRef<[u8]> + ?Sized,
    {
        Self {
            reference: reference.as_ref(),
            query: query.as_ref(),
        }
    }

    /// Calculate [Levenshtein](https://en.wikipedia.org/wiki/Levenshtein_distance) distance
    pub fn levenshtein(&self) -> u32 {
        levenshtein(self.reference, self.query)
    }

    pub fn levenshtein_simd(&self) -> u32 {
        bio::alignment::distance::simd::levenshtein(self.reference, self.query)
    }

    /// Calculate [Hamming](https://en.wikipedia.org/wiki/Hamming_distance) distance
    pub fn hamming_distance(&self) -> u64 {
        hamming(self.reference, self.query)
    }

    pub fn hamming_distance_simd(&self) -> u64 {
        bio::alignment::distance::simd::hamming(self.reference, self.query)
    }

    /// Pairwise alignment using Smith Waterman algorithm
    pub fn pairwise_aligner(&self, gap: GapPanelty, score: Score) -> bio::alignment::Alignment {
        bio::alignment::pairwise::Aligner::with_capacity(
            self.reference.len(),
            self.query.len(),
            gap.open,
            gap.extend,
            score,
        )
        .semiglobal(self.reference, self.query)
    }

    /// CAUTION: Use for small sequence only as its running time complexity is
    /// `O(N^2 * L^2)`, where `N` is the number of sequences and `L` is the length of each sequence.
    pub fn partial_order_alignment<F, S>(
        &self,
        scoring: Scoring<F>,
        references: Option<Vec<S>>,
    ) -> bio::alignment::poa::Alignment
    where
        F: MatchFunc,
        S: AsRef<[u8]>,
    {
        let mut aligner = bio::alignment::poa::Aligner::new(scoring, self.reference);
        references.into_iter().flatten().for_each(|reference| {
            aligner.global(reference.as_ref()).add_to_graph();
        });
        aligner.global(self.query.as_ref()).alignment()
    }
}
