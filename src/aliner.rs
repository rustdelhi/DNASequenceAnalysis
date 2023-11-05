//! This module is used to align two or more DNA/RNA sequences
//! to "align" them, see: https://en.wikipedia.org/wiki/Sequence_alignment

use bio::alignment::{
    distance::{hamming, levenshtein},
    pairwise::{MatchFunc, Scoring},
};

pub struct GapScore {
    pub open: i32,
    pub extend: i32,
}

impl GapScore {
    pub fn new(open: i32, extend: i32) -> Self {
        assert!(open < 0, "Gap open apnelty cant be positive");
        assert!(extend < 0, "Gap extend apnelty cant be positive");
        Self { open, extend }
    }
}

impl From<(i32, i32)> for GapScore {
    fn from(value: (i32, i32)) -> Self {
        Self::new(value.0, value.1)
    }
}

pub struct DiffStat<'seq> {
    reference: &'seq [u8],
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

    pub fn levenshtein(&self) -> u32 {
        levenshtein(self.reference, self.query)
    }

    /// [Self::levenshtein] with SIMD
    pub fn levenshtein_simd(&self) -> u32 {
        bio::alignment::distance::simd::levenshtein(self.reference, self.query)
    }

    pub fn hamming_distance(&self) -> u64 {
        hamming(self.reference, self.query)
    }

    pub fn hamming_distance_simd(&self) -> u64 {
        bio::alignment::distance::simd::hamming(self.reference, self.query)
    }

    pub fn pairwise_aligner<F>(&self, gap: GapScore, score: F) -> bio::alignment::Alignment
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
