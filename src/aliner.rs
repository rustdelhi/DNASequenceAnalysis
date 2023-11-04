//! This module is used to align two or more DNA/RNA sequences
//! to "align" them, see: https://en.wikipedia.org/wiki/Sequence_alignment

use bio::alignment::{
    distance::{hamming, levenshtein},
    pairwise::MatchFunc,
    Alignment,
};

pub struct Gap {
    pub open: i32,
    pub extend: i32,
}

impl Gap {
    pub fn new(open: i32, extend: i32) -> Self {
        Self { open, extend }
    }
}

impl From<(i32, i32)> for Gap {
    fn from(value: (i32, i32)) -> Self {
        Self {
            open: value.0,
            extend: value.1,
        }
    }
}

pub struct DiffStat<'seq> {
    reference: &'seq [u8],
    query: &'seq [u8],
}

impl<'seq> DiffStat<'seq> {
    pub fn new<T>(reference: &'seq T, query: &'seq T) -> Self
    where
        T: AsRef<[u8]>,
    {
        Self {
            reference: reference.as_ref(),
            query: query.as_ref(),
        }
    }

    /// Use levenshtein method to find distance from start of reference to where the refrence sequence
    /// is matched in query tech
    /// see: https://en.wikipedia.org/wiki/Levenshtein_distance
    /// Complexity: O(refrence * query)
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

    pub fn pairwise_aligner<F>(&self, gap: Gap, score: F) -> Alignment
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
}
