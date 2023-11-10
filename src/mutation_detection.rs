use bio::alignment::{pairwise::MatchFunc, Alignment};

use crate::aliner::DiffStat;

struct MutationStats {
    matches: usize,
    gaps: usize,
    miss_match: usize,
    subs: usize,
}

struct Muatation<'m, F>
where
    F: MatchFunc + Clone,
{
    diffstat: &'m DiffStat<'m, F>,
    pairwise_alignment: Option<Alignment>,
}

impl<'m, F> Muatation<'m, F>
where
    F: MatchFunc + Clone,
{
    pub fn from<D>(diffstat: &'m D) -> Self
    where
        D: AsRef<DiffStat<'m, F>>,
        F: MatchFunc,
    {
        Self {
            diffstat: diffstat.as_ref(),
            pairwise_alignment: None,
        }
    }
}
