use bio::alignment::Alignment;

use crate::aliner::DiffStat;

struct MutationStats {
    matches: usize,
    gaps: usize,
    miss_match: usize,
    subs: usize,
}

struct Muatation<'m> {
    diffstat: &'m DiffStat<'m>,
    pairwise_alignment: Option<Alignment>,
}

impl<'m> Muatation<'m> {
    pub fn from<D>(diffstat: &'m D) -> Self
    where
        D: AsRef<DiffStat<'m>>,
    {
        Self {
            diffstat: diffstat.as_ref(),
            pairwise_alignment: None,
        }
    }
}
