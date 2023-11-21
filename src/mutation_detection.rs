use std::fmt::Display;

use bio::alignment::pairwise::MatchFunc;
use tabled::{Table, Tabled};

use crate::aliner::DiffStat;

#[derive(Debug, Default, Tabled)]
pub struct MutationStats {
    r#match: usize,
    miss_match: usize,
    substitution: usize,
    insertions: usize,
    deletions: usize,
    total: usize,
}

impl Display for MutationStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", Table::new(vec![self]))
    }
}

impl MutationStats {
    pub fn inc_match(&mut self) {
        self.r#match += 1;
        self.total += 1
    }

    pub fn inc_miss_match(&mut self) {
        self.miss_match += 1;
        self.total += 1
    }

    pub fn inc_substitution(&mut self) {
        self.substitution += 1;
        self.inc_miss_match();
        self.total += 1
    }

    pub fn inc_insertions(&mut self) {
        self.insertions += 1;
        self.inc_miss_match();
        self.total += 1
    }

    pub fn inc_deletions(&mut self) {
        self.deletions += 1;
        self.inc_miss_match();
        self.total += 1;
    }
}

#[derive(Debug)]
pub struct Muatation<'m, F>
where
    F: MatchFunc + Clone + Display,
{
    diffstat: &'m DiffStat<'m, F>,
}

impl<'m, F> Muatation<'m, F>
where
    F: MatchFunc + Clone + Display,
{
    pub fn from<D>(diffstat: &'m D) -> Self
    where
        D: AsRef<DiffStat<'m, F>>,
        F: MatchFunc,
    {
        assert!(
            diffstat.as_ref().alignment().is_some(),
            "DiffStat is not aligned, please use pairwise alignment on Diffstast before using Muatation::from()"
        );
        Self {
            diffstat: diffstat.as_ref(),
        }
    }

    pub fn mutastion_score(&self) -> Option<MutationStats> {
        tracing::info!("Calcualting mutation score");
        self.diffstat.alignment().map(|alignment| {
            alignment
                .operations
                .iter()
                .fold(MutationStats::default(), |mut ms, operation| {
                    match operation {
                        bio::alignment::AlignmentOperation::Match => ms.inc_match(),
                        bio::alignment::AlignmentOperation::Subst => ms.inc_substitution(),
                        bio::alignment::AlignmentOperation::Del => ms.inc_deletions(),
                        bio::alignment::AlignmentOperation::Ins => ms.inc_insertions(),
                        _ => (),
                    }
                    ms
                })
        })
    }
}

#[cfg(test)]
mod test {
    use crate::aliner::{DiffStat, Score};

    use super::Muatation;

    #[test]
    #[should_panic]
    fn make_mutation_without_alignemnt() {
        let diffstat = DiffStat::new(&[], &[], (-1, -1), Into::<Score>::into((1, -1)));
        let _md = Muatation::from(diffstat.as_ref());
    }

    #[test]
    fn mutation_score_accuracy() {
        // Refer: https://docs.rs/bio/1.4.0/bio/alignment/struct.Alignment.html#method.pretty
        let mut diffstat = DiffStat::new(
            "CCGTCCGGCAAGGG",
            "AAAAACCGTTGACGGCCAA",
            (-1, -1),
            Into::<Score>::into((1, -1)),
        );
        diffstat.pairwise_aligner_global();
        let pretty = diffstat.pretty_string(120).expect("Unable to pretty print");

        #[rustfmt::skip]
        // Do not format this code, as raw strings are weird, and any unwanted spaces in expected
        // string will result in falure of this test
        //
        // This is expected string(3 trailing newlines are added intentionally):
        // -----CCGT--CCGGCAAGGG
        // xxxxx||||xx\||||\|++\
        // AAAAACCGTTGACGGCCA--A
        //
        //
        //
        let expected = r#"-----CCGT--CCGGCAAGGG
xxxxx||||xx\||||\|++\
AAAAACCGTTGACGGCCA--A


"#;
        assert_eq!(pretty, expected);
    }
}
