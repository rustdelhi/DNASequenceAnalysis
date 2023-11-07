use std::time::Instant;

use aliner::DiffStat;
use bio::io::fasta::Reader;

use crate::aliner::Score;

mod aliner;
mod mutation_detection;

// TODO:
// 1. Comment aligner.rs
// 2. Make [Score] as different struct that scores match, unmatch, gap open, gap extend
// 3. Align all other fasta files and save them in accordance to a master refrence

fn main() {
    let covid_beta = Reader::from_file("./assets/SARS-beta.fasta").unwrap();
    let covid_delta = Reader::from_file("./assets/SARS-delta.fasta").unwrap();

    let beta_record = covid_beta.records().next().unwrap().unwrap();
    let beta_seq = beta_record.seq();
    let delta_record = covid_delta.records().next().unwrap().unwrap();
    let delta_seq = delta_record.seq();

    let diff = DiffStat::new(beta_seq, delta_seq);
    let score = Score::new(1, -1);
    let time = Instant::now();
    let distance = diff.pairwise_aligner((-5, -1).into(), score);

    println!("{:?}", distance);

    println!("time taken: {:?}", time.elapsed());
}
