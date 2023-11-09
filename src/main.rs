use std::time::Instant;

use aliner::DiffStat;
use bio::io::fasta::Reader;

use crate::aliner::Score;

mod aliner;
mod mutation_detection;

// TODO:
// 1. Mutation detection scoring (stats) [[mutation_detection.rs]]
// 2. Pretty print is not considering '\n' as newline
// 3. Fasta struct to abstract reading from file
//

// Initilize tracing
fn init() -> tracing_appender::non_blocking::WorkerGuard {
    let (non_blocking, guard) = tracing_appender::non_blocking(std::io::stdout());
    tracing_subscriber::fmt().with_writer(non_blocking).init();
    guard
}

fn main() {
    let guard = init();
    let covid_beta = Reader::from_file("./assets/SARS-beta.fasta").unwrap();
    let covid_delta = Reader::from_file("./assets/SARS-delta.fasta").unwrap();

    let beta_record = covid_beta.records().next().unwrap().unwrap();
    let beta_seq = beta_record.seq();
    let delta_record = covid_delta.records().next().unwrap().unwrap();
    let delta_seq = delta_record.seq();

    let mut diff = DiffStat::new(beta_seq, delta_seq);
    let score = Score::new(1, -1);
    let time = Instant::now();
    diff.pairwise_aligner_semiglobal((-5, -1).into(), score);

    println!("{:?}", diff.pretty(120).unwrap());

    println!("time taken: {:?}", time.elapsed());
}
