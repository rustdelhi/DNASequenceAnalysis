use std::time::Instant;

use aliner::DiffStat;

use crate::{
    aliner::{GapPanelty, Score},
    reader::FastaReader,
};

mod aliner;
mod mutation_detection;
mod reader;

// TODO:
// 1. Mutation detection scoring (stats) [[mutation_detection.rs]]
// 2. Pretty print is not considering '\n' as newline
// 3. Add config file
//

// Initilize tracing
fn init() -> tracing_appender::non_blocking::WorkerGuard {
    let (non_blocking, guard) = tracing_appender::non_blocking(std::io::stdout());
    tracing_subscriber::fmt().with_writer(non_blocking).init();
    guard
}

fn main() -> anyhow::Result<()> {
    let _guard = init();
    let covid_beta = FastaReader::from_file("./assets/SARS-beta.fasta")?;
    let covid_delta = FastaReader::from_file("./assets/SARS-delta.fasta")?;
    let covid_omicron = FastaReader::from_file("./assets/SARS-omicron.fasta")?;
    let covid_gamma = FastaReader::from_file("./assets/SARS-gamma.fasta")?;
    let covid_zeta = FastaReader::from_file("./assets/SARS-zeta.fasta")?;

    let beta_record = covid_beta.records().next().unwrap()?;
    let beta_seq = beta_record.seq();
    let delta_record = covid_delta.records().next().unwrap()?;
    let delta_seq = delta_record.seq();

    let score = Score::new(1, -1);
    let gap = GapPanelty::new(-5, -1);

    let mut diff = DiffStat::new(beta_seq, delta_seq, gap, score);

    let time = Instant::now();
    diff.pairwise_aligner_semiglobal();

    println!("{:?}", diff.pretty(120).unwrap());

    println!("time taken: {:?}", time.elapsed());
    Ok(())
}
