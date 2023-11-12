
use dna_sequence_analysis::{
    aliner::{DiffStat, GapPanelty, Score},
    init_logging,
    mutation_detection::Muatation,
    reader::FastaReader,
};
use std::time::Instant;

fn main() -> anyhow::Result<()> {
    let _guard = init_logging();
    let covid_beta = FastaReader::from_file("./assets/SARS-beta.fasta")?;
    let covid_delta = FastaReader::from_file("./assets/SARS-delta.fasta")?;

    let beta_record = covid_beta.records().next().unwrap()?;
    let beta_seq = beta_record.seq();
    let delta_record = covid_delta.records().next().unwrap()?;
    let delta_seq = delta_record.seq();

    let score = Score::new(1, -1);
    let gap = GapPanelty::new(-5, -1);

    let mut diff = DiffStat::new(beta_seq, delta_seq, gap, score);

    let time = Instant::now();
    diff.pairwise_aligner_semiglobal();

    diff.pretty_print(120);
    let ms = Muatation::from(&diff);

    println!("Score: {}", ms.mutastion_score().unwrap());
    println!("time taken: {:?}", time.elapsed());
    Ok(())
}
