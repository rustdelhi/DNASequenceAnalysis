use clap::{arg, Parser};
use dna_sequence_analysis::{
    aliner::{DiffStat, GapPanelty, Score},
    init_logging,
    mutation_detection::Muatation,
    reader::FastaReader,
};
use std::{path::PathBuf, time::Instant};

enum _CliError {
    Generic(String),
}

#[derive(Parser, Debug)]
struct Cli {
    /// Reference (master) FASTA file
    #[arg(short, long, value_name = "FILE")]
    reference: PathBuf,

    /// Query FASTA file, the sequence which will be aligned
    #[arg(short, long, value_name = "FILE")]
    query: PathBuf,

    /// Print the alignemnt
    #[arg(short, long)]
    print: bool,
}

fn find_mutation(reference: PathBuf, query: PathBuf, print: bool) -> anyhow::Result<()> {
    let reference = FastaReader::from_file(reference)?;
    let query = FastaReader::from_file(query)?;

    // FASTA files contain only 1 sequence
    let reference_record = reference.records().next().unwrap()?;
    let reference_seq = reference_record.seq();
    let query_record = query.records().next().unwrap()?;
    let query_seq = query_record.seq();

    // Default score
    let score = Score::new(1, -1);
    let gap = GapPanelty::new(-5, -1);

    let mut diff = DiffStat::new(reference_seq, query_seq, gap, score);

    let time = Instant::now();
    diff.pairwise_aligner_semiglobal();

    print.then(|| diff.pretty_print(120));

    let ms = Muatation::from(&diff);

    println!("Score: {}", ms.mutastion_score().unwrap());
    println!("time taken: {:?}", time.elapsed());

    Ok(())
}

fn run_cli() -> anyhow::Result<()> {
    let args = Cli::parse();

    find_mutation(args.reference, args.query, args.print)?;
    Ok(())
}

fn main() -> anyhow::Result<()> {
    let _ = init_logging();
    run_cli()?;

    Ok(())
}
