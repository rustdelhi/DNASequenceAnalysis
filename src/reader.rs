use std::{fs::File, io::BufReader, path::Path};

use bio::io::fasta::Records;

#[derive(Debug, thiserror::Error)]
pub enum FastaReaderError {
    #[error("Error: {0}")]
    Generic(String),
}

#[derive(Debug)]
pub struct FastaReader {
    inner: Records<BufReader<File>>,
}

impl FastaReader {
    pub fn from_file<P>(file_path: P) -> Result<Self, FastaReaderError>
    where
        P: AsRef<Path>,
    {
        let fasta_reader = bio::io::fasta::Reader::from_file(file_path.as_ref())
            .map_err(|err| FastaReaderError::Generic(err.to_string()))?;
        Ok(Self {
            inner: fasta_reader.records(),
        })
    }

    pub fn records(self) -> Records<BufReader<File>> {
        self.inner
    }
}
