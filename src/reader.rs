use std::{fmt::Display, fs::File, io::BufReader, path::Path};

use bio::io::fasta::{Record, Records};

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
        P: AsRef<Path> + Display,
    {
        tracing::info!("Fasta reader for file {}", file_path.to_string());
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

impl IntoIterator for FastaReader {
    type Item = Record;

    type IntoIter = FastaReaderIter;

    fn into_iter(self) -> Self::IntoIter {
        todo!()
    }
}

pub struct FastaReaderIter {
    inner: Records<BufReader<File>>,
}

impl Iterator for FastaReaderIter {
    type Item = Record;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().and_then(|rec| rec.ok())
    }
}
