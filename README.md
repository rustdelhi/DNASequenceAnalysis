# DNA Sequence Analysis in Rust

## Introduction

This Rust project aims to provide a robust and efficient solution for basic DNA sequence analysis,
including sequence alignments and mutation detection. Leveraging the power of Rust's performance
and safety features.

## Features

1. _Pairwise Sequence Alignment_: Utilize the Smith Waterman algorithm for global,
   semiglobal, and local pairwise sequence alignments.

2. _Partial Order Alignment (POA)_: It is particularly useful when dealing with sequences that exhibit
   significant variations or insertions and deletions (indels). POA extends traditional pairwise alignment
   methods to handle multiple sequences simultaneously.

3. _Mutation Detection_: Detect mutations in DNA sequences, allowing for the identification
   of variations and differences between sequences.

## Getting Started

### Prerequisites

1. Rust: Ensure that Rust is installed on your system.
   If not, you can install it from https://www.rust-lang.org/tools/install.

### Installation

1. Clone the repository:

```bash
git clone https://github.com/rustdelhi/DNASequenceAnalysis.git
```

2. Navigate to the project directory:

```bash
cd DNASequenceAnalysis
```

3. Build the project:

```bash
cargo build --release
```

### Usage

1. Run the project with 2 FASTA files that contains only 1 sequence each
   (Example files included in assets folder):

```bash
cargo run --release -- --reference ./assets/SARS-beta.fasta --query ./assets/SARS-delta.fasta --print
```

Example output:

```
Score:
+---------+------------+--------------+------------+-----------+-------+
| match   | miss_match | substitution | insertions | deletions | total |
+---------+------------+--------------+------------+-----------+-------+
| 29707   | 156        | 85           | 13         | 58        | 30019 |
+---------+------------+--------------+------------+-----------+-------+
time taken: 11.037773246s
```

Run `cargo run --release -- --help` to know more about CLI usage

### Citations

1. [Sequence Alignmnt](https://en.wikipedia.org/wiki/Sequence_alignment)
2. [Smith Waterm Algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
3. [NCBI](https://www.ncbi.nlm.nih.gov/)
