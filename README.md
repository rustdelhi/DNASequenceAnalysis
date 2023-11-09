# DNA Sequence Analysis in Rust

## Introduction

This Rust project aims to provide a robust and efficient solution for DNA sequence analysis,
including sequence alignments and mutation detection. Leveraging the power of Rust's performance
and safety features, this project is designed to handle large-scale DNA data for various
bioinformatics applications.

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

1. Adjust the parameters and sequences in the source code to fit your specific use case.

2. Run the project:

```bash
cargo run --release
```
