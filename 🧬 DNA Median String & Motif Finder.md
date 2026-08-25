# 🧬 DNA Median String & Motif Finder

A C++ bioinformatics project that analyzes DNA sequences from a FASTA file to find the **top 5 median strings** and their corresponding **consensus motifs**.

## 🚀 Features

- 🧬 Reads multiple DNA sequences from a FASTA file
- 🔍 Generates all possible L-length DNA strings (L-mers)
- 📏 Calculates Hamming distance between DNA sequences
- 🏆 Finds the top 5 median strings with the lowest total distance
- 🧩 Builds a consensus motif for each median string
- 📍 Reports motif positions and matching substrings
- 📊 Calculates total-distance and consensus scores
- 🛡️ Handles invalid FASTA characters by replacing them with `N`

## 🧠 Algorithms Used

### 📏 Hamming Distance
Measures how many positions are different between two DNA strings.

Example:

`ACGT`  
`ACCT`

The Hamming distance is **1** because only one position is different.

### 🔬 Brute-Force L-mer Search
The program generates every possible DNA string of length `L` using:

`A`, `C`, `G`, `T`

For example, if `L = 2`:

`AA, AC, AG, AT, CA, CC ...`

Each possible string is compared against the DNA sequences to find the best median strings.

### 🧬 Consensus Motif
After finding a median string, the program examines the best matching regions from the DNA sequences and selects the most common nucleotide at each position to create a consensus sequence.

## 🛠️ Technologies

- 💻 C++
- 🧬 Bioinformatics Algorithms
- 📂 FASTA File Processing
- 📏 Hamming Distance
- 🔍 Brute-Force Search
- 📊 Sequence Analysis

## ⚙️ Compile

```bash
g++ -o median_finder Prog3.cpp
```

## ▶️ Run

```bash
./median_finder HMP-part.fa 6 results.txt
```

### Arguments

- `HMP-part.fa` — Input FASTA file
- `6` — Motif length (`L`)
- `results.txt` — Output file

## 📤 Output

The program writes:

- 🏆 Top 5 median strings
- 📏 Total distance for each median string
- 🧬 Consensus motif strings
- 📊 Consensus scores
- 📍 Motif positions and matching DNA substrings

## 🎯 Project Purpose

This project demonstrates the use of **C++ algorithms and data structures to solve a bioinformatics problem**. It combines file processing, exhaustive search, DNA sequence comparison, motif discovery, and consensus analysis.

## 👨‍💻 Skills Demonstrated

`C++` • `Algorithms` • `Bioinformatics` • `File I/O` • `String Processing` • `Data Structures` • `Brute-Force Search` • `DNA Sequence Analysis`