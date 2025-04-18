Center‑Star Multiple Sequence Alignment - Görkem Kadir Solun 22003214

This project implements a center‑star multiple sequence alignment (MSA) algorithm with affine‑gap penalties in C++. It reads a FASTA file of sequences, aligns them in PHYLIP format.

Requirements

C++17 (or later) compiler (e.g., g++, clang++).
Make utility.

Compilation

Simply run:

make

Usage

./hw3 -i <input.fasta> -o <output.phy> -s <M:Mm:Go:Ge>

Where:

-i <input.fasta> specifies the FASTA file containing multiple sequences.
-o <output.phy> specifies the output file in PHYLIP format.
-s <M:Mm:Go:Ge> sets scoring parameters:
M  = match score (e.g., 5)
Mm = mismatch score (e.g., -4)
Go = gap opening penalty (e.g., -16)
Ge = gap extension penalty (e.g., -4)

Example

./hw3 -i input.fasta -o output.phy -s 5:-4:-16:-4