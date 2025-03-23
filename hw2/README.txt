Görkem Kadir Solun 22003214 CS481 HW2

This project performs sequence alignment (global using Needleman–Wunsch and local using Smith–Waterman) 
between pairs of sequences read from FASTA files. The program computes the alignment score, 
and constructs the CIGAR and MD:Z strings using a traceback of the alignment matrix. 
The pair with the longest overlap (for global alignment) or highest alignment score (for local alignment)
is written to the output file.

Algorithms are implemented using the course slides and Wikipedia.

Usage:
  make
  ./hw2 -g|-l -p <patterns.fasta> -t <texts.fasta> -o <output.txt> -s <match> <mismatch> <gap>

Example:
  ./hw2 -g -p patterns.fasta -t texts.fasta -o global.txt -s 1 -1 -1
