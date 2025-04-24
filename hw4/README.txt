UPGMA Phylogenetic Tree Construction - Görkem Kadir Solun (22003214) - Bioinformatics Algorithms CS481

This program implements the UPGMA algorithm to construct a phylogenetic tree from DNA sequences provided in a FASTA file. 
It uses the Needleman-Wunsch algorithm to compute pairwise alignment distances, then builds a distance matrix, clusters sequences, 
and outputs the resulting tree in Newick format.

Requirements:
C++17 compiler (e.g., g++)
Make

Compilation:
make

Usage:
./hw4 -i <input.fasta> -t <output_tree.txt> -s   
Options:
-i <input.fasta>     Input FASTA file containing DNA sequences.
-t <output_tree.txt> Output file for the Newick-formatted UPGMA tree.
-s   
Scoring parameters: match score, mismatch penalty, gap penalty.

Example:
./hw4 -i sequences.fasta -t upgma_tree.txt -s 1 -1 -1

Clean:
make clean
