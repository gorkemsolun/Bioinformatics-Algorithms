Görkem Kadir Solun 22003214 HW 1 CS481 kadir.solun@ug.bilkent.edu.tr

To prepare the executable: make
To clean the executable: make clean
To run the executable: ./hw1 -r <reference_fasta> -p <pattern_fasta> -o <output_prefix> -d

A complete solution for the multiple pattern matching assignment.
This implementation builds a single suffix tree that holds all reference sequences
(each separated by a unique terminator) and then searches for each pattern.

Command-line arguments to run:
    -r <reference_fasta> : FASTA file with reference sequences.
    -p <pattern_fasta>   : FASTA file with patterns.
    -o <output_prefix>   : The prefix for output files (a .txt file is always produced,
                                                 and if -d is specified, a DOT file is produced).
    -d                   : Optional; if provided, outputs the suffix tree in DOT format.

The implementation uses Ukkonen's algorithm for suffix tree construction.
It reports pattern occurrences per reference (using 0-indexed positions).

The node leaf naming logic is slightly modified as the following to be suitable for single suffix tree
<reference_name>:<1-based_start_index_in_that_reference>:<global_concat_index>