#!/bin/bash

for i in src/rst_caller/test_seqs/*.fna; do
    echo "***"
    echo $i
    rst_caller -i $i -o test_out
    echo ""
    echo "Output File Content:"
    echo "RST_Type:"
    cat test_out/$(basename ${i%.*})_RST_TYPE.txt
    echo "Fasta Header:"
    grep ">" test_out/$(basename ${i%.*})_AMPLICON.fna
    echo "Fragment Lengths:"
    cat test_out/$(basename ${i%.*})_FRAGMENT_LENGTHS.txt
    echo ""
done