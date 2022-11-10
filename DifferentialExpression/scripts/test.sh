#!/bin/bash

module load rmarkdown
module load deseqanalysis

  Rscript DESeq2.R \
    -i $PWD/../inputs/rsem \
    -o $PWD/../outputs/ \
    -m $PWD/../inputs/TEST_metadata.txt
