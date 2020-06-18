#!/bin/bash
FC="${4:-4}"
padj="${5:-2}"

Rscript combine_htseq.R $2
#Rscript translate.R $1
Rscript analysis.R $1 $2 $3 $FC $padj