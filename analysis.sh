#!/bin/bash
time Rscript combine_htseq.R $1
#time Rscript translate.R $1
time Rscript analysis.R $1 $2