# RNA-Seq_Tool

## Usage
```
./analysis.sh [mode] [input] [output]

Arguments:
    [mode]: 
        -all:     perform all analysis
        -deseq2:  perform only DESeq2 analysis
        -outrider perform only OUTRIDER analysis
    [input]:
        If a directory: path to a directory contaning HTSeq data
        If a file:      existing DESeq2/OUTRIDER rds file
    [output]:
        Directory to store program output
        The tool should produce:
          1. DESeq data set object (dds.rds) / OUTRIDER data set object (ods.rds)
          2. DESeq2_result.csv / OUTRIDER_result.csv
          3. Plot and graphs of analysis's result
```
The script will first execute *combine_htseq.R* to merge all htseq-count files into one *htseq_counts_all.csv*, then it will run *analysis.R* to perform RNA-Seq analysis.

## Note
Run translate.R beforehand if desired genes are not found in translation table. This process may take a long time as it will query for translation.
