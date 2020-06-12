# RNA-Seq_Tool

## Usage
```
./analysis.sh [mode] [input] [output] [log2(FoldChange) cutoff] [log10(padj) cutoff]

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
          
    [log2(FoldChange) cutoff] / [log10(padj) cutoff]:
        The customized cutoff when plotting volcano plot, should be provided as a positive number.
        Default is 4 and 2, respectively.
        
Sample:
    ./analysis.sh -all path/to/htseq-files/  path/to/outputDir/ 4 2
```
The script will first execute *combine_htseq.R* to merge all htseq-count files into one *htseq_counts_all.csv*, then it will run *analysis.R* to perform RNA-Seq analysis.

## Note
Run translate.R beforehand if desired genes are not found in the translation table. This process may take a long time as it needs to query [BioTools.fr](https://biotools.fr/human/ensembl_symbol_converter) if the ensembl isn't present in the downloaded HUGO database (hugo.txt).
