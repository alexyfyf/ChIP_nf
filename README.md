# ChIP_nf
Nextflow pipeline for ChIP-seq data analysis


# ChIP-seq pipeline

# General
Ref to ENCODE specification

# Details
## Pre-analysis
### pre-QC: 
fastqc
### trim (skipped): 
trimmomatic?
### alignment: 
bwa mem
### post processing: 
picard MarkDuplicates, 
picard CollectAlignmentSummaryMetrics,
picard CollectInsertSizeMetrics (paired-end)
samtools flagstat
samtools idxstats
samtools low quality removal, unmapped/unpaired/not proper paired removal

### post qc: 
fastqc
### generate summary metrics (in progress): 
% mapped, % chrM, % dup, % after all filtering

## Core analysis (in progress)
### peak calling: 
macs2 
Homer 
HMMRATAC

## Advanced analysis (in progress)
### peak anno and comparison: 
upsetplot, 
ChIPseeker
### Diff peak analysis:
DiffBind
### motif scan: 
FIMO, 
Homer


