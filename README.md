# ChIP_nf
Nextflow pipeline for ChIP-seq data analysis

# Table of contents

- [General](#general)
- [Details](#details)
  - [Pre-analysis](#pre-analysis)
    - [pre-QC:](#pre-qc)
    - [trim (skipped):](#trim-skipped)
    - [alignment:](#alignment)
    - [post processing:](#post-processing)
    - [post qc:](#post-qc)
    - [generate summary metrics (in progress):](#generate-summary-metrics-in-progress)
  - [Core analysis (in progress)](#core-analysis-in-progress)
    - [peak calling:](#peak-calling)
  - [Advanced analysis (in progress)](#advanced-analysis-in-progress)
    - [peak anno and comparison:](#peak-anno-and-comparison)
    - [Diff peak analysis:](#diff-peak-analysis)
    - [motif scan:](#motif-scan)


# General
Ref to [ENCODE specification](https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#)

# Details
## Pre-analysis
### pre-QC: 
fastqc
### trim (skipped): 
trimmomatic?
### alignment: 
bwa mem
### post processing: 
picard MarkDuplicates  
picard CollectAlignmentSummaryMetrics  
picard CollectInsertSizeMetrics (paired-end)  
samtools flagstat  
samtools idxstats  
samtools low quality removal, unmapped/unpaired/not proper paired removal  
calculate PBC1, PBC2, NRF

### post qc: 
fastqc
### generate summary metrics (in progress): 
% mapped, % chrM, % dup, % after all filtering

## Core analysis (in progress)
### peak calling: 
macs2  
Homer  

## Advanced analysis (in progress)
### peak anno and comparison: 
upsetplot  
ChIPseeker
### Diff peak analysis:
DiffBind
### motif scan: 
FIMO, 
Homer


