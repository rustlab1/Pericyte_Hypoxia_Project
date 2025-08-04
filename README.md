## Pericyte hypoxia RNA seq preprocessing and analysis

paper_title: Brain pericytes derived from human pluripotent stem cells retain vascular and phagocytic functions under hypoxia
doi: 10.1101/2025.04.10.648232
data_source: "GEO GSE304315; BioProject PRJNA1300358; SRA runs SRR34830030 to SRR34830045"


This repository provides an end to end workflow for RNA seq analysis of human iPSC derived pericytes comparing hypoxia and control.
It includes two parts

Preprocessing in shell
downloads SRA runs, builds sample FASTQs, runs QC and trimming, aligns to GRCh38, and produces gene level counts

Analysis in R Markdown
performs differential expression with DESeq2, PCA and volcano plots, a heatmap, and selected gene expression comparisons

#1. Data sources
Public GEO GSE304315 and BioProject PRJNA1300358 human iPSC derived pericyte single end RNA seq

Groups and runs used

Hypoxia SRR34830038 SRR34830039 SRR34830040 SRR34830041 and SRR34830042 SRR34830043 SRR34830044 SRR34830045

Normoxia SRR34830030 SRR34830031 SRR34830032 SRR34830033 and SRR34830034 SRR34830035 SRR34830036 SRR34830037

Reference

Genome GRCh38 HISAT2 grch38_tran index

Annotation GENCODE v44 GTF

#2. How to download
Install the tools with conda and fetch FASTQs from SRA using prefetch and fasterq dump. See workflow scripts for exact commands.

#3. Pre processing and subsampling
No subsampling is applied by default. Raw FASTQs are used as is. The dataset includes hypoxia and normoxia pericyte samples.
Quality control is performed with FastQC

#4. How the workflow works
The shell script executes the following steps. Save your script and run it from a terminal. It creates and uses ./data_pre_processing.

Step 1 Quality control

Purpose check base quality and GC content
Tools FastQC and optional MultiQC
Inputs FASTQs in fastq
Outputs per sample HTML and zip in qc

Step 2 Trimming

Purpose remove adapters and low quality bases
Tools Trimmomatic with TruSeq3 SE adapters
Inputs FASTQs in fastq
Outputs trimmed FASTQs in trimmed and trimming logs in logs

Step 3 Reference setup

Purpose obtain the GRCh38 HISAT2 transcriptome aware index
Tools curl and tar
Outputs index in hisat2_index/grch38_tran

Step 4 Alignment

Purpose align single end reads to GRCh38 and sort BAM
Tools HISAT2 and samtools
Inputs trimmed FASTQs in trimmed
Outputs BAM and BAI in aligned plus logs in logs

Step 5 Post alignment QC

Purpose review mapping statistics
Tools samtools

Step 6 Quantification

Purpose count reads per gene using the GENCODE GTF
Tools featureCounts
Inputs BAM files and GENCODE v44 GTF
Outputs counts in counts and a final matrix data_pre_processing/counts/final_counts_symbols.tsv

Step 7 Analysis

Purpose differential expression and figures
Tools R DESeq2 ggplot2 pheatmap
Inputs featureCounts output and sample metadata
Outputs PCA and volcano plots, QC tables, and selected gene comparisons

In R

run workflow/HypoxiaRNAseq_analysis_clean.Rmd to load data and perform the main analysis

Notes
You can skip raw FASTQ processing and start from the counts matrix if you already have one.