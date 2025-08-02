# Pericyte Hypoxia RNA-seq — Preprocessing & Analysis

This repository contains an end-to-end workflow for RNA-seq analysis of human iPSC-derived pericytes comparing **Hypoxia vs Normoxia** conditions.  
It includes two main components:

1. **Preprocessing** (`workflow/preprocess.sh`)  
   Downloads SRA runs, builds 4 sample FASTQs, performs QC, trimming, alignment, and gene-level counting.

2. **Analysis** (`workflow/HypoxiaRNAseq_analysis.Rmd`)  
   Performs differential expression analysis using DESeq2, generates PCA and volcano plots, and visualizes selected gene panels (as reported in the paper). 

---

## Data

- **SRA BioProject:** [`PRJNA1300358`](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1300358) (Homo sapiens, single-end RNA-seq; 4 biological samples formed by concatenating 4 runs each)
  - **Normoxia:** GSM9147330 (`SRR34830030–33`), GSM9147329 (`SRR34830034–37`)  
  - **Hypoxia:** GSM9147328 (`SRR34830038–41`), GSM9147327 (`SRR34830042–45`)

- **Reference files:**
  - **Transcript annotation:** GENCODE v44 GTF  
  - **Alignment index:** HISAT2 `grch38_tran`

**Related publication:**  
[PMID: 40737487](https://pubmed.ncbi.nlm.nih.gov/40737487/)

---

## Directory Output

The preprocessing script creates a structured working directory `data_pre_processing/` containing:

```
raw/      # downloaded and processed SRA FASTQs  
fastq/    # concatenated FASTQs (one per sample)  
trimmed/  # adapter-trimmed reads  
aligned/  # HISAT2-aligned BAMs  
counts/   # gene-level count matrix  
logs/     # logs from trimming, alignment, counting  
qc/       # FastQC output  
hisat2_index/  # reference index for alignment  
```

---

## Notes

- You can skip raw FASTQ processing and go straight to analysis using `data_RNAcounts/final_counts_symbols.tsv`.
- Compatible with macOS and Linux; no system-specific dependencies in the code beyond common bioinformatics tools.
- Code to generate reproducible answers to questions are in a separate .rmd file under the workflow folder
