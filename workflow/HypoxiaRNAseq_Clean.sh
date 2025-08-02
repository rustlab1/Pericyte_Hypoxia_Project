#!/bin/bash

# RNA-seq mini-pipeline (SE reads): QC → trimming → alignment → counting

# -------------------- Setup --------------------


# Define project structure relative to current location
PROJECT_DIR="./data_pre_processing"
mkdir -p "${PROJECT_DIR}"/{raw,fastq,trimmed,aligned,counts,logs,qc,hisat2_index}


cd "${PROJECT_DIR}/raw"


# Group SRA run IDs by biological sample (4 runs each)
NORM1=(SRR34830030 SRR34830031 SRR34830032 SRR34830033)   # GSM9147330
NORM2=(SRR34830034 SRR34830035 SRR34830036 SRR34830037)   # GSM9147329
HYP1=(SRR34830038 SRR34830039 SRR34830040 SRR34830041)    # GSM9147328
HYP2=(SRR34830042 SRR34830043 SRR34830044 SRR34830045)    # GSM9147327

# -------------------- Download & Convert --------------------

# Download .sra files
for r in "${NORM1[@]}" "${NORM2[@]}" "${HYP1[@]}" "${HYP2[@]}"; do
  prefetch "$r"
done

# Convert to gzipped FASTQ
for r in "${NORM1[@]}" "${NORM2[@]}" "${HYP1[@]}" "${HYP2[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done

# Concatenate per-sample FASTQs
cat "${HYP1[@]/%/.fastq.gz}"  > Hyp1.fastq.gz
cat "${HYP2[@]/%/.fastq.gz}"  > Hyp2.fastq.gz
cat "${NORM1[@]/%/.fastq.gz}" > Norm1.fastq.gz
cat "${NORM2[@]/%/.fastq.gz}" > Norm2.fastq.gz

# Move to fastq/ folder
mv Hyp*.fastq.gz Norm*.fastq.gz ../fastq/

# -------------------- QC --------------------

cd ../fastq
fastqc Hyp1.fastq.gz Hyp2.fastq.gz Norm1.fastq.gz Norm2.fastq.gz \
  -o ../qc --threads 16

# -------------------- Trimming --------------------

cd ..
curl -L -o TruSeq3-SE.fa https://github.com/timflutre/trimmomatic/raw/master/adapters/TruSeq3-SE.fa

cd fastq
for sample in Hyp1 Hyp2 Norm1 Norm2
do
  trimmomatic SE -threads 16 -phred33 \
    ${sample}.fastq.gz ../trimmed/${sample}_trimmed.fastq.gz \
    ILLUMINACLIP:../TruSeq3-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    &> ../logs/${sample}_trimming.log
done

# -------------------- Alignment (HISAT2) --------------------

cd ../hisat2_index
curl -O ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz
tar -xzf grch38_tran.tar.gz

cd ../trimmed
for sample in Hyp1 Hyp2 Norm1 Norm2
do
  hisat2 -p 16 \
    -x ../hisat2_index/grch38_tran/genome_tran \
    -U ${sample}_trimmed.fastq.gz \
    2> ../logs/${sample}_hisat2.log | \
    samtools sort -@ 16 -o ../aligned/${sample}.bam
  samtools index ../aligned/${sample}.bam
done

# -------------------- Quantification (featureCounts) --------------------

cd ..
curl -L -o gencode.v44.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip -f gencode.v44.annotation.gtf.gz

featureCounts -T 16 -t exon -g gene_name \
  -a gencode.v44.annotation.gtf \
  -o counts/raw_counts_gene_sym.txt aligned/*.bam \
  &> logs/featureCounts_gene_sym.log

# Format counts matrix
{ printf "GeneSymbol\t"; head -n 2 counts/raw_counts_gene_sym.txt | tail -n 1 | cut -f7-; } > counts/final_counts_symbols.tsv
tail -n +3 counts/raw_counts_gene_sym.txt | \
  awk -v OFS="\t" '{ out=$1; for(i=7;i<=NF;i++) out=out OFS $i; print out }' >> counts/final_counts_symbols.tsv

sed -i '' '1 s|aligned/||g; 1 s|\.bam||g' counts/final_counts_symbols.tsv

# Done
echo "Pipeline complete. Output saved in: ${PROJECT_DIR}/counts/final_counts_symbols.tsv"
