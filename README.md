# 🧬 Genomics & Variant Calling Pipeline

This pipeline provides a comprehensive workflow for human genomic data analysis, focusing on Germline Variant Calling (SNPs and Indels) using the GRCh38 reference genome.

---

##  Table of Contents
1. [Data Acquisition & Reference Preparation](#1-reference-genome-setup)
2. [Quality Control & Preprocessing](#2-quality-control)
3. [Read Alignment](#3-read-alignment)
4. [Post-Alignment Processing](#4-alignment-post-processing)
5. [Variant Calling](#5-variant-calling)
6. [Variant Filtration](#6-variant-filtration)
7. [Variant Annotation](#7-variant-annotation)
---


## 1. Data Acquisition & Reference Preparation
### 1.1 Downloading the Reference Genome (GRCh38)
We use the GRCh38 "No-ALT" Analysis Set. This version is the industry standard because it removes redundant "alternate" sequences that confuse alignment tools, ensuring reads map to the correct chromosomes.

```Bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna reference.fa
```
### 1.2 Indexing the Reference
Create index files so that tools can jump to specific genomic coordinates:

```Bash
bwa index reference.fa
samtools faidx reference.fa
gatk CreateSequenceDictionary -R reference.fa
```
## 2. Quality Control & Preprocessing
### 2.1 Quality Check
Run FastQC or MultiQC to assess sequencing quality, checking for low-quality bases, adapter contamination, or low-quality cycles before starting the analysis:

```Bash
mkdir fastqc_out
fastqc *.fastq.gz -o fastqc_out/
multiqc fastqc_out/
```
### 2.2 Trimming 
Remove synthetic adapter sequences and cuts off low-quality bases from the ends of the reads to prevent false variant calls using Trimmomatic:

```Bash
trimmomatic PE -threads 8 sample_R1.fastq.gz sample_R2.fastq.gz \
    R1_paired.fq.gz R1_unpaired.fq.gz \
    R2_paired.fq.gz R2_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
## 3. Read Alignment
Align reads using BWA-MEM to map the processed reads to their original location on the GRCh38 human genome. The @RG (Read Group) info is added here to label the sample for GATK:

```Bash
bwa mem -t 8 -R "@RG\tID:sample01\tLB:lib1\tPL:ILLUMINA\tSM:sample01" \
    reference.fa R1_paired.fq.gz R2_paired.fq.gz > aligned.sam
```
## 4. Post-Alignment Processing
### 4.1 Sorting and Marking Duplicates
Convert SAM to BAM and sort by coordinate and filter for properly-paired reads (Q20):
* **-f 0x02:** Keeps only Properly Paired reads
* **-q 20:** Removes reads with low mapping quality

```Bash
# Align with BWA-MEM (including Read Group info)
bwa mem -t 8 -R "@RG\tID:sample1\tLB:lib1\tPL:ILLUMINA\tSM:sample1" \
    reference.fa R1_p.fq.gz R2_p.fq.gz > aligned.sam

# Filter and Convert to BAM
samtools view -h -f 0x02 -q 20 -b aligned.sam > aligned_pp.bam

# Sort and index the cleaned data
samtools sort aligned_pp.bam -o aligned_pp_sorted.bam
samtools index aligned_pp_sorted.bam
```
### 4.2 Marking Duplicates
Identify and flag PCR/optical duplicates to avoid false confidence in variant calling:

```Bash
gatk MarkDuplicates \
    -I aligned_pp_sorted.bam \
    -O aligned_marked.bam \
    -M metrics.txt \
    --CREATE_INDEX true
```
### 4.3 Base Quality Score Recalibration (BQSR)
Check known variant sites (like dbSNP) to identify and fix systematic errors made by the sequencing machine:

```Bash
# Analyze errors using multiple known resource files (from gatk)
gatk BaseRecalibrator \
    -I aligned_marked.bam \
    -R reference.fa \
    --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O recal_data.table

# 2. Apply corrections to the quality scores
gatk ApplyBQSR \
    -I aligned_marked.bam \
    -R reference.fa \
    --bqsr-recal-file recal_data.table \
    -O final_ready.bam
```
## 5. Variant Calling
Using GATK HaplotypeCaller to produce a Genomic VCF (gVCF) for individual samples:

```Bash
# Call Variants (HaplotypeCaller)
gatk HaplotypeCaller \
    -R reference.fa \
    -I final_ready.bam \
    -O output.g.vcf.gz \
    -ERC GVCF
```
## 6. Variant Filtration
Apply Hard Filtering based on GATK recommendations (QD, FS, SOR, MQ):
* **QD (Quality by Depth):** Low confidence despite high depth (< 2.0)
* **FS (Fisher Strand):** Evidence only on one DNA strand (> 60.0)
* **SOR (Strand Odds Ratio):** Mathematical bias at ends of reads (> 3.0)
* **MQ (Mapping Quality):** Reads map to multiple locations (< 40.0)

```Bash
#Convert gVCF to VCF
gatk GenotypeGVCFs \
    -R reference.fa \
    -V output.g.vcf.gz \
    -O raw_variants.vcf

# Filter
gatk VariantFiltration \
    -R reference.fa \
    -V raw_variants.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "filter_fail" \
    -O filtered_variants.vcf
```
## 7. Functional Annotation
Annotate filtered variants to determine their biological impact (e.g., missense, nonsense, intronic) using SnpEff or VEP:

```Bash
# Using VEP for GRCh38
vep -i filtered_variants.vcf \
    -o annotated.vcf \
    --assembly GRCh38 \
    --cache \
    --dir_cache /path/to/vep_cache \
    --refseq \
    --offline \
    --everything
```

## As an sbatch file in HPC

```bash
#!/bin/bash

#SBATCH --job-name=Human_Genomics_Pipeline
#SBATCH --output=genomics_pipeline_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00         
#SBATCH --partition=compute
#SBATCH --mail-user=your_email@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

# ---------------------------------------------------------
# 0. Load Required Modules
# (Update these names according to your HPC's module system)
# ---------------------------------------------------------
module load bwa/0.7.17
module load samtools/1.19
module load gatk/4.5.0
module load fastqc/0.11.9
module load trimmomatic/0.39
# module load vep/110 # Or use your conda environment

# ---------------------------------------------------------
# 1. Configuration & Variables
# ---------------------------------------------------------
REF="reference.fa"
KNOWN_DBSNP="Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_INDELS="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
CACHE_DIR="/path/to/vep_cache"
THREADS=$SLURM_CPUS_PER_TASK

echo "=========================================================="
echo "PIPELINE START TIME: $(date)"
echo "RUNNING ON NODE: $SLURMD_NODENAME"
echo "USING $THREADS CORES"
echo "=========================================================="

# ---------------------------------------------------------
# 2. Quality Control & Preprocessing
# ---------------------------------------------------------
echo "STEP 2.1: Running FastQC..."
mkdir -p fastqc_out
fastqc *.fastq.gz -o fastqc_out/
echo "FastQC finished."

echo "STEP 2.2: Trimming adapters with Trimmomatic..."
trimmomatic PE -threads $THREADS sample_R1.fastq.gz sample_R2.fastq.gz \
    R1_p.fq.gz R1_u.fq.gz R2_p.fq.gz R2_u.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
echo "Trimming finished."

# ---------------------------------------------------------
# 3. Read Alignment
# ---------------------------------------------------------
echo "STEP 3: Aligning reads with BWA-MEM..."
bwa mem -t $THREADS -R "@RG\tID:sample1\tLB:lib1\tPL:ILLUMINA\tSM:sample1" \
    $REF R1_p.fq.gz R2_p.fq.gz > aligned.sam
echo "Alignment finished."

# ---------------------------------------------------------
# 4. Post-Alignment Processing
# ---------------------------------------------------------
echo "STEP 4.1: Filtering for Properly Paired reads (MAPQ > 20)..."
samtools view -h -f 0x02 -q 20 -b aligned.sam > aligned_pp.bam
echo "Filtering finished."

echo "STEP 4.2: Sorting and Indexing BAM..."
samtools sort -@ $THREADS aligned_pp.bam -o aligned_pp_sorted.bam
samtools index aligned_pp_sorted.bam
echo "Sorting and indexing finished."

echo "STEP 4.3: Marking PCR Duplicates..."
gatk MarkDuplicates \
    -I aligned_pp_sorted.bam \
    -O aligned_marked.bam \
    -M metrics.txt \
    --CREATE_INDEX true
echo "Marking duplicates finished."

echo "STEP 4.4: Base Quality Score Recalibration (BQSR)..."
gatk BaseRecalibrator \
    -I aligned_marked.bam \
    -R $REF \
    --known-sites $KNOWN_DBSNP \
    --known-sites $KNOWN_INDELS \
    -O recal_data.table

gatk ApplyBQSR \
    -I aligned_marked.bam \
    -R $REF \
    --bqsr-recal-file recal_data.table \
    -O final_ready.bam
echo "BQSR finished."

# ---------------------------------------------------------
# 5. Variant Calling
# ---------------------------------------------------------
echo "STEP 5: Variant Calling with HaplotypeCaller (gVCF mode)..."
gatk HaplotypeCaller \
    -R $REF \
    -I final_ready.bam \
    -O output.g.vcf.gz \
    -ERC GVCF
echo "Variant calling finished."

# ---------------------------------------------------------
# 6. Variant Filtration
# ---------------------------------------------------------
echo "STEP 6.1: Genotyping gVCF to VCF..."
gatk GenotypeGVCFs \
    -R $REF \
    -V output.g.vcf.gz \
    -O raw_variants.vcf

echo "STEP 6.2: Applying Hard Filters..."
gatk VariantFiltration \
    -R $REF \
    -V raw_variants.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "filter_fail" \
    -O filtered_variants.vcf
echo "Filtration finished."

# ---------------------------------------------------------
# 7. Functional Annotation
# ---------------------------------------------------------
echo "STEP 7: Annotating variants with VEP..."
vep -i filtered_variants.vcf \
    -o annotated.vcf \
    --assembly GRCh38 \
    --cache \
    --dir_cache $CACHE_DIR \
    --refseq \
    --offline \
    --everything
echo "Annotation finished."

echo "=========================================================="
echo "PIPELINE FINISH TIME: $(date)"
echo "COMPLETED SUCCESSFULLY."
echo "=========================================================="
```
* Save the file as human_variant_calling.sh

* Submit the job: Use sbatch human_variant_calling.sh

* Monitor the progress: Use squeue -u your_username to see if it's running

* Check the logs: You can watch the echo messages in real-time by running: tail -f genomics_pipeline_*.log