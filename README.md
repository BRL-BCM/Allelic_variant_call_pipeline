# Allelic_variant_call_pipeline
Pipeline for allele-specific variant calling from RNAseq. This pipeline has been optimized to work with ERCC exRNA samples on BRL-BCM HPC.

## Prerequisite tools and databases
1. Construction of the diploid genomes requires the following dependencies:
  - gnomAD frequency database hg38_gnomad30_genome.txt available at [Annovar Main Package](https://annovar.openbioinformatics.org/en/latest/user-guide/download/#annovar-main-package)
  - [vcf2diploid](https://github.com/abyzovlab/vcf2diploid)
  - [Crossmap](https://github.com/liguowang/CrossMap)

2. The variant calling pipeline requires the following dependencies:
  - STAR, gatk, picard, samtools, BEDTools
  - Hg38 reference file,reference annotation file and known variants sites files such as 1000G_phase1.snps.high_confidence.hg38.vcf.gz, Mills_and_1000G_gold_standard.indels.hg38.vcf.gz, all available at [GATK resource bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).

## Running the pipeline
### Step 1: Keep only genome regions covered by at least one ERCC samples
1. getBedCoveragePerStudy.pbs.txt: generate unioned coverage depth bedfile across ERCC studies
2. getSumBedCoveragePerStudy.pbs.txt: utilize unionbedgByStudy.R to sumarize covered regions by study
3. 

Convert unioned ERCC genome coverage bedgraphs from version hg19 to hg38
### Step 2: Construct (artificial) diploid genome
### Step 3: Call extracellularly expressed variants
