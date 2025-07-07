#!bin/bash 

#GATK Germline Variant Calling Training Project (HG00096)-Germline variant calling-SRR070802
#sample information link: https://www.ebi.ac.uk/ena/browser/view/SRR070802
#https://gatk.broadinstitute.org


#-------------------------------------------------------------------------

# download raw data --> paired end R1,R2
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070802/SRR070802_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070802/SRR070802_2.fastq.gz


echo "Run files..."

#------------------------------------------------------------------------

# download refrence 
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz 


# Create .fai and .dict for FASTA to enable fast access and ensure GATK-compatible reference metadata

#index of refrence file
samtools faix /Users/azadehnikouee/Desktop/GATK_germline_variant_calling/reference/hg38.fa


#dictionary of refrence file 
gatk CreateSequenceDictionary R=/Users/azadehnikouee/Desktop/GATK_germline_variant_calling/reference/hg38.fa   o=/Users/azadehnikouee/Desktop/GATK_germline_variant_calling/reference/hg38.dict


#known SNP sites (VCF + index) used by GATK BaseRecalibrator to avoid penalizing real variants during BQSR
#download known sites files for BQSR from GATK resource bundle

wget -P /Users/azadehnikouee/Desktop/GATK_germline_variant_calling/reference  https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /Users/azadehnikouee/Desktop/GATK_germline_variant_calling/reference   https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

#============================================================ End of prepare needed data to start variant calling


#directories :
ref="~/Desktop/GATK_germline_variant_calling/reference/hg38.fa"
reads="~/Desktop/GATK_germline_variant_calling/reads"
FASTQC="~/Desktop/GATK_germline_variant_calling/FASTQC"
aligned_reads="~/Desktop/GATK_germline_variant_calling/aligned_read"
known_var="~/Desktop/GATK_germline_variant_calling/known_var/Homo_sapiens_assembly38.dbsnp138.vcf"
results="~/Desktop/GATK_germline_variant_calling/result"

#=============================== 1- quality control or QC
#Input files = two fastq files (R1 ,R2) --> output HTML of QC summary results     --> RUN: fastqc

fastqc ${reads}/SRR070802_1.fastq.gz  -o ${FASTQC}/
fastqc ${reads}/SRR070802_2.fastq.gz  -o ${FASTQC}/

#quality good enough for continue and no need triming or using fastp

#=============================2-Map to reference --> run BWA-MEM
#make index from refrence , then BWA alignement , -t = number of threads -R= add groups record
sample_id="SRR070802"
rg="@RG\tID:${sample_id}\tPL:ILLUMINA\tSM:${sample_id}" 

# ===Index the reference (only needed once) ===
bwa index ${ref}


# === Run alignment ===
bwa mem -t 4 -R "${rg}" ${ref} \
  ${reads}/${sample_id}_1.fastq.gz \
  ${reads}/${sample_id}_2.fastq.gz \
  > ${aligned_reads}/${sample_id}.sam


echo "Alignment complete"

#================================= 3-Mark Duplicates and Sort bam files


# == pathway and folders==

aligned_reads="/Users/azadehnikouee/Desktop/GATK_germline_variant_calling/aligned_read"
logs_dir="/Users/azadehnikouee/Desktop/GATK_germline_variant_calling/logs"


# === Create output dir if needed ===
mkdir -p ${logs_dir}


# === Define file paths ===
input_sam="${aligned_reads}/SRR070802.sam"
output_bam="${aligned_reads}/SRR070802_sorted_dedup_reads.bam"
metrics_file="${logs_dir}/SRR070802.markdup.metrics.txt"



# === Run MarkDuplicatesSpark ===
gatk MarkDuplicatesSpark \
  -I ${input_sam} \
  -O ${output_bam} \
  --metrics-file ${metrics_file}


echo"MarkDuplicatesSpark completed. Output: ${output_bam}"


# After alignment read get Sam file theb check by samtools for some quality of this:

## samtools view *.sam | less 
## samtools flagstat *.sam 



#============================== 4-Base quality recalibration


#A-Build the recalibration model using known SNPs
gatk BaseRecalibrator \
  -I ${aligned_reads}/SRR070802_sorted_dedup_reads.bam \
  -R ${ref} \
  --known-sites ${known_var} \
  -O ${results}/recal_data.table

# B-Apply the BQSR model to adjust base quality scores
gatk ApplyBQSR \
  -I ${aligned_reads}/SRR070802_sorted_dedup_reads.bam \
  -R ${ref} \
  --bqsr-recal-file ${results}/recal_data.table \
  -O ${aligned_reads}/SRR070802_sorted_dedup_bqsr_reads.bam



#============================== 5-Collect Alignment & Insert Size Metrics

# A - Collect general alignment summary metrics
gatk CollectAlignmentSummaryMetrics \
  R=${ref} \
  I=${aligned_reads}/SRR070802_sorted_dedup_bqsr_reads.bam \
  O=${aligned_reads}/alignment_metrics.txt

# B - Collect insert size metrics and generate a histogram
gatk CollectInsertSizeMetrics \
  INPUT=${aligned_reads}/SRR070802_sorted_dedup_bqsr_reads.bam \
  OUTPUT=${aligned_reads}/insert_size_metrics.txt \
  HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf


#====================== 6- variant calling by haplotype caller (GATK)

# Call raw variants (VCF) using HaplotypeCaller
gatk HaplotypeCaller \
  -R ${ref} \
  -I ${aligned_reads}/SRR070802_sorted_dedup_bqsr_reads.bam \
  -O ${results}/raw_variants.vcf

echo "raw VCF completed

#====================== 7- EXtract SNP and Indel

# Extract SNPs
gatk SelectVariants \
  -R ${ref} \
  -V ${results}/raw_variants.vcf \
  --select-type SNP \
  -O ${results}/raw_snps.vcf

# Extract INDELs
gatk SelectVariants \
  -R ${ref} \
  -V ${results}/raw_variants.vcf \
  --select-type INDEL \
  -O ${results}/raw_indels.vcf


echo "SNP and Indel extracted"

#==================== MultiQC to summerize

cd /Users/azadehnikouee/Desktop/GATK_germline_variant_calling/

multiqc . -o multiqc_report








