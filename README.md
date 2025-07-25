GATK Germline Variant Calling Pipeline – Training Project (HG00096)
This document provides a complete guide and description of the shell script used for performing GATK Best Practices Germline Variant Calling using the sample SRR070802 from the 1000 Genomes Project. The pipeline goes from raw sequencing reads to SNP and INDEL variant extraction.
-Sample Information
Sample ID: SRR070802
Source: https://www.ebi.ac.uk/ena/browser/view/SRR070802
Reference Genome: Human genome hg38
GATK Resource Bundle: https://gatk.broadinstitute.org
-Requirements
Make sure the following tools are installed and available in your PATH:
- wget
- samtools
- bwa
- gatk
- fastqc
- multiqc
- Pipeline Overview
1. Download paired-end FASTQ files (SRR070802) from ENA.
2. Download and index hg38 reference genome.
3. Download dbSNP known sites VCF for BQSR.
4. Run FastQC on FASTQ reads.
5. Align reads with BWA-MEM.
6. Mark duplicates and sort BAM using GATK.
7. Perform Base Quality Score Recalibration (BQSR).
8. Collect alignment and insert size metrics.
9. Call variants using GATK HaplotypeCaller.
10. Extract SNPs and INDELs from raw VCF.
11. Generate a final summary report using MultiQC.
-Directory Structure

GATK_germline_variant_calling/

├── reads/                       # FASTQ input files

├── reference/                   # hg38 reference and known sites

├── aligned_read/               # Aligned BAM files

├── logs/                       # Log and metrics files

├── result/                     # VCF outputs and recalibration table

├── FASTQC/                     # FastQC reports

└── multiqc_report/             # MultiQC HTML summary

- Notes

This script is designed for training and educational purposes using a single human sample. Update the file paths to match your environment as needed. Adjust thread (-t) settings based on available resources.
## Citation
If you use this pipeline, please cite:
- GATK Best Practices: https://gatk.broadinstitute.org/
- 1000 Genomes Project: https://www.internationalgenome.org/
- ENA Sample SRR070802: https://www.ebi.ac.uk/ena/browser/view/SRR070802

##===============================
Additional Analysis: Ti/Tv Ratio and Substitution Patterns
This step uses `bcftools` and `plot-vcfstats` to evaluate the quality of variant calls by calculating the Transition/Transversion (Ti/Tv) ratio and visualizing base substitution patterns.

## Requirements:
- `bcftools`
- `plot-vcfstats` (included with bcftools)
- (Optional) `tectonic` or `pdflatex` for generating PDF reports from LaTeX

## What It Does:
1. Runs `bcftools stats` on the filtered SNP VCF (`raw_snps.vcf`)
2. Generates plots such as:
   - Ti/Tv ratio vs. quality
   - Base substitution types
   - Depth distribution
3. (Optional) Converts LaTeX summary (`summary.tex`) into a PDF if `tectonic` is installed

## Output:
- “result/vcf_stats.txt”: bcftools stats output
- “result/vcf_plots/*.png”: Visual summary of variant quality
- “result/vcf_plots/summary.tex”: LaTeX version of report
- “result/vcf_plots/summary.pdf”: Optional PDF (if LaTeX is installed)

This analysis provides a useful QC check. A Ti/Tv ratio around 2.0–3.0 indicates good-quality SNP calling for human data.  





