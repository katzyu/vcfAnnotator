# VCF annotation tool

This repository contains a VCF annotation tool written for the Tempus Bioinformatics Technical Challenge.  

## Requirements
R version 3.6.2 or later and the following R packages:
* httr (>= 1.4.2)
* jsonlite (>= 1.7.2)

## Installation

Clone repo to local computer
```bash
git clone https://github.com/katzyu/vcfAnnotator.git
```
Install the package dependencies in R:

```bash
install.packages(c("httr", "jsonlite"))
```

## Usage

The variant annotation tool takes in a VCF file and extracts the following annotated features from the original VCF file and queries to the Broad Institute's [ExAC browser](http://exac.hms.harvard.edu/):
1. Chromosome
2. Position
3. Reference allele(s)
4. Alternate (non-reference) allele(s)
5. Type of variant (e.g. SNP, MNP, deletion, insertion)
6. Total read depth at locus
7. Number of reads supporting the variant.
8. Percent of reads supporting variant
9. Effect of variant (e.g. intergenic_variant, transcript_ablation)
10. Allele frequency of variant from ExAC API

When choosing the most deleterious effect, [Ensembl's Variant Effect Predictor (VEP)](https://uswest.ensembl.org/info/docs/tools/vep/index.html) terms were used to rank the variant effects by their potential impact. The table of ranked VEP terms was copied from [Ensembl's Calculated variant consequences page]( http://uswest.ensembl.org/info/genome/variation/prediction/predicted_data.html) and saved as a text file in the data directory (ensembl_vep_table.txt). 

### Annotate VCF file

To run the variant annotation tool, a vcf file and output path must be provided by the user. An example command is provided:
```bash
cd vcfAnnotator
Rscript annotate_vcf.R 'data/Challenge_data_(1).vcf' data/Challenge_data_annotated.txt
```
