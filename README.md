# VCF annotation tool

This repository contains a VCF annotation tool written for the Tempus Bioinformatics Technical Challenge.  

## Requirements
R version 3.6.2 or later and the following R packages:
* httr (>= 1.4.2)
* jsonlite (>= 1.7.2)

## Installation

Clone repo to local computer
```bash
git clone https://github.com/katharineyu/VCFannotator.git
```
Install the package dependencies in R:

```bash
install.packages(c("httr", "jsonlite"))
```

## Usage

The variant annotation tool takes in a VCF file and extracts the following annotated features from the original VCF file and queries to the Broad Institute's [ExAC browser](http://exac.hms.harvard.edu/):
1. Type of variation (substitution, insertion, CNV, etc.)
2. Effect of variant (missense, silent,
intergenic, etc.). If there are multiple effects, the most deleterious possibility is chosen (see below).
3. Depth of sequence coverage at the site of variation.
4. Number of reads supporting the variant.
5. Percentage of reads supporting the variant versus those supporting reference reads.
6. Allele frequency of variant from ExAC API

When choosing the most deleterious effect, Ensembl's Variant Effect Predictor (VEP) terms were used to rank the variant effects by their potential impact. The table of ranked VEP terms was copied from [Ensembl's Calculated variant consequences page]( http://uswest.ensembl.org/info/genome/variation/prediction/predicted_data.html) and saved as a text file in the data directory (ensembl_vep_table.txt). 

### Annotate VCF file

To run the variant annotation tool, a vcf file and output path should be provided by the user. An example command is provided:
```bash
cd VCFannotator
Rscript annotate_vcf.R '../data/Challenge_data_(1).vcf' ../data/Challenge_data_annotated.txt
```