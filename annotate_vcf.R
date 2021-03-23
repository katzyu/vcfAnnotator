# VCF Variant Annotation Tool
#
# Taking a VCF file as input, each variant will be annotated with the following
# pieces of information:
#  1. Type of variation (substitution, insertion, CNV, etc.) and their effect (missense,
#  silent, intergenic, etc.). If there are multiple effects, annotate with the
#  most deleterious possibility.
#  2. Depth of sequence coverage at the site of variation.
#  3. Number of reads supporting the variant.
#  4. Percentage of reads supporting the variant versus those supporting reference reads.
#  5. Allele frequency of variant from ExAC API

#-----------------------------------------------------------------------
#                            Load packages
#-----------------------------------------------------------------------

library(httr)
library(jsonlite)

#-----------------------------------------------------------------------
#                        Function definitions
#-----------------------------------------------------------------------

check_vcf_file <- function(input_file){
  if (!file.exists(input_file) | dir.exists(input_file)) {
    stop("File not found. Please enter a VCF file.")
  }
  if (tools::file_ext(input_file) != "vcf") {
    stop("File is not a VCF file. Please check file.")
  }
  vcf_format_line <- readLines(input_file, n = 1)
  if (substr(vcf_format_line, 1, 16) != "##fileformat=VCF") {
    stop("File is not a properly formatted VCF file. Please check file.")
  }
}

check_output_dir <- function(output_file){
  if (!dir.exists(dirname(output_file))) {
    stop("Output directory does not exist.")
  }
}

extract_colnames <- function(vcf_file) {
  # search for line containing #CHROM which indicates a line with
  # the variant column names
  con <- file(vcf_file, "r")
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if (grepl("#CHROM", line)) {
      break
    }
  }
  close(con)
  
  # remove # sign from column names and split string
  line <- gsub("#", "", line)
  line <- strsplit(line, "\t")[[1]]
  return(line)
}

read_vcf <- function(vcf_file) {
  # read in vcf variant data
  data <- read.table(vcf_file, sep = "\t", comment.char = "#", stringsAsFactors = F)
  
  # find column names from vcf file and set colnames
  vcf_colnames <- extract_colnames(vcf_file)
  colnames(data) <- vcf_colnames
  return(data)
}

parse_vcf_info <- function(info_entry) {
  # Parses the entry in the INFO column into a named vector
  info_data <- gsub(".*=", "", strsplit(info_entry, ";")[[1]])
  names(info_data) <- gsub("=.*", "", strsplit(info_entry, ";")[[1]])
  return(info_data)
}

generate_vep_scores <- function(ensembl_vep_table) {
  # Creates a named vector of ensembl VEP scores from VEP table
  # names = effect (e.g. missense_variant), values = score
  # VEP scores reflect the impact of effect (high score = more deleterious)
  # Table was copied and saved as a tsv from:
  # http://uswest.ensembl.org/info/genome/variation/prediction/predicted_data.html
  ensembl_vep_table <- read.table(ensembl_vep_table, sep = "\t", header = T)
  ensembl_vep_table$score <- seq(nrow(ensembl_vep_table), 1, -1)
  ensembl_vep_scores <- c(ensembl_vep_table$score)
  names(ensembl_vep_scores) <- ensembl_vep_table$SO.term
  return(ensembl_vep_scores)
}

query_exac_api_bulk <- function(exac_query) {
  # Submits bulk query to ExAC REST API
  query_results <- POST("http://exac.hms.harvard.edu/rest/bulk/variant/variant", body = toJSON(exac_query))
  
  # check response code status
  if (query_results$status_code != 200){
    stop("ExAC query failed. Please check internet connection or ExAC server.")
  }
  query_results <- fromJSON(rawToChar(query_results$content))
  return(query_results)
}

get_exac_query_name <- function(exac_annotations_raw_entry){
  exac_query <- names(exac_annotations_raw_entry)
  return(exac_query)
}

get_vep_consequence <- function(exac_annotations_raw_entry, ensembl_vep_scores){
  vep_consequence <- exac_annotations_raw_entry[[1]]$vep_annotations$Consequence[1]
  if (is.null(vep_consequence)) {
    # If no consequence is found, replace NULL with NA
    vep_consequence <- NA
  } else {
    # If consequence is found, return most deleterious consequence based on VEP score
    vep_consequence_list <- strsplit(vep_consequence, "&")[[1]]
    vep_consequence <- names(ensembl_vep_scores)[ensembl_vep_scores == max(ensembl_vep_scores[vep_consequence_list])]
  }
  return(vep_consequence)
}

get_allele_freq <- function(exac_annotations_raw_entry){
  allele_frequency <- exac_annotations_raw_entry[[1]]$allele_freq
  if (is.null(allele_frequency)) {
    allele_frequency <- NA
  }
  return(allele_frequency)
}

generate_exac_annotations <- function(exac_queries, ensembl_vep_scores) {
  # Fetches VEP consequences and allele frequencies from ExAC REST API
  # Returns data frame with information from ExAC
  # If multiple VEP consequences are present, chooses the most deleterious

  print("Querying ExAC API")
  exac_annotations_raw <- query_exac_api_bulk(exac_queries)

  # Extract VEP consequences and allele frequency info from ExAC query results
  print("Extracting features of interest from ExAC query results")
  exac_annotations <- data.frame()
  pb <- txtProgressBar(min = 0, max = length(exac_annotations_raw), style = 3)
  for (i in 1:length(exac_annotations_raw)) {
    exac_query <- get_exac_query_name(exac_annotations_raw[i])
    vep_consequence <- get_vep_consequence(exac_annotations_raw[i], ensembl_vep_scores)
    allele_frequency <- get_allele_freq(exac_annotations_raw[i])
    exac_annotations <- rbind(exac_annotations, cbind(exac_query, vep_consequence, allele_frequency))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  # ExAC REST API returns queries out of order; reorder rows to match original queries
  exac_annotations <- exac_annotations[match(exac_queries, exac_annotations$exac_query), ]
  exac_annotations$exac_query <- NULL
  return(exac_annotations)
}

get_variant_type <- function(vcf_info){
  variant_type <- as.character(vcf_info["TYPE"])
  return(variant_type)
}

get_total_read_depth <- function(vcf_info){
  # get total read depth at locus
  total_read_depth <- as.numeric(vcf_info["DP"])
  return(total_read_depth)
}

get_alt_depth <- function(vcf_info){
  # get number of alt reads 
  alt_depth <- as.numeric(strsplit(vcf_info["AO"], ",")[[1]])
  if (length(alt_depth) > 1) {
    # If there are multiple alt variants, sum the alt counts
    alt_depth <- sum(alt_depth)
  }
  return(alt_depth)
}

get_alt_pct <- function(alt_depth, total_read_depth){
  # get percent of alt reads/total reads
  pct_alt <- alt_depth / total_read_depth * 100
  return(pct_alt)
}

annotate_vcf <- function(vcf_file, num_lines = -1) {
  # Main function to create the annotated vcf data frame
  vcf <- read_vcf(vcf_file)

  # Create annoted vcf data frame by extracting features of interest from the vcf file
  print("Extracting features of interest from vcf file:")
  annotated_vcf <- data.frame()
  pb <- txtProgressBar(min = 0, max = nrow(vcf), style = 3)
  for (i in 1:nrow(vcf)) {
    chrom <- vcf$CHROM[i]
    pos <- vcf$POS[i]
    ref <- vcf$REF[i]
    alt <- vcf$ALT[i]
    
    # Get read depth, alt depth, and alt percent from INFO column
    vcf_info <- parse_vcf_info(vcf$INFO[i])
    variant_type <- get_variant_type(vcf_info)
    total_read_depth <- get_total_read_depth(vcf_info)
    alt_depth <- get_alt_depth(vcf_info)
    pct_alt <- get_alt_pct(alt_depth, total_read_depth)
    
    annotated_vcf <- rbind(annotated_vcf, cbind(chrom, pos, ref, alt, variant_type, total_read_depth, alt_depth, pct_alt))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Fetch VEP consequence and allele frequency from ExAC REST API
  exac_queries <- paste(vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT, sep = "-")
  ensembl_vep_scores <- generate_vep_scores("data/ensembl_vep_table.txt")
  exac_annotations <- generate_exac_annotations(exac_queries, ensembl_vep_scores)

  # Combine vcf annotations and ExAC annotations to make final data frame
  final_annotated_vcf <- cbind(annotated_vcf, exac_annotations)
  return(final_annotated_vcf)
}


#-----------------------------------------------------------------------
#                         Annotate VCF file
#-----------------------------------------------------------------------

# Get input file and output file name
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Check user inputs
if (length(args) != 2) {
  msg <- paste("Please supply a vcf input file and output path.\n")
  msg <- paste(msg, "Command should like like:\n")
  msg <- paste(msg, "Rscript annotate_vcf.R <input_vcf> <output_file>")
  stop(msg)
}

# Check vcf file
check_vcf_file(input_file)

# Check existance of output directory
check_output_dir(output_file)

# Annotate variants
final_annotated_vcf <- annotate_vcf(input_file)
write.table(final_annotated_vcf, file = output_file, quote = FALSE, row.names = FALSE, sep = "\t")
print("Done annotating files.")