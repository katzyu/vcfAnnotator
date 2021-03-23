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

read_vcf <- function(file, num_lines){
  # Reads in vcf file and returns data frame of variants
  if(!file.exists(file) | dir.exists(file)){
    stop("File not found. Please enter a VCF file.")
  }
  vcf_format_line <- readLines(file, n = 1)
  if(tools::file_ext(file) != "vcf"){
    stop("File is not a VCF file. Please check file.")
  }
  if(substr(vcf_format_line, 1, 16) != "##fileformat=VCF"){
    stop("File is not a properly formatted VCF file. Please check file.")
  }
  data <- read.table(file, sep = "\t", comment.char = "#", stringsAsFactors = F)
  
  # Search for column names in vcf file. 
  # If vcf file is big, can read just the first <num_lines> 
  vcf_colnames <- strsplit(gsub("#", "", 
                                grep("#CHROM",readLines(file, n = num_lines), value = TRUE)), "\t")[[1]]
  colnames(data) <- vcf_colnames
  return(data)
}

parse_vcf_info<- function(info_entry){
  # Parses the INFO column in vcf file into a named vector
  info_data <- gsub(".*=", "", strsplit(info_entry, ";")[[1]])
  names(info_data) <- gsub("=.*", "", strsplit(info_entry, ";")[[1]])
  return(info_data)
}

generate_vep_scores <- function(ensembl_vep_table){
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

query_exac_api_bulk <- function(exac_query){
  # Submits bulk query to ExAC REST API 
  query_results <- POST("http://exac.hms.harvard.edu/rest/bulk/variant/variant", body = toJSON(exac_query))
  query_results <- fromJSON(rawToChar(query_results$content))
  return(query_results)
}

generate_exac_annotations <- function(exac_queries, ensembl_vep_scores){
  # Fetches VEP consequences and allele frequencies from ExAC REST API
  # Returns data frame with information from ExAC
  # If multiple VEP consequences are present, chooses the most deleterious
  
  print("Querying ExAC API")
  exac_annotations_raw <- query_exac_api_bulk(exac_queries)
  
  # Extract VEP consequences and allele frequency info from ExAC query
  print("Extracting features of interest from ExAC query")
  exac_annotations <- data.frame()
  pb <- txtProgressBar(min = 0, max = length(exac_annotations_raw), style = 3)
  for (i in 1:length(exac_annotations_raw)){
    exac_query <- names(exac_annotations_raw)[i]
    vep_consequence <- exac_annotations_raw[[i]]$vep_annotations$Consequence[1]
    allele_frequency <- exac_annotations_raw[[i]]$allele_freq
    
    if (is.null(vep_consequence)){
      # If no consequence is found, replace NULL with NA
      vep_consequence <- NA
    } else {
      # If consequence is found, return most deleterious consequence based on VEP score
      vep_consequence_list <- strsplit(vep_consequence, "&")[[1]]
      vep_consequence = names(ensembl_vep_scores)[ensembl_vep_scores == max(ensembl_vep_scores[vep_consequence_list])]
    }
    
    if (is.null(allele_frequency)){
      allele_frequency <- NA
    }
    exac_annotations <- rbind(exac_annotations, cbind(exac_query, vep_consequence, allele_frequency))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # ExAC REST API returns queries out of order; reorder rows to match original queries
  exac_annotations <- exac_annotations[match(exac_queries, exac_annotations$exac_query), ]
  exac_annotations$exac_query <- NULL
  return(exac_annotations)
}

annotate_vcf <- function(vcf_file, num_lines = -1){
  # Main function to create the annotated vcf data frame
  
  vcf <- read_vcf(vcf_file, num_lines)
  ensembl_vep_scores <- generate_vep_scores("../data/ensembl_vep_table.txt")
  
  # Create annoted vcf data frame by extracting features of interest from the vcf file
  print("Extracting features of interest from vcf file:")
  annotated_vcf <- data.frame()
  pb <- txtProgressBar(min = 0, max = nrow(vcf), style = 3)
  for (i in 1:nrow(vcf)){
    chrom <- vcf$CHROM[i]
    pos <- vcf$POS[i]
    ref <- vcf$REF[i]
    alt <- vcf$ALT[i]
    vcf_info <- parse_vcf_info(vcf$INFO[i])
    total_read_depth <- as.numeric(vcf_info["DP"])
    alt_depth <- as.numeric(strsplit(vcf_info["AO"], ",")[[1]])
    if(length(alt_depth)>1){
      alt_depth <- sum(alt_depth)
    }
    pct_alt <- alt_depth/total_read_depth * 100
    variant_type <- as.character(vcf_info["TYPE"])
    annotated_vcf <- rbind(annotated_vcf, cbind(chrom, pos, ref, alt, variant_type, total_read_depth, alt_depth, pct_alt))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Fetch VEP consequence and allele frequency from ExAC REST API
  exac_queries <- paste(vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT, sep = "-")
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
if (length(args) != 2){
  msg <- paste("Please supply a vcf input file and output path.\n")
  msg <- paste(msg, "Command should like like:\n")
  msg <- paste(msg, "Rscript annotate_vcf.R <input_vcf> <output_file>")
  stop(msg)
}

# Annotate variants
final_annotated_vcf <- annotate_vcf(input_file)
write.table(final_annotated_vcf, file = output_file, quote = FALSE, row.names = FALSE, sep = "\t")
print("Done annotating files.")
