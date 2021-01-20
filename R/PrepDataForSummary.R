# -------------------------- #
# Make table with
# Part ID
# Tumour platekey
# Germline platekey
# Filepath for subclones file for passed run
# Filepath for cellularity_ploidy file for passed run
# -------------------------- #

#' @name
#' PrepDataForSummary
#'
#' @title
#' Make summary table with filepaths for required files (subclones and purity)
#'
#' @description
#' Converts snv counts from vcf into input format for DPClust
#'
#' @param
#' Filepath of QualityControl.tsv file
#' Dir of subclones and cellularity files
#'
#' @return
#' A summary table

x=3
