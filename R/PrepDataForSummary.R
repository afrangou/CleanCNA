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

prepDataForSummary <- function(qctabledir, subclonesdir) {

  qc <- read.csv(qctabledir,
                 sep="\t",
                 stringsAsFactors=F)

  outputTable <- qc[,c("participant_id",
                       "tumour_sample_platekey",
                       "germline_sample_platekey")]

  subclonesPaths <- paste0(subclonesdir,
                           outputTable$tumour_sample_platekey,
                           "_subclones.txt")

  purityPaths <- paste0(subclonesdir,
                        outputTable$tumour_sample_platekey,
                        "_cellularity_ploidy.txt")

  qc <- cbind(qc,subclonesPaths,purityPaths)

}

