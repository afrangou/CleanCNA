% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CollateSubclones.R
\name{CollateSubclones}
\alias{CollateSubclones}
\title{Collate all segments in all samples in the cohort}
\usage{
CollateSubclones(qc, segfile_dir, segfile_name)
}
\arguments{
\item{qc}{table from PrepDataForSummary
segfile_dir = where segfile will go
segfile_name = cohort name}
}
\value{
A table with all segments from all subclones files collated
}
\description{
Prepares subclones info to label copy number of segments
}
