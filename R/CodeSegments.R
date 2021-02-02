# -------------------------- #
# Code segments from subclones files
# -------------------------- #

#' @name
#' CodeSegments
#'
#' @title
#' Code all segments in all samples in the cohort
#'
#' @description
#' Codes all segments as a particular CNA type from across the cohort
#'
#' @param
#' segfile_dir = where segfile will go
#' segfile_name = cohort name
#'
#' @return
#' A table for use with bedtools etc
#'
#'



# code all segments into homdel, loh, otherloss, no change, gain, big gain
CodeSegments <- function(segfile_dir,segfile_name) {

  subs <- read.table(paste0(segfile_dir,
                            segfile_name,
                            "_segsfull.txt"),
                     sep="\t",
                     hea=T,
                     stringsAsFactors=F)

  subs$total_cn = rowSums(subs[c("nMajor","nMinor")])

  # code dip,tetra to 2,4
  subs$class[subs$class=="dip"] <- 2
  subs$class[subs$class=="tetra"] <- 4
  colnames(subs)[14]="dip.tetra"
  subs$dip.tetra = as.numeric(subs$dip.tetra)

  # coding each CNA type
  subs$coded_total_cn = NA
  # homdels
  subs$coded_total_cn[which(subs$total_cn==0)] = "homdel"
  # LOH
  subs$coded_total_cn[which(subs$nMajor>0 & subs$nMinor==0)] = "loh"
  # other loss (non homdel, non LOG)
  subs$coded_total_cn[which(subs$total_cn < subs$dip.tetra
                            & subs$total_cn != 0
                            & subs$nMinor != 0)] = "otherloss"
  # no change from ploidy class
  subs$coded_total_cn[which((subs$nMajor==subs$nMinor | subs$nMajor==3 & subs$nMinor==1)
                            & (subs$total_cn == subs$dip.tetra))] = "nochange"
  # gain
  subs$coded_total_cn[which(subs$total_cn>subs$dip.tetra)] = "gain"
  # big gain
  subs$coded_total_cn[which(subs$total_cn>(5*subs$dip.tetra))] = "biggain"

  # change X to 23
  subs$chr[subs$chr=="X"] = 23

  # separate into CNA types
  subs <- subs[order(subs$chr,subs$startpos,subs$endpos),]
  subs$chr = paste0("chr",subs$chr)
  cnv_levels <- c("homdel","loh","otherloss","nochange","gain","biggain")
  cnv_numbers = c(-2,-1,-0.5,2,5,10)
  for (cnv in 1:length(cnv_levels)) {

    cnv_subs = subs[which(subs$coded_total_cn==cnv_levels[cnv]),]
    cnv_subs = cnv_subs[order(cnv_subs$chr,cnv_subs$startpos,cnv_subs$endpos),]

    # write full file
    write.table(cnv_subs,paste0(segfile_dir,segfile_name,"_",cnv_levels[cnv],".full"),
                quote=F,
                row.names=F,
                sep="\t")

    # write as bedfile for bedtools
    subs_bed = cbind(cnv_subs[c("chr","startpos","endpos","sample")])
    colnames(subs_bed)[4]="name"
    print(nrow(subs_bed))
    write.table(subs_bed,paste0(segfile_dir,segfile_name,"_",cnv_levels[cnv],".bed"),
                quote=F,
                row.names=F,
                col.names=F,
                sep="\t")
  }




}
