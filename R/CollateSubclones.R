# -------------------------- #
# Collate all segments from subclones files
# -------------------------- #

#' @name
#' CollateSubclones
#'
#' @title
#' Collate all segments in all samples in the cohort
#'
#' @description
#' Prepares subclones info to label copy number of segments
#'
#' @param
#' qc table from PrepDataForSummary
#' segfile_dir = where segfile will go
#' segfile_name = cohort name
#'
#' @return
#' A table with all segments from all subclones files collated
#'
#'

# collate all segments in subclones files across cohort
CollateSubclones <- function(qc,
                             segfile_dir,
                             segfile_name) {

  # structure of final table
  subsfull = matrix(nrow=0,ncol=13)
  #colnames(subsfull) = c("sample","chr","startpos","endpos","nMajor","nMinor")

  # paths of subclones files
  qc <- read.csv(qc,
                 sep="\t",
                 stringsAsFactors=F)

  samples <- paste0(qc$participant_id,".",
                    qc$tumour_sample_platekey,".",
                    qc$germline_sample_platekey)

  # get subclones paths
  subclones <- list.files(segfile_dir)
  subclonesPaths <- paste0(segfile_dir,subclones[grep("subclones",subclones)])

  # get purity paths
  purity <- list.files(segfile_dir)
  purityPaths <- paste0(segfile_dir,purity[grep("cellularity",purity)])

  # check same number and name of purity and subclones paths

  for (i in 1:length(subclonesPaths)) {

    sub <- read.csv(subclonesPaths[i],
                   hea=T,
                   stringsAsFactors=F,
                   sep="\t")[,c(1:13)]

    samplename <- paste0(qc$participant_id[i],".",
                         qc$tumour_sample_platekey[i],".",
                         qc$germline_sample_platekey[i])

    # switch subclones so biggest subclone first
    sub$subclonal = sub$frac1_A > sub$frac2_A
    firstsubclone = sub[,8:10]
    secondsubclone = sub[,11:13]
    sub[which(sub$subclonal==F),8:10] = secondsubclone[which(sub$subclonal==F),]
    sub[which(sub$subclonal==F),11:13] = firstsubclone[which(sub$subclonal==F),]
    sub=sub[,-c(4:7,14)]

    sub = cbind(as.character(samples[i]),sub)

    colnames(sub)[1:6]=c("sample","chr","startpos","endpos","nMajor","nMinor")

    # add ploidy for sample
    sub$ploidy = signif(read.table(purityPaths[i],hea=T,stringsAsFactors=F)[1,3])

    # length of segment
    sub$size = sub$endpos - sub$startpos +1

    # categorise as diploid or tetraploid (PCAWG eqn)
    # minor CN of 0 after rounding CN
    # use separate table for this
    sub2 <- sub
    sub2 <- sub2[1:max(which(sub2$chr==22)),]
    sub2[is.na(sub2)] <- 0
    # LOH regions
    genome.length <- 2923364639
    sub2fracLOH <- sub2[which(round((sub2$nMinor*sub2$frac1_A)+(sub2$nMin2_A*sub2$frac2_A))==0),]
    sub2fracLOH <- sum(sub2fracLOH$endpos-sub2fracLOH$startpos)/genome.length
    # calculate psi_t, different to average tumour ploidy
    sub$psi_t <- sum(sub2$size * (sub2$nMajor + sub2$nMinor))/sum(sub2$size)
    psi_t <- sub$psi_t
    # classify sample
    sampleploidy <- ifelse(psi_t > (sub2fracLOH*-2)+2.9,
                      "tetra",
                      "dip")
    sub$class <- sampleploidy

    # add to subsfull
    subsfull <- rbind(subsfull,sub)

    # write.table(sub,paste0(segfile_dir,segfile_name,"_",
    #                        qc$tumour_sample_platekey[i],
    #                        "_segsfull.txt"),
    #             sep="\t",
    #             row.names=F,
    #             quote=F)
  }

  # save full table
  write.table(subsfull,paste0(segfile_dir,
                            segfile_name,
                            "_segsfull.txt"),
            sep="\t",
            row.names=F,
            quote=F)

}













