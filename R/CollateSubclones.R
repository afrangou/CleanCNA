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
#' A table for use with bedtools etc
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
  cnv_levels <- c("homdel","lohs","nochange","gain","biggain")
  cnv_numbers = c(-2,-1,-0.5,2,5,10)
  for (cnv in 1:length(cnvs_levels)) {

    cnv_subs = subs[which(subs$coded_total_cn==cnv_levels[cnv]),]
    cnv_subs = cnv_subs[order(cnv_subs$chr,cnv_subs$startpos,cnv_subs$endpos),]

    # write full file
    write.table(cnv_subs,paste0(segfile_dir,cnv_levels[cnv],"_",segfile_name,"_normalised.full"),
                quote=F,
                row.names=F,
                sep="\t")

    # write as bedfile for bedtools
    subs_bed = cbind(cnv_subs[c("chr","startpos","endpos","sample")])
    colnames(subs_bed)[4]="name"
    print(nrow(subs_bed))
    write.table(subs_bed,paste0(segfile_dir,cnv_levels[cnv],"_",segfile_name,"_normalised.bed"),
                quote=F,
                row.names=F,
                col.names=F,
                sep="\t")
  }




}










