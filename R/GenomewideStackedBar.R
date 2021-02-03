# -------------------------- #
# Prepare data for stacked barplot
# -------------------------- #

#' @name
#' GenomewideStackedBar
#'
#' @title
#' Prep data for genome wide stacked barplot
#' Create stacked bar plot
#'
#' @description
#' Take output from bedtools intersect and prepare data for stacked barplot
#' Create stacked barplot for 6 different CNA types
#' homdel, loh, otherloss, nochange, gain, biggain
#'
#' @param
#' filestub = dir (with trailling slash) containing *_partitioned_regions_all_CNA_types_overlapped.out files for all cna types
#' segfile_name = label of cohort, eg 'TGCT'
#'
#' @return
#' Prepped data
#' Stacked barplot
#'

# collate all segments in subclones files across cohort
GenomewideStackedBar <- function(filestub,segfile_name) {

  all = read.table(paste0(filestub,segfile_name,
                   "_partitioned_regions_all_CNA_types_overlapped.out"),
                   sep="\t",
                   stringsAsFactors=F)

  pos = paste(all[,1],all[,2],all[,3],sep="_")
  all = cbind(all,pos)
  pos = unique(pos)
  poschr = sapply(pos,function(x){strsplit(x,"_")[[1]][1]})
  posleft = sapply(pos,function(x){strsplit(x,"_")[[1]][2]})
  posright = sapply(pos,function(x){strsplit(x,"_")[[1]][3]})

  forplot = cbind(poschr,posleft,posright)
  homdel = rep(0,nrow(forplot))
  loh = rep(0,nrow(forplot))
  otherloss = rep(0,nrow(forplot))
  nochange = rep(0,nrow(forplot))
  gain = rep(0,nrow(forplot))
  biggain = rep(0,nrow(forplot))

  forplot = as.data.frame(cbind(forplot,homdel,loh,otherloss,nochange,gain,biggain))
  for (i in 2:ncol(forplot)) {forplot[,i]=as.integer(as.character(forplot[,i]))}
  cnas = c("homdel","loh","otherloss","nochange","gain","biggain")

  for (cna in 1:length(cnas)) {

    # reduce table to cnas of a specific type and only the genomic regions where they exist
    allcna = all[which(all[,9]==cnas[cna] & all[,8]!=0),]

    # find these regions in the forplot table and fill in for current cna type
    test = match(rownames(forplot),allcna$pos)
    forplot_rows = which(!is.na(test))
    allcna_rows = test[-which(is.na(test))]
    forplot[forplot_rows,cnas[cna]] = allcna[cna_rows,8]

    #     for (i in 1:nrow(allcna)) {
    #   row = which(rownames(forplot)==allcna$pos[i])
    #   if (length(row)>0) {
    #     forplot[row,cnas[cna]] = allcna[i,8]
    #   }
    #
    # }

  }

  write.table(forplot,paste0(filestub,"_",segfile_name,
                             "_partitioned_regions_all_CNA_types_overlapped_forplot.out"),
              quote=F,
              sep="\t")
}
