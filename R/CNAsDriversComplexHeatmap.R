# -------------------------- #
# Make complexheatmap
# -------------------------- #

#' @name
#' CNAsDriversComplexHeatmap
#'
#' @title
#' Heatmap to show CNA status in a set of driver genes across a cohort of tumour samples
#' Dendrogram automatically goes alongside
#' Bars classifying samples in lots of ways produced alongside heatmap
#'
#' @description
#' Create heatmap showing CNA status at driver genes using complexheatmap
#' homdel, loh, otherloss, nochange, gain, biggain
#'
#' @param
#' segfile_dir = dir with all *_subclones.txt files are
#' segfile_name = label of cohort, eg 'TGCT'
#' bed_file = tab sep file, 7 cols, chr,start,end,strand,ens,ens,genename
#' driver_file = tab sep file, colnames SYMBOL (of gene) and TIER (1,2,etc) existing
#'
#' @return
#' Heatmap, dendrogram of CNAs v drivers, with classifying bars
#'

# collate all segments in subclones files across cohort
CNAsDriversComplexHeatmap <- function(segfile_dir,
                                         segfile_name,
                                         bed_file,
                                         driver_file) {

  # label drivers as Tier 1 or 2
  # get drivers and chr positions
  bed = read.table(bed_file,
                   stringsAsFactors = F)
  dri=read.table(driver_file,
                 hea=T,
                 stringsAsFactors = F,
                 fill=T)
  # select drivers in tier 1 (highest confidence) and 2 (rescue)
  dri = dri[which(dri$TIER==1 | dri$TIER==2),]
  bed=bed[match(dri$SYMBOL,bed[,7]),]
  positions_chr = as.character(bed[,1])
  positions_chr = as.character(gsub("chr","",positions_chr))
  positions_start = as.character(bed[,2])
  positions_end = as.character(bed[,3])
  positions_gene = as.character(bed[,7])
  genes = as.data.frame(cbind(positions_gene,positions_chr,positions_start,positions_end))
  for (i in 1:ncol(genes)){genes[,i]=as.character(genes[,i])}
  genes$positions_chr[which(genes$positions_chr=="X")]="23"

  # load tomelt, used to make previous heatmap
  tomelt <- as.matrix(read.csv(paste0(segfile_dir,segfile_name,"_dendrogram.csv"),
                     stringsAsFactors = F,
                     sep=" "))

  # remove colnames
  samples = rownames(tomelt)
  tumour_platekeys = sapply(samples,function(x){strsplit(samples,"\\.")[[1]][2]})
  rownames(tomelt)=NULL


  # get MSI/MSS status
  # get MSI classification
  msi=read.table("/re_gecip/cancer_ovarian/endometrial_current/analyses_v8/MSINGS/output/combined_results.txt",
                 stringsAsFactors = F,
                 hea=T)
  msistatus = msi$msi_status[match(tumour_platekeys,msi$Position)]

  pdf("")
  Heatmap(tomelt, name = "ENDOMETRIAL", row_km = 3, column_km=3,
          top_annotation = HeatmapAnnotation(foo1 = 1:48),
          right_annotation = rowAnnotation(msi=msistatus))
  dev.off()

  # get purity



  # get ploidy



  # Duke's stage



  # location of tumour



  # get TMB (tumour mutation burder)



  # get number SVs (structural variants)




}
