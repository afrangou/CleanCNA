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
#' with various bars down the side
#' for bars down side, feed vectors to function
#'
#' @param
#' segfile_dir = dir with all *_subclones.txt files are
#' segfile_name = label of cohort, eg 'TGCT'
#' plot_name = name for plot
#' tomelt = saved from CNAsDriversHeatmap function
#' bed_file =
#' driver_file = tab sep file, colnames SYMBOL (of gene) and TIER (1,2,etc) existing
#'
#' @return
#' Heatmap, dendrogram of CNAs v drivers, with classifying bars
#'

CNAsDriversComplexHeatmap <- function(segfile_dir,
                             segfile_name,
                             plot_name,
                             bed_file,
                             driver_file,
                             tomelt) {

  library(ComplexHeatmap)

  bed = read.table(bed_file, stringsAsFactors = F)
  dri = read.table(driver_file, hea = T, stringsAsFactors = F,
                   fill = T)

  dri = dri[which(dri$TIER == 1 | dri$TIER == 2), ]
  bed = bed[match(dri$SYMBOL, bed[, 7]), ]
  bed = cbind(bed,dri$TIER)
  positions_chr = as.character(bed[, 1])
  positions_chr = as.character(gsub("chr", "", positions_chr))
  positions_start = as.character(bed[, 2])
  positions_end = as.character(bed[, 3])
  positions_gene = as.character(bed[, 7])
  genes = as.data.frame(cbind(positions_gene, positions_chr,
                              positions_start, positions_end))
  genes$tier =
    for (i in 1:ncol(genes)) {
      genes[, i] = as.character(genes[, i])
    }
  genes$positions_chr[which(genes$positions_chr == "X")] = "23"
  tomelt <- as.matrix(read.csv(paste0(segfile_dir, segfile_name,
                                      "_dendrogram.csv"), stringsAsFactors = F, sep = " "))
  samples = rownames(tomelt)
  rownames(tomelt)=NULL
  tumour_platekeys = sapply(samples, function(x) {
    strsplit(samples, "\\.")[[1]][2]
  })


  Cairo:::CairoPDF(paste0(segfile_dir,segfile_name,"_",plot_name),width = 10,height=6)
  Heatmap(tomelt,name="CN",
          use_raster=T,raster_device="CairoPNG",
          rect_gp = gpar(col = "white", lwd = 1),
          # row_split = subtype,
          # row_split = semornon,
          # clustering_method_rows = "single",#"ward.D",#"complete", ,#"average",
          #row_km = 3,
          top_annotation = HeatmapAnnotation(Tier = bed[,8],
                                             col = list(Tier = c(`1` = "darkseagreen2",
                                                                 `2` = "darkseagreen4")),
                                             annotation_legend_param = list(Tier = list(title = "Tier",
                                                                                        title_gp = gpar(fontsize = 8),
                                                                                        labels_gp = gpar(fontsize = 8)))),
          # cluster_rows = row_dend,
          # draw(lgd,just = c("right", "top"))
          right_annotation = rowAnnotation(subtype = subtype,
                                           age = age,
                                           col = list(subtype = c("?" = "gray52","?Other" = "gray77",
                                                                  "EmbCarc" = "gold2",
                                                                  "NonSem" = "darkorange2",
                                                                  "Sem" = "deeppink3","SemClas" = "deeppink4",
                                                                  "Tera" = "darkslategray2","TeraUn" = "darkslategray3")),
                                           # semornon = c("Sem" = "deeppink3",
                                           #              "NonSem" = "darkorange2",
                                           #              "?" = "gray52","?Other" = "gray77")),
                                           annotation_legend_param = list(subtype = list(title = "subtype",
                                                                                         title_gp = gpar(fontsize = 10),
                                                                                         labels_gp = gpar(fontsize = 8))),
                                           # semornon = list(title = "Sem/NS",
                                           #                 title_gp = gpar(fontsize = 10),
                                           #                 labels_gp = gpar(fontsize = 8))),
                                           tmb = anno_barplot(tmb, baseline = "min")),
          left_annotation = rowAnnotation(semornon=semornon,
                                          type=type,
                                          ploidy=ploidy,
                                          col=list(semornon = c("Sem" = "deeppink3",
                                                                "NonSem" = "darkorange2",
                                                                "?" = "gray52","?Other" = "gray77"),
                                                   type = c("METASTASES" = "deeppink",
                                                            "PRIMARY" = "goldenrod4")),
                                          annotation_legend_param = list(semornon = list(title = "semornon",
                                                                                         title_gp = gpar(fontsize = 10),
                                                                                         labels_gp = gpar(fontsize = 8),
                                                                                         title = "Type",
                                                                                         title_gp = gpar(fontsize = 10),
                                                                                         labels_gp = gpar(fontsize = 8)))),
          column_names_gp = grid::gpar(fontsize = 8))

  dev.off()

}

