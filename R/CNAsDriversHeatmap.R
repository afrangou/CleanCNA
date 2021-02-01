# -------------------------- #
# Make heatmap
# -------------------------- #

#' @name
#' CNAsDriversHeatmapDendrogram
#'
#' @title
#' Heatmap to show CNA status in a set of driver genes across a cohort of tumour samples
#' Dendrogram to go alongside - heatmap ordered to match dendrogram
#'
#' @description
#' Create heatmap showing CNA status at driver genes using ggplot2
#' homdel, loh, otherloss, nochange, gain, biggain
#'
#' @param
#' segfile_dir = dir with all *_subclones.txt files are
#' segfile_name = label of cohort, eg 'TGCT'
#' bed_file = tab sep file, 7 cols, chr,start,end,strand,ens,ens,genename
#' driver_file = tab sep file, colnames SYMBOL (of gene) and TIER (1,2,etc) existing
#'
#' @return
#' Heatmap and associated dendrogram (in same order) of CNAs v drivers
#'

# collate all segments in subclones files across cohort
CNAsDriversHeatmapDendrogram <- function(segfile_dir,
                               segfile_name,
                               bed_file,
                               driver_file) {

    # libs
    library(dplyr)
    library(ggplot2)
    library(reshape)
    library(ggdendro)
    library(reshape2)
    library(grid)

    # subclones nf dir
    subclones_dir = segfile_dir
    subclones = system(paste0("ls ",subclones_dir,"*subclones.txt"),intern=T)
    samples = sapply(subclones,function(x){strsplit(x,"_subclones.txt")[[1]][1]})
    samples = sapply(samples,function(x){strsplit(x,subclones_dir)[[1]][2]})

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

    genes = genes[order(as.integer(genes$positions_chr)),]
    genes$positions_chr = as.integer(genes$positions_chr)
    genes$positions_start = as.integer(genes$positions_start)
    genes$positions_end = as.integer(genes$positions_end)

    # get subs files from first step
    subs = read.table(paste0(segfile_dir,segfile_name,"_segsfull.txt"),
                      sep="\t",
                      stringsAsFactors=F,
                      hea=T)

    subs$class[subs$class=="dip"]=2
    subs$class[subs$class=="tetra"]=4
    subs$class = as.integer(subs$class)
    subs$coded_total_cn=NA
    subs$total_cn = subs$nMajor + subs$nMinor

    # homdels
    subs$coded_total_cn[which(subs$total_cn==0)]= "homdel"
    # LOH
    subs$coded_total_cn[which(subs$nMajor>0 & subs$nMinor==0)] = "loh"
    # non homdel non LOH loss
    subs$coded_total_cn[which(subs$total_cn < subs$class
                              & subs$total_cn != 0
                              & subs$nMinor != 0)] = "otherloss"
    # no change
    subs$coded_total_cn[which((subs$nMajor==subs$nMinor | subs$nMajor==3 & subs$nMinor==1)
                              & (subs$total_cn == subs$class))] = "nochange"
    # gain
    subs$coded_total_cn[which(subs$total_cn>subs$class)] ="gain"
    # big gain
    subs$coded_total_cn[which(subs$total_cn>(5*subs$class))] = "biggain"
    # change X to 23
    subs$chr[subs$chr=="X"]=23

    # classify these in same way as was done for whole genome plot
    # homdel = read.csv(paste0(segfile_dir,"ENDOMETRIAL_homdel_subs.full"),
    #                   stringsAsFactors=F,
    #                   sep="\t")

    # make data table with rows = samples, columns = gene names of importance
    tomelt = matrix(NA,nrow=length(unique(subs$sample)),ncol=nrow(genes))
    rownames(tomelt) = unique(subs$sample)
    colnames(tomelt) = genes$positions_gene

    subs$sample = sapply(subs$sample,function(x){strsplit(x,"\\.")[[1]][2]})
    subs$chr=as.integer(subs$chr)
    for (i in 1:length(samples)) {
      sub = subs[subs$sample==samples[i],]
      for (g in 1:nrow(genes)) {
        chr = sub[sub$chr==genes$positions_chr[g],]
        row = which((chr$startpos <= genes$positions_start[g]) & (chr$endpos >= genes$positions_end[g]))
        if (length(row)==0) {
          minrows = max(which(chr$startpos <= genes$positions_start[g]))
          maxrows = min(which(chr$endpos >= genes$positions_end[g]))
          rows = minrows:maxrows
          longestregion = which(chr[rows,3]-chr[rows,2] == max(chr[rows,3]-chr[rows,2]))
          row = rows[longestregion]
          tomelt[i,g] = sub$coded_total_cn[row]
        } else {
          tomelt[i,g] = sub$coded_total_cn[row]
        }
      }
    }

    # save unmelted matrix to use for dendrogrm
    tomelt[tomelt=="homdel"]=0
    tomelt[tomelt=="loh"]=1
    tomelt[tomelt=="otherloss"]=2
    tomelt[tomelt=="nochange"]=3
    tomelt[tomelt=="gain"]=4
    tomelt[tomelt=="biggain"]=5
    for (i in 1:ncol(tomelt)) {tomelt[,i]=as.numeric(tomelt[,i])}
    write.table(tomelt,paste0(segfile_dir,segfile_name,"_dendrogram.csv"),
                              quote=F)
    # write melted table for heatmap
    melted = melt(tomelt)
    colnames(melted) = c("tumour", "driver","normalised_total_cn")
    write.table(melted,paste0(segfile_dir,segfile_name,"_heatmap.csv"),
                row.names=F,
                quote=F)


    # make dendrogram
    library(dendextend)
    todend <- as.data.frame(tomelt)
    for (i in 1:ncol(todend)) {todend[,i]=as.numeric(todend[,i])}
    todend.scaled <- scale(todend)
    # hierarchical clustering
    todend.scaled.clustered <- hclust(dist(todend.scaled))
    # make dendrogram
    todend.dend <-as.dendrogram(todend.scaled.clustered)
    # plot
    dendplot <- ggdendrogram(data = todend.dend, rotate=T,axis.text.y = element_text(size = 1))
    pdf(paste0(segfile_dir,segfile_name,"_dendrogram.pdf"))
      dendplot
    dev.off()

    # get order for heatmap from dendrogram
    heatmap.order = order.dendrogram(todend.dend)

    # make heatmap
    toplot = read.csv(paste0(segfile_dir,segfile_name,"_heatmap.csv"),
                      stringsAsFactors=F,
                      sep=" ")

    toplot$normalised_total_cn = as.factor(toplot$normalised_total_cn)
    # label genes with a factor based on order in genome so plots correctly
    toplot$order = as.character(sort(rep(1:length(unique(toplot$driver)),length(unique(toplot$tumour)))))
    toplot$driver <- factor(toplot$driver, levels=unique(toplot$driver)[order(toplot$order)])
    # label tumours with a factor based on order of dendrogram so plots correctly
    toplot = cbind(toplot,order2=rep(1:length(unique(toplot$tumour)),length(unique(toplot$driver))))
    toplot = cbind(toplot,use=match(toplot$order2,heatmap.order))
    toplot = toplot[order(as.integer(toplot$order),toplot$use),]
    toplot$tumour <- factor(toplot$tumour, levels = unique(toplot$tumour)[order(as.numeric(toplot$use))])

    plot1 <- ggplot(toplot,aes(driver, tumour, fill = normalised_total_cn)) +
      geom_tile(size = 0.2) +#, color = "white") +
      theme(text = element_text(size=8),
            axis.title.x = element_blank(),
            axis.title.y=element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle = 80, vjust = 0.5),
            panel.border = element_blank()) +
      scale_fill_manual(values=c("black","deepskyblue4","deepskyblue2", "darkolivegreen2", "deeppink","deeppink4")) +
      theme(legend.position = "left")

    pdf(paste0(segfile_dir,segfile_name,"_heatmap.pdf"))
      plot1
    dev.off()

}
