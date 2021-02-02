# -------------------------- #
# Prepare data for stacked barplot
# -------------------------- #

#' @name
#' GenomewideStackedBarPlot
#'
#' @title
#' Create stacked bar plot
#'
#' @description
#' Create stacked barplot for 6 different CNA types
#' homdel, loh, otherloss, nochange, gain, biggain, split into loss, no change, and gain
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
GenomewideStackedBarPlot <- function(filestub,
                                 segfile_name,
                                 number_samples) {
  # load libs
  library(scales)
  library(ggplot2)

  # load file
  toplot = read.table(paste0(filestub,segfile_name,"_partitioned_regions_all_CNA_types_overlapped.out"),
                      hea=T,
                      stringsAsFactors=F,
                      fill=T)

  # centromeres for hg38
  centro = data.frame(pos=c(123400000, 93900000, 90900000,
                            50000000, 48750000,  60550000,
                            60100000, 45200000, 43850000,
                            39800000, 53400000, 35500000,
                            17700000, 17150000, 19000000,
                            36850000, 25050000, 18450000,
                            26150000, 28050000, 11950000, 15550000),poschr=c(1:22))

  # colours
  colours <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
  nochangecol = colours[3]
  gaincol = colours[2]
  biggainpcol = colours[1]
  otherlosscol = colours[6]
  lohcol = colours[5]
  homdelcol = colours[4]

  # make totals of nochange, loss, gain, for separate plots
  toplot$valuenochange = toplot$nochange
  toplot$valuegain = toplot$gain
  toplot$valuegain2 = toplot$gain + toplot$biggain
  toplot$valueloss = toplot$otherloss
  toplot$valueloss2 = toplot$otherloss + toplot$loh
  toplot$valueloss3 = toplot$otherloss + toplot$loh + toplot$homdel

  # change to proportions
  for (column in 4:15) {toplot[,column]=(toplot[,column]/number_samples)*100}

  # order chromosomes numerically
  toplot$poschr = as.integer(gsub("chr","",toplot$poschr))
  all = toplot[which(toplot$poschr==1),]
  for (chr in 2:23) {

    sub = toplot[which(toplot$poschr==chr),]
    all = rbind(all,sub)

  }
  toplot=all

  # gains plot
  gainstoplot = toplot[,c(1:3,8:9,11:12),]

  pdf(paste0(filestub,segfile_name,"_genomewide_stackedbar_gains.pdf"),
      width=10,height=2)

    ggplot(gainstoplot) +
      # above x axis
      geom_rect(aes(ymin = 0, xmin = posleft, xmax = posright, ymax = valuegain, fill = gaincol)) +
      geom_rect(aes(ymin = valuegain, xmin = posleft, xmax = posright, ymax = valuegain2, fill = biggainpcol)) +
      # ylabel
      ylab("% tumour samples") +
      # split into chrs
      facet_grid(~poschr,scales = "free_x",space="free_x") +
      # remove x axis labels
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      # control y axis scaling
      coord_cartesian(ylim=c(0,50)) +
      scale_y_continuous(breaks = c(0, 10, 20, 30, 40,50),
                         labels = c(0, 10, 20, 30, 40,50)) +
      #add black line for chr and dot for centromere
      scale_x_continuous(breaks = NULL) +
      geom_point(data=centro, aes(x=pos,y=0), size=1.5) +
      geom_hline(yintercept=0, colour="black", size=0.5) +
      # add legend
      scale_fill_manual(values = c(gaincol,biggainpcol),
                        labels =c("gain","biggain"))

  dev.off()


  # losses plot
  lossestoplot = toplot[,c(1:6,13:15),]

  pdf(paste0(filestub,segfile_name,"_genomewide_stackedbar_losses.pdf"),
      width=10,height=2)
  ggplot(lossestoplot) +

    # above x axis
    geom_rect(aes(ymin = 0, xmin = posleft, xmax = posright, ymax = valueloss, fill = otherlosscol)) +
    geom_rect(aes(ymin = valueloss, xmin = posleft, xmax = posright, ymax = valueloss2, fill = lohcol)) +
    geom_rect(aes(ymin = valueloss2, xmin = posleft, xmax = posright, ymax = valueloss3, fill = homdelcol)) +
    # ylabel
    ylab("% tumour samples") +
    # split into chrs
    facet_grid(~poschr,scales = "free_x",space="free_x") +
    # remove x axis labels
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    # control y axis scaling
    coord_cartesian(ylim=c(0,50)) +
    scale_y_continuous(breaks = c(0, 10, 20, 30, 40,50),
                       labels = c(0, 10, 20, 30, 40,50)) +
    #add black line for chr and dot for centromere
    scale_x_continuous(breaks = NULL) +
    geom_point(data=centro, aes(x=pos,y=0), size=1.5) +
    geom_hline(yintercept=0, colour="black", size=0.5) +
    # add legend
    scale_fill_manual(values = c(otherlosscol,lohcol,homdelcol),
                      labels =c("homdel","loh","otherloss"))
  dev.off()



  # no change plot
  nochangetoplot = toplot[,c(1:3,7,10),]
  nochangetoplot$value = nochangetoplot$nochange

  pdf(paste0(filestub,segfile_name,"_genomewide_stackedbar_nochange.pdf"),
      width=10,height=2)
  ggplot(nochangetoplot) +
    # above x axis
    geom_rect(aes(ymin = 0, xmin = posleft, xmax = posright, ymax = valuenochange, fill = nochangecol)) +
    # geom_rect(aes(ymin = valueneg, xmin = posleft, xmax = posright, ymax = valueneg2, fill = lohcol)) +
    # geom_rect(aes(ymin = valueneg2, xmin = posleft, xmax = posright, ymax = valueneg3, fill = homdelcol)) +
    # below x axis
    # geom_rect(aes(ymin = 0, xmin = posleft, xmax = posright, ymax = valueneg, fill = otherlosscol)) +
    # geom_rect(aes(ymin = 0, xmin = posleft, xmax = posright, ymax = valueneg, fill = lohcol)) +
    # geom_rect(aes(ymin = valueneg, xmin = posleft, xmax = posright, ymax = valueneg2, fill = homdelcol)) +
    # ylabel
    ylab("% tumour samples") +
    # split into chrs
    facet_grid(~poschr,scales = "free_x",space="free_x") +
    # remove x axis labels
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    # control y axis scaling
    coord_cartesian(ylim=c(0,100)) +
    scale_y_continuous(breaks = c(0, 50,100),
                       labels = c(0, 50,100)) +
    #add black line for chr and dot for centromere
    scale_x_continuous(breaks = NULL) +
    geom_point(data=centro, aes(x=pos,y=0), size=1.5) +
    geom_hline(yintercept=0, colour="black", size=0.5) +
    # add legend
    scale_fill_manual(values = nochangecol,
                      labels ="no change")

  dev.off()


}

