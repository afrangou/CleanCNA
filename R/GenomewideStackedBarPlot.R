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

  library(scales)
  library(ggplot2)

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

  # make plot
  pdf(paste0(filestub,segfile_name,"_genomewide_stacked_bar.pdf"),
      width=10,height=6)
    ggplot(toplot) +
      # above x axis
      geom_rect(aes(ymin = 0, xmin = posleft, xmax = posright, ymax = value, fill = nochangecol)) +
      geom_rect(aes(ymin = value, xmin = posleft, xmax = posright, ymax = value2, fill = gaincol)) +
      geom_rect(aes(ymin = value2, xmin = posleft, xmax = posright, ymax = value3, fill = biggainpcol)) +
      # geom_rect(aes(ymin = value3, xmin = posleft, xmax = posright, ymax = value4, fill = pentapcol)) +
      # below x axis
      geom_rect(aes(ymin = 0, xmin = posleft, xmax = posright, ymax = valueneg, fill = otherlosscol)) +
      geom_rect(aes(ymin = 0, xmin = posleft, xmax = posright, ymax = valueneg, fill = lohcol)) +
      geom_rect(aes(ymin = valueneg, xmin = posleft, xmax = posright, ymax = valueneg2, fill = homdelcol)) +
      # ylabel
      ylab("% tumour samples") +
      # split into chrs
      facet_grid(~poschr,scales = "free_x",space="free_x") +
      # remove x axis labels
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      # control y axis scaling
      coord_cartesian(ylim=c(-100,100)) +
      scale_y_continuous(breaks = c(50, 100, -50, -100, 0),
                         labels = c(50, 100, 50, 100, 0)) +
      # add black line for chr and dot for centromere
      scale_x_continuous(breaks = NULL) +
      geom_point(data=centro, aes(x=pos,y=0), size=1.5) +
      geom_hline(yintercept=0, colour="black", size=0.5) +
      # add legend
      scale_fill_manual(values = c(nochangecol,gaincol,biggainpcol,otherlosscol,lohcol,homdelcol),
                        labels =c("nochange","gain","biggain","otherloss","loh","homdel"))

  dev.off()

}

