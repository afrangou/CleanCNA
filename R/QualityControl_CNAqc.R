# complete QC of Battenberg output
# written by Anna Frangou and Alex J. Cornish
#
# usage:
# Rscript QC.R [RUN] [DIR_BATTENBERG] [SUBDIR_CALLSUBCLONES] [SUBDIR_POSTPROCESSING] [SUBDIR_DPCLUST] [SUBDIR_VAFPEAKS] [FILENAME_SAMPLELIST] [FILENAME_CONFIG] [FILENAME_CHRSIZES] [FILENAME_QC]
#
#### subroutines ####
pasteu <- function(...) paste(..., sep="_")

pastedot <- function(...) paste(..., sep=".")

log.message <- function(x, level=0) message(paste0("[", Sys.time(), "] ", paste(rep(" ", times=level * 4), collapse=""), "- ", x))

# estimate.new.ploidy <- function(rho.old, psi.old, rho.new) (rho.old * (psi.old - 2) + 2 * rho.new) / rho.new # from eq. S5, Van Loo et al. (PNAS, 2010)

estimate.new.ploidy <- function(rho.old, psi.old, rho.new) {
  if (psi.old == "tetra") {ploidytype = "tetra"} else {ploidytype = "dip"}
  if (ploidytype=="dip") {psi.new = ((rho.old * psi.old) + 2*(rho.new - rho.old)) / rho.new}
  if (ploidytype=="tetra") {psi.new = ((rho.old * psi.old) + 4*(rho.new - rho.old)) / rho.new}
  return(psi.new)
}

switch.ploidy.get.purity <- function(rho.old,psi.old) {
  if (psi.old == "tetra") {ploidytype = "tetra"} else {ploidytype = "dip"}
  if (ploidytype == "tetra") {rho.new = 2*rho.old / (rho.old+1)}
  if (ploidytype == "dip")  {rho.new = rho.old / (2-rho.old)}
  return(rho.new)
}

switch.ploidy.get.ploidy <- function(rho.old,psi.old) {
  if (psi.old == "tetra") {ploidytype = "tetra"} else {ploidytype = "dip"}
  if (ploidytype == "tetra") {psi.new = psi.old/2}
  if (ploidytype == "dip")  {psi.new = psi.old*2}
  return(psi.new)
}

estimate.dpclust.purity <- function(clusters, rho.old) {
  # identify clonal cluster as cluster with highest CCF (noisy clusters should have already been removed)
  # reestimate purity using position of this clonal cluster
  # if this new purity is >1, exclude cluster and retry
  # if purity >1 after all clusters considered, return old purity
  if (rho.old > 1) stop("old rho should be <1")
  consider <- rep(T, nrow(clusters))
  rho.new <- Inf
  while (rho.new > 1) {
    if (sum(consider)) {
      # select cluster with greatest CCF to compute new purity
      clonalcluster <- which(clusters$location == max(clusters$location[consider]))[1]
      rho.new <- rho.old * clusters$location[clonalcluster] # compute new purity
      consider[clonalcluster] <- F # ensure this cluster is not considered in any later rounds
    } else {
      # if rho.new is >1 and there are no more clusters to consider, return old purity
      warning("failed to compute new purity")
      rho.new <- rho.old
    }
  }
  rho.new
}

get.worst.filter <- function(filter.list) {
  # identify the worst filter result from a list of filter results (FAIL worse than FLAG worse than PASS)
  worst.filter <- "PASS"
  if ("FLAG" %in% filter.list) worst.filter <- "FLAG"
  if ("FAIL" %in% filter.list) worst.filter <- "FAIL"
  worst.filter
}

larger.all.later.numbers <- function(x) sapply(1:length(x), function(i) ifelse(i == length(x), T, x[i] > max(x[(i+1):length(x)])))

qc_CNAqc <- function(
  run,
  dir.battenberg,
  subdir.callsubclones,
  subdir.postprocessing,
  subdir.dpclust,
  subdir.vafpeaks,
  filename.samplelist,
  filename.config,
  filename.chr.sizes,
  filename.qc
) {
  # complete QC

  # required libraries
  library(stringr)

  source(filename.config) # source config file
  time.pretty <- str_replace_all(str_replace_all(Sys.time(), " ", "_"), ":", "-")
  run.name <- paste0("run", run) # name of this run

  log.message("reading in data")
  chr.sizes.raw <- read.delim(filename.chr.sizes, sep="\t", header=F, stringsAsFactors=F)
  chr.sizes <- structure(chr.sizes.raw[[2]], names=chr.sizes.raw[[1]])
  samplelist <- read.delim(filename.samplelist, sep="\t", header=T, stringsAsFactors=F)
  print(rownames(samplelist))
  rownames(samplelist) <- ids <- paste0("tumo", samplelist$tumour_sample_platekey, "_norm", samplelist$germline_sample_platekey)
  print(ids)

  # directories and filenames
  filenames.battenberg.segs <- sapply(ids, function(id) file.path(
    dir.battenberg,"intermediateFiles",samplelist[id, "participant_id"],id,
    subdir.postprocessing,
    pasteu(samplelist[id, "tumour_sample_platekey"], "subclones.txt")))
  filenames.battenberg.purity <- sapply(ids, function(id) file.path(
    dir.battenberg,"intermediateFiles",samplelist[id, "participant_id"],id,
    subdir.callsubclones,
    pasteu(samplelist[id, "tumour_sample_platekey"], "cellularity_ploidy.txt")))
  filenames.battenberg.plot <- sapply(ids, function(id) file.path(
    dir.battenberg,"intermediateFiles",samplelist[id, "participant_id"],id,
    subdir.callsubclones,
    pasteu(samplelist[id, "tumour_sample_platekey"], "BattenbergProfile_average.pdf")))
  filenames.dpclust.optima <- sapply(ids, function(id) file.path(
    dir.battenberg,"intermediateFiles",samplelist[id, "participant_id"],id,
    subdir.dpclust,
    pasteu(samplelist[id, "tumour_sample_platekey"], "optimaInfo.txt")))
  filenames.peaks <- sapply(ids, function(id) file.path(
    dir.battenberg,"intermediateFiles",samplelist[id, "participant_id"],id,
    subdir.vafpeaks,
    pastedot(samplelist[id, "participant_id"], id, "peak_data.RData")))

  log.message("reading in sample-level data")
  # gets IDs for samples that are being QC'd
  ids <- ids[file.exists(filenames.battenberg.purity) & file.exists(filenames.battenberg.segs) & file.exists(filenames.dpclust.optima)]
  # reformat subclones file
  segs <- sapply(ids, function(id) read.delim(filenames.battenberg.segs[[id]], sep="\t", header=T, stringsAsFactors=F)[,1:13], simplify=F)
  # add 'chr'
  for (id in ids) segs[[id]]$chr <- paste0("chr", segs[[id]]$chr)
  # get battenberg purity from cellularity_ploidy file
  battenberg.purity <- do.call("rbind", sapply(ids, function(id) read.delim(filenames.battenberg.purity[id], sep="\t", header=T, stringsAsFactors=F), simplify=F))
  # get dpclust optimaInfo file
  dpclust <- sapply(ids, function(id) read.delim(filenames.dpclust.optima[id], sep="\t", header=T, stringsAsFactors=F), simplify=F)
  # calculate number of mutations for this sample that constitutes the threshold for superclonal clusters or 0.5 clusters during qc
  fivepcclusters <- sapply(dpclust,function(x) round(sum(x$no.of.mutations) * thres.propmuts.superclonal.or.tetra))

  # remove clusters smaller than thres.propmuts (% of total mutations)
  log.message("removing small clusters from DPClust optima")
  dpclust <- sapply(dpclust, function(x) x[x$no.of.mutations > (sum(x$no.of.mutations) * thres.propmuts), ], simplify=F)

  # read in or setup QC file - this contains info for all runs, for all samples on the sample list
  if (run > 1) {
    log.message("reading existing QC file")
    qc <- read.delim(filename.qc, sep="\t", header=T, stringsAsFactors=F)
    rownames(qc) <- qc$tumour_normal_id
  } else {
    log.message("creating QC file")
    qc <- samplelist[ids, c("participant_id", "sex", "tumour_sample_platekey", "germline_sample_platekey", "somatic_small_variants_qced_vcf")]
    qc$tumour_normal_id <- ids
  }

  # if not the first run, filter down to only deal with samples that failed previous runs
  if (sum(grepl(paste0(run.name, "_"), colnames(qc)))) warning("columns already found for run ", run)
  if (run > 1) ids <- ids[which(qc[ids, pasteu(paste0("run", run - 1), "filter_overallfilter")] == "FAIL")] # identify samples not passing previous run
  if (length(ids) == 0) stop("no tumour-normal pairs to QC")

  # add purity and ploidy estimates from BB (cellularity file) and dpclust (optimaInfo)
  log.message("add purity and ploidy estimates")
  qc[ids, pasteu(run.name, "battenberg_purity")] <- battenberg.purity[ids, ifelse("cellularity" %in% colnames(battenberg.purity),'cellularity','purity')] 
  qc[ids, pasteu(run.name, "battenberg_ploidy")] <- battenberg.purity[ids, "psi"]
  qc[ids, pasteu(run.name, "dpclust_purity")] <- sapply(ids, function(id) estimate.dpclust.purity(dpclust[[id]], qc[id, pasteu(run.name, "battenberg_purity")]))
  qc[ids, pasteu(run.name, "dpclust_ploidy")] <- sapply(ids, function(id) estimate.new.ploidy(qc[id, pasteu(run.name, "battenberg_purity")], qc[id, pasteu(run.name, "battenberg_ploidy")], qc[id, pasteu(run.name, "dpclust_purity")]))
  qc[ids, pasteu(run.name, "fivepc_cluster_size")] <- fivepcclusters[ids]
  for (str in c("vafpeaks_purity", "vafpeaks_ploidy", "vafpeaks_score", "reestimated_purity", "reestimated_ploidy")) qc[[pasteu(run.name, str)]] <- NA

  # compute statistics required for QC for each sample
  log.message("computing QC statistics")
  for (id in ids) {
    log.message(paste0("sample ", which(ids == id), "/", length(ids)), level=1)

    # add segment sizes, extract autosome segments and get segment copy number in major and minor subclones
    segs[[id]]$size <- segs[[id]]$endpos - segs[[id]]$startpos + 1
    segs.auto <- segs[[id]][segs[[id]]$chr %in% autosomes, ]
    # remove segments with missing major/minor copy number (should only be telomeres and maybe centromeres)
    segs.auto <- segs.auto[!is.na(segs.auto$nMaj1_A) & !is.na(segs.auto$nMin1_A), ]
    # put copy number of largest subclone first in the segs table
    flip.clones <- rep(F, nrow(segs.auto))
    consider.flip <- !is.na(segs.auto$frac2_A)
    flip.clones[consider.flip] <- segs.auto$frac1_A[consider.flip] < 0.5
    segs.auto[flip.clones, c("nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")] <- segs.auto[flip.clones, c("nMaj2_A", "nMin2_A", "frac2_A", "nMaj1_A", "nMin1_A", "frac1_A")]
    # total copy number for each subclone
    segs.auto$ctClone1 <- segs.auto$nMaj1_A + segs.auto$nMin1_A
    segs.auto$ctClone2 <- segs.auto$nMaj2_A + segs.auto$nMin2_A

    # identify failing chromosomes - ie those with less than 50% of the chromosome with a CNV call
    chrs.fail <- list()
    chrs.fail$chrmissing <- sort(setdiff(autosomes, segs.auto$chr))
    chrs.fail$chrsizeincorrect <- chrs[(sapply(autosomes, function(chr) sum(segs.auto$size[segs.auto$chr == chr])) / chr.sizes[autosomes]) < thres.chrsizeincorrect.tol]

    # add information on which chromosomes fail and number of failures
    for (condition in names(chrs.fail)) {
      qc[id, pasteu(run.name, "which", condition)] <- paste(chrs.fail[[condition]], collapse=";")
      qc[id, pasteu(run.name, "n", condition)] <- length(chrs.fail[[condition]])
    }

    # identify segments satisfying various conditions used to apply filters
    segs.conditions <- list()

    segs.conditions$all <- segs.auto # all segments, needs to be the first condition

    segs.conditions$copynumber11 <- segs.auto[
      (segs.auto$nMaj1_A == 1 & segs.auto$nMin1_A == 1) &
        segs.auto$frac1_A == 1,
      ] # copy number is 1:1 and there is no evidence of sub-clonal change

    segs.conditions$copynumber10 <- segs.auto[
      ((segs.auto$nMaj1_A == 1 & segs.auto$nMin1_A == 0) | (segs.auto$nMaj1_A == 0 & segs.auto$nMin1_A == 1)) &
        segs.auto$frac1_A == 1,
      ] # copy number is 1:0 and there is no evidence of sub-clonal change

    #minorCNround = round(segs.auto$nMin1_A*segs.auto$frac1_A+segs.auto$nMin2_A*segs.auto$frac2_A)
    segs.auto2 = segs.auto
    segs.auto2[is.na(segs.auto2)] = 0
    segs.conditions$minorCNzero <- segs.auto2[
      (which(round((segs.auto2$nMin1_A*segs.auto2$frac1_A)+(segs.auto2$nMin2_A*segs.auto2$frac2_A))==0)),
      ] # minor CN of 0 after rounding copy number, allows us to define a sample as diploid or tetraploid

    segs.conditions$copynumber21 <- segs.auto[
      ((segs.auto$nMaj1_A == 2 & segs.auto$nMin1_A == 1) | (segs.auto$nMaj1_A == 1 & segs.auto$nMin1_A == 2)) &
        segs.auto$frac1_A == 1,
      ] # copy number is 2:1 and there is no evidence of sub-clonal change

    segs.conditions$copynumber22 <- segs.auto[
      (segs.auto$nMaj1_A == 2 & segs.auto$nMin1_A == 2) &
        segs.auto$frac1_A == 1,
      ] # copy number is 2:2 and there is no evidence of sub-clonal change

    segs.conditions$copynumber32 <- segs.auto[
      ((segs.auto$nMaj1_A == 3 & segs.auto$nMin1_A == 2) | (segs.auto$nMaj1_A == 2 & segs.auto$nMin1_A == 3)) &
        segs.auto$frac1_A == 1,
      ] # copy number is 3:2 and there is no evidence of sub-clonal change

    segs.conditions$homodelall <- segs.auto[
      segs.auto$ctClone1 == 0 |
        sapply(as.double(segs.auto$ctClone2), identical, 0.0),
      ] # homozygous deletion of either major or (if evidence of sub-clonality) minor clone

    segs.conditions$homodellargest <- if (nrow(segs.conditions$homodelall)) segs.conditions$homodelall[which(segs.conditions$homodelall$size == max(segs.conditions$homodelall$size))[1], ] else segs.conditions$homodelall

    segs.conditions$cnodd <- segs.auto[
      (((segs.auto$nMaj1_A == 2 & segs.auto$nMin1_A == 1) | (segs.auto$nMaj1_A == 1 & segs.auto$nMin1_A == 2)) |
         ((segs.auto$nMaj1_A == 1 & segs.auto$nMin1_A == 0) | (segs.auto$nMaj1_A == 0 & segs.auto$nMin1_A == 1)) |
         ((segs.auto$nMaj1_A == 3 & segs.auto$nMin1_A == 2) | (segs.auto$nMaj1_A == 2 & segs.auto$nMin1_A == 3))) &
        segs.auto$frac1_A == 1,
      ] # copy number is 2:1 or 1:0 or 3:2 and there is no evidence of sub-clonal change

    segs.conditions$aroundpoint5.narrow <- segs.auto[
      segs.auto$frac1_A <= thres.50pcpeak.upper & segs.auto$frac1_A >= thres.50pcpeak.lower,
      ] # copy number fraction around 0.5, narrow boundary

    segs.conditions$aroundpoint5.wide <- segs.auto[
      segs.auto$frac1_A <= thres.50pcpeak.upper.wide & segs.auto$frac1_A >= thres.50pcpeak.lower.wide,
      ] # copy number fraction around 0.5, wide boundary

    # add information on the number of segments, length of segments and fraction of considered genome segments represent
    for (condition in names(segs.conditions)) {
      qc[id, pasteu(run.name, "nsegs", condition)] <- nrow(segs.conditions[[condition]])
      qc[id, pasteu(run.name, "lsegs", condition)] <- sum(segs.conditions[[condition]]$size)
      qc[id, pasteu(run.name, "fgenome", condition)] <- qc[id, pasteu(run.name, "lsegs", condition)] / sum(segs.conditions$all$size)
    }

    # add tumour ploidy, ie psi_t (different to average ploidy)
    qc[id, pasteu(run.name, "psi_t")] <- sum(segs.auto$size * (segs.auto$nMaj1_A + segs.auto$nMin1_A))/sum(segs.auto$size)
    # add fraction of LOH *-2 +2.9 (from PCAWG equation)
    qc[id, pasteu(run.name, "minorCNzerotimesminus2plus2.9")] <- (qc[id,pasteu(run.name,"fgenome_minorCNzero")]*-2)+2.9

    # define sample as 'narrow' if there is one cluster within 0.95-1.05 and 0.45-0.55 (or as per config file)
    # or 'wide' if there are two clusters in either of those boundaries
    # classify sample as narrow or wide
    qc[id, pasteu(run.name, "narrow.or.wide")] <- ifelse(((sum(dpclust[[id]]$location >= thres.clonalpeak.lower.wide &
                                                                 dpclust[[id]]$location <= thres.clonalpeak.upper.wide)>1) |
                                                            (sum(dpclust[[id]]$location >= thres.50pcpeak.lower.wide &
                                                                   dpclust[[id]]$location <= thres.50pcpeak.upper.wide)>1)),"wide","narrow")


    # classify sample as currently diploid or tetraploid
    # if qc$psi_t > qc$minorCNzerotimesminus2plus2.9, then tetraploid, otherwise diploid
    qc[id, pasteu(run.name, "dip.or.tetra")] <- ifelse(qc[id,pasteu(run.name,"psi_t")] > qc[id,pasteu(run.name,"minorCNzerotimesminus2plus2.9")],
                                                       "tetra",
                                                       "dip")

    # peak closest to clonal
    qc[id, pasteu(run.name, "peak.closest.to.clonal")] <- min(abs(1-dpclust[[id]]$location))

    # add DPClust peak subclonal, superclonal, 50% location information,
    # this is dependent on whether the sample is defined as narrow or wide
    # only consider clusters comprising a certain proportion of mutations
    qc[id, pasteu(run.name, "nsubclonalpeaks")] <- ifelse(qc[id,pasteu(run.name,"narrow.or.wide")]=="narrow",
                                                          sum(dpclust[[id]]$location < thres.clonalpeak.lower),
                                                          sum(dpclust[[id]]$location < thres.clonalpeak.lower.wide))
    qc[id, pasteu(run.name, "nclonalpeaks")] <- ifelse(qc[id,pasteu(run.name,"narrow.or.wide")]=="narrow",
                                                       sum(dpclust[[id]]$location >= thres.clonalpeak.lower &
                                                             dpclust[[id]]$location <= thres.clonalpeak.upper),
                                                       sum(dpclust[[id]]$location >= thres.clonalpeak.lower.wide &
                                                             dpclust[[id]]$location <= thres.clonalpeak.upper.wide))
    qc[id, pasteu(run.name, "nsuperclonalpeaks")] <- ifelse(qc[id,pasteu(run.name,"narrow.or.wide")]=="narrow",
                                                            sum(dpclust[[id]]$location > thres.clonalpeak.upper &
                                                                  dpclust[[id]]$no.of.mutations > qc[id, pasteu(run.name, "fivepc_cluster_size")]),
                                                            sum(dpclust[[id]]$location > thres.clonalpeak.upper.wide &
                                                                  dpclust[[id]]$no.of.mutations > qc[id, pasteu(run.name, "fivepc_cluster_size")]))
    qc[id, pasteu(run.name, "n50pcpeaks")] <- ifelse(qc[id,pasteu(run.name,"narrow.or.wide")]=="narrow",
                                                     sum(dpclust[[id]]$location > thres.50pcpeak.lower &
                                                           dpclust[[id]]$location < thres.50pcpeak.upper &
                                                           dpclust[[id]]$no.of.mutations > qc[id, pasteu(run.name, "fivepc_cluster_size")]),
                                                     sum(dpclust[[id]]$location > thres.50pcpeak.lower.wide  &
                                                           dpclust[[id]]$location < thres.50pcpeak.upper.wide &
                                                           dpclust[[id]]$no.of.mutations > qc[id, pasteu(run.name, "fivepc_cluster_size")]))

    # add number of mutations in each of the subclonal, clonal, and superclonal categories
    qc[id, pasteu(run.name, "nmutssubclonal")] <-	ifelse(qc[id,pasteu(run.name,"narrow.or.wide")]=="narrow",
                                                         sum(dpclust[[id]]$no.of.mutations[dpclust[[id]]$location < thres.clonalpeak.lower]),
                                                         sum(dpclust[[id]]$no.of.mutations[dpclust[[id]]$location < thres.clonalpeak.lower.wide]))


    qc[id, pasteu(run.name, "nmutsclonal")] <- ifelse(qc[id,pasteu(run.name,"narrow.or.wide")]=="narrow",
                                                      sum(dpclust[[id]]$no.of.mutations[dpclust[[id]]$location >= thres.clonalpeak.lower &
                                                                                          dpclust[[id]]$location <= thres.clonalpeak.upper]),
                                                      sum(dpclust[[id]]$no.of.mutations[dpclust[[id]]$location >= thres.clonalpeak.lower.wide &
                                                                                          dpclust[[id]]$location <= thres.clonalpeak.upper.wide]))

    qc[id, pasteu(run.name, "nmutssuperclonal")] <- ifelse(qc[id,pasteu(run.name,"narrow.or.wide")]=="narrow",
                                                           sum(dpclust[[id]]$no.of.mutations[dpclust[[id]]$location > thres.clonalpeak.upper &
                                                                                               dpclust[[id]]$no.of.mutations > qc[id, pasteu(run.name, "fivepc_cluster_size")]]),
                                                           sum(dpclust[[id]]$no.of.mutations[dpclust[[id]]$location > thres.clonalpeak.upper.wide &
                                                                                               dpclust[[id]]$no.of.mutations > qc[id, pasteu(run.name, "fivepc_cluster_size")]]))

                                               
    # uses VAF peak data calculated in previous step (peaks)
    load(filenames.peaks[id])                                      
    #expected.locs <- peaks$summary$expected.locs
    expected.locs <- peaks$peaks_analysis$matches
    if (is.null(nrow(expected.locs))) {
      # if there are no matched peaks then 2:2 peak matching does not occur
      if (peaks$peaks_analysis$QC=="FLAG") {
        qc[id, pasteu(run.name, "vafpeaks_tetraploid")] <- TRUE  
      } else {
        qc[id, pasteu(run.name, "vafpeaks_tetraploid")] <- FALSE
      }
    } else {
      # select peaks to consider and determine whether they are within the threshold and purity distance off
      #expected.locs <- expected.locs[order(abs(expected.locs$peak.diff)), ]
      expected.locs <- expected.locs[order(abs(expected.locs$offset)), ]

      # find mutations in 2:2 regions with multiplicity 1 and multiplicity 2
      mult1 <- which(expected.locs$karyotype=="2:2" & expected.locs$mutation_multiplicity==1)
      #which(expected.locs$state == "2:2" & expected.locs$multiplicity == 1)
      mult2 <- which(expected.locs$karyotype=="2:2" & expected.locs$mutation_multiplicity==2)
      # if these mutations exist, the position of this peak must be within a certain distance of expected peak
      # in order to qualify as these 2:2 peaks existing - this is evidence of proper tetraploidy as there are
      # mutations at multiplicity 1 and 2 in 2:2 regions
      mult1.match <- ifelse(length(mult1)>0, abs(expected.locs$offset[mult1]) < thres.incorrecttetraploid.peak.diff, FALSE)
      mult2.match <- ifelse(length(mult2)>0, abs(expected.locs$offset[mult2]) < thres.incorrecttetraploid.peak.diff, FALSE)
      qc[id, pasteu(run.name, "vafpeaks_tetraploid")] <- mult1.match & mult2.match

      # add other information from vafpeaks
      #qc[id, pasteu(run.name, "vafpeaks_purity")] <- peaks$summary$purity.new
      qc[id, pasteu(run.name, "vafpeaks_purity")] <- qc[id, pasteu(run.name, "battenberg_purity")]-peaks$peaks_analysis$score
      # if this purity goes below 5%, fix it at 5%, if it goes above 95%, fix it at 95%
      if (qc[id, pasteu(run.name, "vafpeaks_purity")]<0.05) {qc[id, pasteu(run.name, "vafpeaks_purity")]=0.05}
      if (qc[id, pasteu(run.name, "vafpeaks_purity")]>0.95) {qc[id, pasteu(run.name, "vafpeaks_purity")]=0.95}
      qc[id, pasteu(run.name, "vafpeaks_ploidy")] <- estimate.new.ploidy(qc[id, pasteu(run.name, "battenberg_purity")],
                                                                         qc[id, pasteu(run.name, "battenberg_ploidy")],
                                                                         qc[id, pasteu(run.name, "vafpeaks_purity")])
      #qc[id, pasteu(run.name, "vafpeaks_score")] <- peaks$summary$eta
      qc[id, pasteu(run.name, "vafpeaks_score")] <-peaks$peaks_analysis$score
      # add in peaks$peaks_analysis$QC to qc table
      qc[id, pasteu(run.name, "vafpeaks")] <-peaks$peaks_analysis$QC
    }


    # give a call for whether the ploidy is judged to be correct by our criteria
    # if diploid, else if tetraploid
    #qc[id, pasteu(run.name, "ploidy_type_accepted")] <- if (qc[id,pasteu(run.name,"battenberg_ploidy")]<=2.6) {
    qc[id, pasteu(run.name, "ploidy_type_accepted")] <- if (qc[id,pasteu(run.name,"dip.or.tetra")]=="dip") {
      if (qc[id,pasteu(run.name,"narrow.or.wide")]=="narrow") {
        ifelse(qc[id,pasteu(run.name,"fgenome_aroundpoint5.narrow")] <= thres.incorrectdiploid.aroundpoint5 |
                 qc[id,pasteu(run.name,"n50pcpeaks")]==0,"PASS","FAIL")
      } else if (qc[id,pasteu(run.name,"narrow.or.wide")]=="wide") {
        ifelse(qc[id,pasteu(run.name,"fgenome_aroundpoint5.wide")] <= thres.incorrectdiploid.aroundpoint5 |
                 qc[id,pasteu(run.name,"n50pcpeaks")]==0,"PASS","FAIL")
      }
    #} else if (qc[id,pasteu(run.name,"battenberg_ploidy")] > 2.6) {
    } else if (qc[id,pasteu(run.name,"dip.or.tetra")] == "tetra") {
      ifelse(qc[id,pasteu(run.name,"fgenome_cnodd")] >= thres.incorrecttetraploid.cnodd,"PASS","FAIL")
    }

    # calculate the new ploidy and purity for IF the ploidy is deemed incorrect - not necessarily used
    qc[id, pasteu(run.name, "newploidy_purity")] <- switch.ploidy.get.purity(qc[id, pasteu(run.name,"battenberg_purity")],
                                                                             qc[id, pasteu(run.name,"battenberg_ploidy")])


    qc[id, pasteu(run.name, "newploidy_ploidy")] <- switch.ploidy.get.ploidy(qc[id, pasteu(run.name,"battenberg_purity")],
                                                                             qc[id, pasteu(run.name,"battenberg_ploidy")])

    # add pass or fail for basic filters
    qc[id, pasteu(run.name, "basic_filters_passed")] <- ifelse(qc[id, pasteu(run.name, "n_chrmissing")] != 0 |
                                                                 qc[id, pasteu(run.name, "n_chrsizeincorrect")] != 0 |
                                                                 qc[id, pasteu(run.name, "lsegs_homodellargest")] > thres.homodel.homodellargest |
                                                                 qc[id, pasteu(run.name, "nclonalpeaks")] == 0 |
                                                                 qc[id, pasteu(run.name, "nsuperclonalpeaks")] != 0 |
                                                                 #is.na(qc[id, pasteu(run.name, "vafpeaks_score")]) | # flagging this as unassessable instead of failing it through insufficient mutations.
                                                                 qc[id, pasteu(run.name, "vafpeaks_score")] > thres.purity.diff,"FAIL","PASS")

    # select reestimated purity/ploidy (these are only used if sample fails with filters)
    # if sample's ploidy is deemed wrong, make these the new parameters
    # if not, use vafpeaks params, unless it doesn't exist, in which case use DPClust params
    # if it's the last run (4), use DPClust params, unless the sample has switched ploidy twice, in which case use VAFPeaks                                        
                                               
    # if ploidy type fails, run with new ploidy params
    if (qc[id, pasteu(run.name, "ploidy_type_accepted")]=="FAIL") {
      qc[id, pasteu(run.name, "reestimated_purity")] <- qc[id, pasteu(run.name, "newploidy_purity")]
      qc[id, pasteu(run.name, "reestimated_ploidy")] <- qc[id, pasteu(run.name, "newploidy_ploidy")]
    }
    # if ploidy type passes, and we have vafpeaks params, run with those (if we don't have vafpeaks here, use dpclust)
    if (qc[id, pasteu(run.name, "ploidy_type_accepted")]=="PASS" &
        # (qc[id, pasteu(run.name, "vafpeaks_score")] > thres.purity.diff &
        !is.na(qc[id, pasteu(run.name, "vafpeaks_score")])) {
      qc[id, pasteu(run.name, "reestimated_purity")] <- qc[id, pasteu(run.name, "vafpeaks_purity")]
      qc[id, pasteu(run.name, "reestimated_ploidy")] <- qc[id, pasteu(run.name, "vafpeaks_ploidy")]
    } else if (qc[id, pasteu(run.name, "ploidy_type_accepted")]=="PASS" &
               is.na(qc[id, pasteu(run.name, "vafpeaks_score")])) {
      qc[id, pasteu(run.name, "reestimated_purity")] <- qc[id, pasteu(run.name, "dpclust_purity")]
      qc[id, pasteu(run.name, "reestimated_ploidy")] <- qc[id, pasteu(run.name, "dpclust_ploidy")]
    }
    # if ploidy passes, and we're on the last run, overwrite vafpeaks params with dpclust params
    # if ploidy has failed, it reruns with new ploidy as above (this could flipflop)
    if (qc[id, pasteu(run.name, "ploidy_type_accepted")]=="PASS" & run.name=="run3") {
      qc[id, pasteu(run.name, "reestimated_purity")] <- qc[id, pasteu(run.name, "dpclust_purity")]
      qc[id, pasteu(run.name, "reestimated_ploidy")] <- qc[id, pasteu(run.name, "dpclust_ploidy")]
    }
                                               
    # if homdels are too large, use dpclust parameters not vafpeaks parameters
    if (qc[id, pasteu(run.name, "lsegs_homodelall")] >= 100000000) {
      qc[id, pasteu(run.name, "reestimated_purity")] <- qc[id, pasteu(run.name, "dpclust_purity")]
      qc[id, pasteu(run.name, "reestimated_ploidy")] <- qc[id, pasteu(run.name, "dpclust_ploidy")]
    }

  }

  log.message("applying filters")
  filters <- list()
  # all chrs must be listed at least once in the subclones file
  filters$chrmissing <- ifelse(qc[ids, pasteu(run.name, "n_chrmissing")] != 0, "FAIL", "PASS")
  # chr size must pass the required threshold listed in the qc confi file
  filters$chrsizewrong <- ifelse(qc[ids, pasteu(run.name, "n_chrsizeincorrect")] != 0, "FAIL", "PASS")
  # no single homdel can be larger than the threshold listed in the qc config file
  filters$homodeletions <- ifelse(qc[ids, pasteu(run.name, "lsegs_homodellargest")] > thres.homodel.homodellargest, "FAIL", "PASS")

  # there must exist a clonal peak between the boundaries in the qc config file, dependent on narrow or wide status of sample
  filters$noclonalpeak <- ifelse(qc[ids, pasteu(run.name, "nclonalpeaks")] == 0, "FAIL", "PASS")
  # there must be no superclonal peaks that are of a size larger than the threshold in the qc config file
  filters$superclonalpeaks <- ifelse(qc[ids, pasteu(run.name, "nsuperclonalpeaks")] != 0, "FAIL", "PASS")

  # vafpeaks filter must have passed or be FLAG in order for the sample to pass (ie if not enough mutations to assess, we only use other metrics to assess the call)
  filters$vafpeaks <- ifelse(peaks$peaks_analysis$QC=="FAIL", "FAIL", "PASS")                                            

  # the ploidy call must be deemed correct, so if diploid, must have a ploidy classified as diploid through PCAWG eqn, and
  # can have no more than qc.config threshold of the genome between 0.5 boundary (if narrow sample (same boundary), or
  # can have a cluster of mutations between narrow 0.5 boundary if narrow sample (or wide boundary if wide sample)
  # but can't have both
  # and if tetraploid, must have a ploidy classified as tetraploid through PCAWG eqn (psi_t>LOH*-2+2.9)=tetra, else dip), and
  # must have at least 10% of the genome at 1+0, 2+1, or 3+2 states, and
  # must have peaks of mutations in 2:2 regions with multiplicity 1 and 2
  filters$ploidytype <- ifelse(qc[ids, pasteu(run.name,"ploidy_type_accepted")] == "PASS", "PASS", "FAIL")
  # rest of the filters must be deemed ok
  filters$basicfilterspassed <- ifelse(filters$chrmissing=="PASS" &
                                         filters$chrsizewrong=="PASS" &
                                         filters$homodeletions=="PASS" &
                                         filters$noclonalpeak=="PASS" &
                                         filters$superclonalpeaks=="PASS" &
                                         (filters$vafpeaks=="PASS" | filters$vafpeaks=="FLAG"),"PASS","FAIL")

  # classify the sample as having passed or failed
  filters <- data.frame(filters, stringsAsFactors=F)
  filters$overallfilter <- apply(filters, 1, get.worst.filter)
  filters.qc <- filters
  dimnames(filters.qc) <- list(ids, pasteu(run.name, "filter", colnames(filters.qc)))
  qc <- data.frame(qc, filters.qc[rownames(qc), ])

  log.message("writing table of QC results")
  write.table(qc, file=filename.qc, quote=F, sep="\t", row.names=F, col.names=T)

}
