# compute peaks in copy-number-state-specific VAF distributions
# written by Anna Frangou using code from Alex J. Cornish & Giulio Caravagna
#
# usage:
# Rscript computePeaksVAF.R [FILENANE_SMALL_VARIANT_VCF] [FILENAME_BATTENBERG_SEGS] [FILENAME_BATTENBERG_PURITY_PLOIDY] [FILENAME_CONFIG] [FILENAME_OUTPUT]
# Requires R/3.4.0 (or maybe any 3.4.X)


ComputePeaksVAF_nf <- function(filename.ssm,
                            filename.segs,
                            filename.purity.ploidy,
                            filename.config,
                            filename.output,
                            output_dir) {
  # required libraries
  # library(peakPick)
  # library(dplyr)
  # library(stringr)
  # library(VariantAnnotation)

  # source configuration
  source(filename.config)
  message("- sourced config file")

  # read in required data
  message("- reading in data ",filename.ssm)
  # message(filename.ssm)
  # chrs = 1:22
  # message(chrs)
  snvs <- CleanCNA:::read.snvs(filename.ssm, chrs,nucleotides=c("A", "C", "G", "T"))
  message(paste(filename.ssm," - vcf read"))
  segs <- CleanCNA:::read.segs(filename.segs, chrs)
  message(paste(filename.segs," - subclones read"))
  purity.ploidy <- read.delim(filename.purity.ploidy, sep="\t", header=T, stringsAsFactors=F)
  message(paste(filename.purity.ploidy," - purity read"))
  purity <- purity.ploidy[1, "cellularity"]

  # setup output
  peaks <- list(states=list(), summary=list())

  # compute peaks
  message("- computing copy-number-state-specific VAF distribution peaks")
  for (state in intersect(states, segs$state)) {
    vafs <- snvs$vaf[queryHits(findOverlaps(snvs, segs[segs$state == state, ]))] # get VAF of SNVs in regions of this copy number state
    n <- length(vafs) # number of SNVs in region of this copy number state
    if (n > 1) {
      # only complete peak calling if segments contains >1 SNVs, otherwise error produced
      peaks$states[[state]] <- list()
      peaks$states[[state]]$vafs <- vafs
      peaks$states[[state]]$n <- n

      # add expected location of peaks
      copy.numbers <- as.numeric(unlist(str_split(state, ":")))
      possible.multiplicities <- 1:max(copy.numbers)
      peaks$states[[state]]$expected.locs <- data.frame(
        state=state,
        ploidy=sum(as.numeric(unlist(str_split(state, ":")))),
        multiplicity=possible.multiplicities,
        peak.expected=sapply(possible.multiplicities, computed.expected.peak.loc, purity, sum(copy.numbers)),
        peak.matched=NA,
        peak.diff=NA,
        purity.expected=purity,
        purity.reestimated=NA,
        purity.diff=NA
      )

      # identify peaks in segments
      peaks$states[[state]]$den <- density(peaks$states[[state]]$vafs, kernel="gaussian", n=2048) # run KDE
      peaks$states[[state]]$peaks <- CleanCNA:::call.peaks(peaks$states[[state]]$den)
      peaks.avail <- peaks$states[[state]]$peaks[peaks$states[[state]]$peaks$y > thres.peaks.y, ] # consider peaks with y value greater than threshold available for matching

      for (i in nrow(peaks$states[[state]]$expected.locs):1) {
        # if clonal peak, match to peak with highest VAF, otherwise match to peak with closest VAF
        multiplicity <- peaks$states[[state]]$expected.locs$multiplicity[i]
        match.clonal <- state %in% c("1:0", "1:1") | multiplicity == 2

        match.successful <- F
        while (!match.successful & nrow(peaks.avail)) {
          # continue until match is successful or we run out of peaks

          # whether we are trying to find the clonal peak or not
          if (match.clonal) {
            best.peak <- which(peaks.avail$x == max(peaks.avail$x)) # select peak with highest VAF
          } else {
            peaks.avail$diff <- abs(peaks.avail$x - peaks$states[[state]]$expected.locs$peak.expected[i])
            best.peak <- which(peaks.avail$diff == min(peaks.avail$diff)) # select peak with closest location
          }

          # reestimate purity using best-matching peak, and claim success if the reestimated purity is less than 1
          peak.matched <- peaks.avail$x[best.peak]
          purity.reestimated <- CleanCNA:::reestimate.purity(peak.matched, multiplicity, peaks$states[[state]]$expected.locs$ploidy[i])
          match.successful <- purity.reestimated < 1

          # remove as available peak if the match was unsuccessful (previously also removed successful matched if other peaks were available, but now stopped doing this. should look out to see whether this reinstroduced the balanced error issue)
          if (!match.successful) peaks.avail <- peaks.avail[-best.peak, ]
        }

        # if match was successful, add information to table of expected locations
        if (match.successful) {
          peaks$states[[state]]$expected.locs$peak.matched[i] <- peak.matched
          peaks$states[[state]]$expected.locs$purity.reestimated[i] <- purity.reestimated
        }
      }

      # compute difference in expected/observed peak position and purity
      peaks$states[[state]]$expected.locs$peak.diff <- peaks$states[[state]]$expected.locs$peak.matched - peaks$states[[state]]$expected.locs$peak.expected
      peaks$states[[state]]$expected.locs$purity.diff <- peaks$states[[state]]$expected.locs$purity.reestimated - peaks$states[[state]]$expected.locs$purity.expected
    }
  }

  # identify copy-number-states to consider when re-estimating sample purity
  n.mutations <- sapply(peaks$states, function(i) i$n)

  # combine all expected peak locations, and remove those expected peaks that were not matched
  if (length(n.mutations) > 0 & !all(n.mutations == 0)) {
    states.consider <- names(n.mutations)[n.mutations / sum(n.mutations) > thres.peak.prop.n] # identify copy-number-states containing greater number of mutations
    peaks$summary$expected.locs <- do.call("rbind", lapply(states.consider, function(state) peaks$states[[state]]$expected.locs)) # combine expected peaks from copy-number states we are considering
    peaks$summary$expected.locs <- peaks$summary$expected.locs[!is.na(peaks$summary$expected.locs$purity.diff), ] # remove any expected peaks that were not matched
    if (nrow(peaks$summary$expected.locs) == 0) peaks$summary$expected.locs <- NA # if no expected peaks were matched, set to NA
  } else {
    peaks$summary$expected.locs <- NA
  }

  # re-estimate purity if possible
  peaks$summary$purity.old <- purity
  if (is.null(nrow(peaks$summary$expected.locs))) {
    # if there are no matched peaks then it is not possible to reestimate purity, so set to NA
    for (element in c("purity.new", "eta")) peaks$summary[[element]] <- NA
  } else {
    # if there is more than 1 matched peak, remove the most poorly matched peak
    if (nrow(peaks$summary$expected.locs) > 1) peaks$summary$expected.locs <- peaks$summary$expected.locs[order(abs(peaks$summary$expected.locs$purity.diff), decreasing=T), ][-1, ]

    # compute weights for each matched peak
    q <- table(peaks$summary$expected.locs$state)[peaks$summary$expected.locs$state]
    N <- sum(n.mutations[unique(as.character(peaks$summary$expected.locs$state))])
    peaks$summary$expected.locs$w <- n.mutations[as.character(peaks$summary$expected.locs$state)] / (N * q)
    if (abs(sum(peaks$summary$expected.locs$w) - 1) > 1E-6) stop("expected peak weights do not sum to 1")

    # re-estimate purity and compute difference with old purity estimate
    peaks$summary$purity.new <- sum(peaks$summary$expected.locs$w * peaks$summary$expected.locs$purity.reestimated)
    peaks$summary$eta <- abs(peaks$summary$purity.new - peaks$summary$purity.old)
  }

  # save output
  message("- saving peak data to ", filename.output)
  save(peaks, file=filename.output)

  print(peaks$summary)

  # write a PASS or FAIL file saying whether this call passes or fails the peak calling assessment
  if (is.na(peaks$summary$eta)) {
    write.csv("FAILPEAKS",paste0(output_dir,"/FAILPEAKS"),quote=F,row.names = F)
  } else if (peaks$summary$eta > 0.05) {
    write.csv("FAILPEAKS",paste0(output_dir,"/FAILPEAKS"),quote=F,row.names = F)
  } else if (peaks$summary$eta <= 0.05) {
    write.csv("PASSPEAKS",paste0(output_dir,"/PASSPEAKS"),quote=F,row.names = F)
  }

  # print(paste("Written indicator peak calling pass/fail file for ",sampledir))

}
