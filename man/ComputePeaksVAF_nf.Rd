\name{ComputePeaksVAF_nf}
\alias{ComputePeaksVAF_nf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Peaks VAF (for nf)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ComputePeaksVAF_nf(filename.ssm, filename.segs, filename.purity.ploidy, filename.config, filename.output, output_dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename.ssm}{
%%     ~~Describe \code{filename.ssm} here~~
}
  \item{filename.segs}{
%%     ~~Describe \code{filename.segs} here~~
}
  \item{filename.purity.ploidy}{
%%     ~~Describe \code{filename.purity.ploidy} here~~
}
  \item{filename.config}{
%%     ~~Describe \code{filename.config} here~~
}
  \item{filename.output}{
%%     ~~Describe \code{filename.output} here~~
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (filename.ssm, filename.segs, filename.purity.ploidy,
    filename.config, filename.output, output_dir)
{
    source(filename.config)
    message("- sourced config file")
    message("- reading in data ", filename.ssm)
    snvs <- read.snvs(filename.ssm, chrs, nucleotides = c("A",
        "C", "G", "T"))
    message(paste(filename.ssm, " - vcf read"))
    segs <- read.segs(filename.segs, chrs)
    message(paste(filename.segs, " - subclones read"))
    purity.ploidy <- read.delim(filename.purity.ploidy, sep = "\t",
        header = T, stringsAsFactors = F)
    message(paste(filename.purity.ploidy, " - purity read"))
    purity <- purity.ploidy[1, "cellularity"]
    peaks <- list(states = list(), summary = list())
    message("- computing copy-number-state-specific VAF distribution peaks")
    for (state in intersect(states, segs$state)) {
        vafs <- snvs$vaf[queryHits(findOverlaps(snvs, segs[segs$state ==
            state, ]))]
        n <- length(vafs)
        if (n > 1) {
            peaks$states[[state]] <- list()
            peaks$states[[state]]$vafs <- vafs
            peaks$states[[state]]$n <- n
            copy.numbers <- as.numeric(unlist(str_split(state,
                ":")))
            possible.multiplicities <- 1:max(copy.numbers)
            peaks$states[[state]]$expected.locs <- data.frame(state = state,
                ploidy = sum(as.numeric(unlist(str_split(state,
                  ":")))), multiplicity = possible.multiplicities,
                peak.expected = sapply(possible.multiplicities,
                  computed.expected.peak.loc, purity, sum(copy.numbers)),
                peak.matched = NA, peak.diff = NA, purity.expected = purity,
                purity.reestimated = NA, purity.diff = NA)
            peaks$states[[state]]$den <- density(peaks$states[[state]]$vafs,
                kernel = "gaussian", n = 2048)
            peaks$states[[state]]$peaks <- call.peaks(peaks$states[[state]]$den)
            peaks.avail <- peaks$states[[state]]$peaks[peaks$states[[state]]$peaks$y >
                thres.peaks.y, ]
            for (i in nrow(peaks$states[[state]]$expected.locs):1) {
                multiplicity <- peaks$states[[state]]$expected.locs$multiplicity[i]
                match.clonal <- state \%in\% c("1:0", "1:1") |
                  multiplicity == 2
                match.successful <- F
                while (!match.successful & nrow(peaks.avail)) {
                  if (match.clonal) {
                    best.peak <- which(peaks.avail$x == max(peaks.avail$x))
                  }
                  else {
                    peaks.avail$diff <- abs(peaks.avail$x - peaks$states[[state]]$expected.locs$peak.expected[i])
                    best.peak <- which(peaks.avail$diff == min(peaks.avail$diff))
                  }
                  peak.matched <- peaks.avail$x[best.peak]
                  purity.reestimated <- reestimate.purity(peak.matched,
                    multiplicity, peaks$states[[state]]$expected.locs$ploidy[i])
                  match.successful <- purity.reestimated < 1
                  if (!match.successful)
                    peaks.avail <- peaks.avail[-best.peak, ]
                }
                if (match.successful) {
                  peaks$states[[state]]$expected.locs$peak.matched[i] <- peak.matched
                  peaks$states[[state]]$expected.locs$purity.reestimated[i] <- purity.reestimated
                }
            }
            peaks$states[[state]]$expected.locs$peak.diff <- peaks$states[[state]]$expected.locs$peak.matched -
                peaks$states[[state]]$expected.locs$peak.expected
            peaks$states[[state]]$expected.locs$purity.diff <- peaks$states[[state]]$expected.locs$purity.reestimated -
                peaks$states[[state]]$expected.locs$purity.expected
        }
    }
    n.mutations <- sapply(peaks$states, function(i) i$n)
    if (length(n.mutations) > 0 & !all(n.mutations == 0)) {
        states.consider <- names(n.mutations)[n.mutations/sum(n.mutations) >
            thres.peak.prop.n]
        peaks$summary$expected.locs <- do.call("rbind", lapply(states.consider,
            function(state) peaks$states[[state]]$expected.locs))
        peaks$summary$expected.locs <- peaks$summary$expected.locs[!is.na(peaks$summary$expected.locs$purity.diff),
            ]
        if (nrow(peaks$summary$expected.locs) == 0)
            peaks$summary$expected.locs <- NA
    }
    else {
        peaks$summary$expected.locs <- NA
    }
    peaks$summary$purity.old <- purity
    if (is.null(nrow(peaks$summary$expected.locs))) {
        for (element in c("purity.new", "eta")) peaks$summary[[element]] <- NA
    }
    else {
        if (nrow(peaks$summary$expected.locs) > 1)
            peaks$summary$expected.locs <- peaks$summary$expected.locs[order(abs(peaks$summary$expected.locs$purity.diff),
                decreasing = T), ][-1, ]
        q <- table(peaks$summary$expected.locs$state)[peaks$summary$expected.locs$state]
        N <- sum(n.mutations[unique(as.character(peaks$summary$expected.locs$state))])
        peaks$summary$expected.locs$w <- n.mutations[as.character(peaks$summary$expected.locs$state)]/(N *
            q)
        if (abs(sum(peaks$summary$expected.locs$w) - 1) > 1e-06)
            stop("expected peak weights do not sum to 1")
        peaks$summary$purity.new <- sum(peaks$summary$expected.locs$w *
            peaks$summary$expected.locs$purity.reestimated)
        peaks$summary$eta <- abs(peaks$summary$purity.new - peaks$summary$purity.old)
    }
    message("- saving peak data to ", filename.output)
    save(peaks, file = filename.output)
    print(peaks$summary$eta)
    if (is.na(peaks$summary$eta)) {
        write.csv("FAILPEAKS", paste0(output_dir, "/FAILPEAKS"),
            quote = F, row.names = F)
    }
    else if (peaks$summary$eta > 0.05) {
        write.csv("FAILPEAKS", paste0(output_dir, "/FAILPEAKS"),
            quote = F, row.names = F)
    }
    else if (peaks$summary$eta <= 0.05) {
        write.csv("PASSPEAKS", paste0(output_dir, "/PASSPEAKS"),
            quote = F, row.names = F)
    }
    print(paste("Written indicator peak calling pass/fail file for ",
        sampledir))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
