read.segs <- function(filename.ssm, chrs, states, pval.subclonal=0.05) {
  # read Battenberg segments from Battenberg output file
  # return GRanges object containing segments with states in X:X format as meta data
  segs <- read.delim(filename.ssm, sep="\t", header=T, stringsAsFactors=F)
  segs$chr <- paste0("chr", segs$chr) # ensure chromosome notation matches
  segs <- segs[segs$chr %in% chrs, ] # keep only segments on chromosomes on interest
  segs <- segs[segs$pval > pval.subclonal, ] # remove segments with evidence of subclonality
  segs <- segs[!is.na(segs$nMaj1_A) & !is.na(segs$nMin1_A), ] # remove segments where either the major or minor copy number is NA
  segs$state <- paste(segs$nMaj1_A, segs$nMin1_A, sep=":")
  GRanges(
    seqnames=segs$chr,
    ranges=IRanges(start=segs$startpos, end=segs$endpos),
    state=segs$state
  )
}
