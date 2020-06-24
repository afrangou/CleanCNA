# read vcf function
read.snvs <- function(filename, chrs, nucleotides=c("A", "C", "G", "T")) {
  message("- inside read.snvs function")
  # read SNVs from VCF containing QCed SSMs
  # return GRanges object containing VAF as meta data
  vcf <- readVcf(filename, genome="hg38")
  message("- read vcf")
  vcf <- vcf[filt(vcf) == "PASS", ] # keep only passing variants
  vcf <- vcf[as.character(seqnames(vcf)) %in% chrs, ] # keep only required chromosomes
  vcf <- vcf[as.character(ref(vcf)) %in% nucleotides & as.character(unlist(alt(vcf))) %in% nucleotides, ] # keep only SNVs
  GRanges(
    seqnames=as.character(seqnames(vcf)),
    ranges=IRanges(start=start(vcf), width=1),
    vaf=info(vcf)$Tumour_BAF
  )
}
