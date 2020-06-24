writebeagle.as.impute=function(vcf, outfile) {
  beagleout <- as.data.frame(data.table::fread(vcf))
  haplotypes <- strsplit(beagleout[, 10], split = "\\|")
  dt <- cbind(paste0("snp_index", 1:nrow(beagleout)),
              paste0("rs_index",1:nrow(beagleout)),
              beagleout[, 2],
              beagleout[, 4],
              beagleout[,5],
              sapply(haplotypes, "[", 1), sapply(haplotypes, "[",2))
  write.table(dt, file = outfile,
              quote = F,
              col.names = F,
              row.names = F,
              sep = "\t")
}
