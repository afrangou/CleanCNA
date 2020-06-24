convert.impute.input.to.beagle.input=function(imputeinput,
                                              chrom,
                                              filepath,
                                              vcfversion = "4.2",
                                              genomereference = "GRCh38") {
  chrom <- if (chrom == "23")
    "X"
  else chrom
  inp <- as.data.frame(data.table::fread(imputeinput))
  coln <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", paste0("SAMP001"))
  vcf <- cbind(rep(chrom, nrow(inp)), inp[, 3], rep(".", nrow(inp)),
               inp[, 4], inp[, 5], rep(".", nrow(inp)), rep("PASS",
                                                            nrow(inp)), rep(".", nrow(inp)), rep("GT", nrow(inp)),
               paste(inp[, 6], inp[, 7], inp[, 8], sep = "-"))
  vcf[vcf[, 10] == "1-0-0", 10] <- "0/0"
  vcf[vcf[, 10] == "0-1-0", 10] <- "0/1"
  vcf[vcf[, 10] == "0-0-1", 10] <- "1/1"
  vcf <- vcf[vcf[, 10] != "0-0-0", ]
  colnames(vcf) <- coln

  cat(paste0("##fileformat=VCFv", vcfversion, "\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##reference=",
             genomereference, "\n"), file = filepath)

  suppressWarnings(write.table(vcf, file = filepath, sep = "\t",
                               col.names = T, row.names = F, quote = F, append = T))

}
