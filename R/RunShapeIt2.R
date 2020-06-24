#---------------------------#
# Function to run ShapeIt2
#---------------------------#

# Run as array job across all chrs
# sampledir = "/home/AFrangou/testicular_new/battenberg_results/tumoX_normY"
# tumourplatekey = "tumoX"
# imputeinfofile = "/home/AFrangou/ALL_100G_phase1integrated_v3_impute/impute_info_shapeit.txt"
# is.male = T
# chrom = 23
# seed = as.integer(Sys.time())
# filename.shapeit2 = "/home/AFrangou/bin/shapeit"
# filename.input = paste0(sampledir,"/E-RunShapeIt/",tumourplatekey,"_shapeit_input_chr",chrom)
# filename.output = paste0(sampledir,"/E-RunShapeIt/",tumourplatekey,"_shapeit_output_chr",chrom)



RunShapeIt2 <- function (sampledir,tumourplatekey,imputeinfofile,
                         is.male, chrom,
                         filename.input, filename.output, filename.shapeit2,
                         seed = as.integer(Sys.time())) {

        # basefile from impute defining where hap/legend/etc is per chromosome
        impute.info = read.table(imputeinfofile, stringsAsFactors = F)
        colnames(impute.info) = c("chrom", "impute_legend", "genetic_map",
                                  "impute_hap", "impute_sample", "start", "end", "is_par")
        # if male, we use X,Y (and not the second X haplotype)
        if (is.male) {
          impute.info = impute.info[impute.info$is_par == 1, ]
        }
        # get chr names
        chr_names = unique(impute.info$chrom)
        # select rows of this table for required chromosome (usually 1 row, if X/Y, then is 2 rows, either XX or XY)
        if (!is.na(chrom)) {
          impute.info = impute.info[impute.info$chrom == chr_names[chrom],]
        }

        # n in loop below is either 1 for autosomes, 2 for male X (use only two diploid regions of X), 3 for female X (use entire X cos all is diploid)
        # this corresponds to the three rows seen in the impute.info file for the X chromosome
        n = nrow(impute.info)

        for (r in 1:n) {
                # if dealing with XX/XY (XY only includes PAR1 and PAR2)
                if (chrom == 23) {
                        # create filename for saving
                        filename.input.ok <- paste0(filename.input, ".filt",r)
                        # legend input file
                        legend <- read.delim(impute.info$impute_legend[r],sep = " ", header = T, stringsAsFactors = F)
                        # make gen file using the output from PLATEKEY_impute_input_chrA.txt
                        impute <- read.table(paste0(sampledir,"/D-GenerateImputeInputFromAlleleFrequencies/",
                                                    tumourplatekey,"_impute_input_chr23.txt"),
                                             stringsAsFactors = F)
                        # filter to just positions in relevant region of X chr
                        impute <- impute[impute[[3]] %in% legend$position, ]
                        indels.inds <- which(nchar(impute[[4]])>1 | nchar(impute[[5]])>1 | impute[[4]]=="-" | impute[[5]]=="-")
                        indels <- impute[indels.inds,]
                        gen <- impute[-indels.inds,]
                        write.table(gen, file = paste0(filename.input.ok, ".gen"),quote = F,sep = " ",row.names = F,col.names = F)
                        write.table(indels, file = paste0(filename.input.ok, ".gen.missing"),quote = F,sep = " ",row.names = F,col.names = F)

                        # change filename for .sample file for X (done in PrepShapeIt2 in bbdpcrecall)
                        #system(paste0("cp ", filename.input, ".sample ",
                        #              filename.input.ok, ".sample"))
                } else { # for all autosomes, just make the gen files
                        filename.input.ok <- filename.input
                        impute <- read.table(paste0(sampledir,"/D-GenerateImputeInputFromAlleleFrequencies/",
                                                    tumourplatekey,"_impute_input_chr",chrom,".txt"),
                                             stringsAsFactors = F)
                        # remove indels
                        indels.inds <- which(nchar(impute[[4]])>1 | nchar(impute[[5]])>1 | impute[[4]]=="-" | impute[[5]]=="-")
                        indels <- impute[indels.inds,]
                        gen <- impute[-indels.inds,]
                        write.table(gen, file = paste0(filename.input.ok, ".gen"),quote = F,sep = " ",row.names = F,col.names = F)
                        write.table(indels, file = paste0(filename.input.ok, ".gen.missing"),quote = F,sep = " ",row.names = F,col.names = F)
                }
                # create shapeit command
                cmd = paste(filename.shapeit2,
                            "-G", filename.input.ok,
                            "-M", impute.info$genetic_map[r],
                            "-R", impute.info$impute_hap[r],
                            impute.info$impute_legend[r],
                            impute.info$impute_sample[r],
                            "-O", ifelse(n == 1, filename.output, paste(filename.output,r, sep = ".")), ifelse(!impute.info$is_par[r],
                                                                                   "--chrX", ""),
                            "--thread 1",
                            "--seed", seed,
                            "--force")
                system(cmd, wait = T)
        }
        if (n > 1) {
                # n is >1 if we are dealing with the X chr
                # in this case, load all the haps files and rbind together to make it into one file
                filenames.haps <- paste(filename.output, 1:n, "haps",sep = ".")
                haps <- do.call("rbind", sapply(filenames.haps, read.delim,
                                                sep = " ",
                                                header = F,
                                                stringsAsFactors = F,
                                                simplify = F))
                haps <- haps[!duplicated(haps[[3]]), ]
                haps <- haps[order(haps[[3]]), ]
                write.table(haps, file = paste0(filename.output, ".haps"),
                            quote = F, sep = " ", row.names = F, col.names = F)
                system(paste0("cp ", filename.output, ".1.sample ", filename.output,
                              ".sample"))
        }
}

# RunShapeIt2(sampledir,tumourplatekey,imputeinfofile,
#             is.male, chrom,
#             filename.input, filename.output, filename.shapeit2,
#             seed = as.integer(Sys.time()))
