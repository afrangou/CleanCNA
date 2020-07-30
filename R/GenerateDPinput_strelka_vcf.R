#------------------------------------------------------------------------------------------------#
# Function to generate the DPinput file from Battenberg output and sample vcfs which is then used in DPClust
# Beagle, NF
#------------------------------------------------------------------------------------------------#

#' @name
#' GenerateDPinput_strelka_vcf
#'
#' @title
#' Generate DP input for strelka vcfs, using tier 1 snvs
#'
#' @description
#' Converts snv counts from vcf into input format for DPClust
#'
#' @param
#' Loci of SNP positions used in CNV call
#' Allele frequencies of same SNP positions
#' File containing purity and ploidy of tumour sample
#' Subclones file from Battenberg containing copy number segments across genome
#' Sex of patient
#'
#' @return
#' A file containing all snv loci,
#' their mutation copy number status,
#' and their CCF (cancer cell fraction)

GenerateDPinput_strelka_vcf <- function(tumourplatekey,
                                        normalplatekey,
                                        gender,
                                        vcffilepath,
                                        rhoandpsifilepath,
                                        subcloneshg38filepath,
                                        output_loci_file,
                                        output_file,
                                        output_DPinput_file) {

  # source(config_file)

  # definitions
  nucleotides = c("A","C","G","T")

  # print samplename to logfile
  samplename=paste0("tumo",tumourplatekey,"_norm",normalplatekey)

  # give some info
  print(paste("Sample name:",samplename))
  #print(paste("This is run",run))
  print(paste("Using",subcloneshg38filepath))

  # check vcf exists and turn it into format for DPClust, else print that it doesn't exist and the sample is exiting
  if (file.exists(vcffilepath)) {

    # some info
    print(paste("VCF is ",vcffilepath))

    # load vcf
    snvs=read.table(vcffilepath,stringsAsFactors=F)
    print(paste("vcfloaded",vcffilepath))

    # select only those variants which have passed all filters and which are SNVs (not indels)
    snvs=snvs[which(snvs[,7]=="PASS" & nchar(snvs[,4])==1 & nchar(snvs[,5])==1),]

    # make new snvs table for dpinput
    newsnvs=snvs[,c(1,2,2,4,5,8,10,11)]
    newsnvs[,2]=newsnvs[,3]-1

    # remove stuff not assigned to a chr
    if ("chrX" %in% snvs[,1]) {
      newsnvs=newsnvs[1:max(which(newsnvs[,1]=="chrX")),]
    }
    # replace chr1 with 1 etc
    newsnvs[,1]=gsub("chr","",newsnvs[,1])

    # turn counts from vcf into required format
    tumA=sapply(sapply(newsnvs[,7],function(x){strsplit(x,":")[[1]][5]}),function(x){strsplit(x,",")[[1]][1]})
    tumC=sapply(sapply(newsnvs[,7],function(x){strsplit(x,":")[[1]][6]}),function(x){strsplit(x,",")[[1]][1]})
    tumG=sapply(sapply(newsnvs[,7],function(x){strsplit(x,":")[[1]][7]}),function(x){strsplit(x,",")[[1]][1]})
    tumT=sapply(sapply(newsnvs[,7],function(x){strsplit(x,":")[[1]][8]}),function(x){strsplit(x,",")[[1]][1]})
    # create file type for tumour
    tumpos=cbind(newsnvs[,c(1,3,4,5),],0,0,0,0,0)
    colnames(tumpos)=c("Chromosome","Position","REF","ALT","count_A","count_C","count_G","count_T","total_depth")
    tumpos[,5]=tumA
    tumpos[,6]=tumC
    tumpos[,7]=tumG
    tumpos[,8]=tumT
    refcol = tumpos$REF
    refcol[which(refcol=="A")]=1
    refcol[which(refcol=="C")]=2
    refcol[which(refcol=="G")]=3
    refcol[which(refcol=="T")]=4
    refcol=as.integer(refcol)+4
    altcol = tumpos$ALT
    altcol[which(altcol=="A")]=1
    altcol[which(altcol=="C")]=2
    altcol[which(altcol=="G")]=3
    altcol[which(altcol=="T")]=4
    altcol=as.integer(altcol)+4
    refcounts=tumpos[cbind(seq_along(refcol),refcol)]
    altcounts=tumpos[cbind(seq_along(altcol),altcol)]
    tumpos[,9]=as.integer(refcounts)+as.integer(altcounts)

    # file containing all SNV positions and counts of all alleles at these positions
    # normalpositioncounts=normpos[,c(1,2,5:9)]
    tumourpositioncounts=tumpos[,c(1,2,5:9)]

    print(paste("Writing tumour position counts to ",
                output_file))
    write.table(tumourpositioncounts,
                output_file,
                row.names=F,
                quote=F,
                sep="\t")
    # file containing just the loci for input into the function below
    tumourloci=tumpos[,c(1:4)]
    print(paste("Writing tumour loci to ",
                output_loci_file))
    write.table(tumourloci,
                output_loci_file,
                row.names=F,
                col.names=F,
                quote=F,
                sep="\t")

    # define filenames for function
    loci_file=output_loci_file
    allele_frequencies_file=output_file
    cellularity_file=rhoandpsifilepath
    subclone_file=subcloneshg38filepath
    output_DPinput_file=output_DPinput_file

    print(paste("Using loci file",loci_file))
    print(paste("Using allele frequencies file",allele_frequencies_file))
    print(paste("Using cellularity file",cellularity_file))
    print(paste("Using subclones file",subclone_file))
    print(paste("Writing to output file",output_DPinput_file))


    if (gender=="FEMALE") {
      gender="female"
    } else if (gender=="MALE") {
      gender="male"
    } else {
      print("Gender unspecified, exiting")
      break
    }

    runGetDirichletProcessInfo(loci_file,
                               allele_frequencies_file,
                               cellularity_file,
                               subclone_file,
                               gender,
                               SNP.phase.file="NA",
                               mut.phase.file="NA",
                               output_file=output_DPinput_file)


  } else {

    print(paste("No vcf available - exiting"))

  }
}

