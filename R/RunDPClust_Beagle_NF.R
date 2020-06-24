#--------------------------#
# Function to run DPClust
# Beagle, NF
#--------------------------#

RunDPClust_Beagle_nf <- function(tumourplatekey,
                       normalplatekey,
                       rhoandpsifilepath,
                       dpinfofilepath,
                       output_folder) {

  # define nucleotides
  nucleotides = c("A","C","G","T")
  # samplename
  samplename=paste0("tumo",tumourplatekey,"_norm",normalplatekey)
  # print sample running
  print(paste0(tumourplatekey," ",normalplatekey," preparing"))

  if(file.exists(rhoandpsifilepath) & file.exists(dpinfofilepath)){

    cellularity = read.table(rhoandpsifilepath,header=T,stringsAsFactors=F,sep="\t")[2,1]

    DP.info = read.table(dpinfofilepath,sep="\t",header=T, stringsAsFactors=F)

    #filter out VAF=0 and mutations without copy number info
    DP.info = DP.info[!is.na(DP.info$subclonal.fraction) & DP.info$subclonal.fraction>0,]

    mutCount = array(DP.info$mut.count,c(nrow(DP.info),1))
    WTCount = array(DP.info$WT.count,c(nrow(DP.info),1))
    totalCopyNumber = array(DP.info$subclonal.CN,c(nrow(DP.info),1))
    copyNumberAdjustment = array(DP.info$no.chrs.bearing.mut,c(nrow(DP.info),1))
    mutation.copy.number = array(DP.info$mutation.copy.number,c(nrow(DP.info),1))

    no.iters=1000

    output_folder=output_folder

    print(paste("Using cellularity file",rhoandpsifilepath))
    print(paste("Using DPinput file",dpinfofilepath))
    print(paste("Writing output to",output_folder))

    #This function does everything. It may be better to run separate functions to perform Gibbs sampling
    # and mutation assignment
    DirichletProcessClustering(mutCount = mutCount,
                               WTCount = WTCount,
                               totalCopyNumber = totalCopyNumber,
                               copyNumberAdjustment = copyNumberAdjustment,
                               mutation.copy.number = mutation.copy.number,
                               cellularity = cellularity,
                               output_folder = output_folder,
                               no.iters = no.iters,
                               no.iters.burn.in = ceiling(no.iters/5),
                               subsamplesrun = samplename,
                               samplename=tumourplatekey,
                               conc_param = 1,
                               cluster_conc = 5,
                               mut.assignment.type = 1,
                               most.similar.mut = NA,
                               mutationTypes = NA,
                               min.frac.snvs.cluster = NA,
                               max.considered.clusters=30)

    print(paste0("Clustering has run for ",tumourplatekey,"_",normalplatekey))

  } else {

    print(paste0(rhoandpsifilepath," or ",dpinfofilepath," not found, exiting."))
  }

}

