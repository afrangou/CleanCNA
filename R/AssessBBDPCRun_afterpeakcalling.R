#--------------------------------------------------------#
# Collect metrics for a run of Battenberg and DPClust
#--------------------------------------------------------#

# Metric collection done per sample
# Probably good to summarise all runs so can quickly identify those which still have a few catch up
# runs to do. This will mean loads of columns though. Hm.

# DPC info collected per sample is:
# And this is done for every run of DPC performed so we can see the improvement in the calls
# number mutations
# number mutations after 1pc removal
# number clusters
# number clusters after 1pc removal
# position clonal cluster (clonal = cluster closest to clonal)
# number mutations clonal cluster
# position top CCF cluster
# number mutations top CCF cluster
# position cluster most mutations
# number mutations cluster most mutations
# position superclonal cluster
# number mutations superclonal cluster
# position biggest subclonal cluster
# number mutations biggest subclonal cluster
# is there a cluster at 50% in dpc
# position of 50% cluster
# number mutations of 50% cluster
# passed dpc

# BB info collected is:
# purity
# ploidy
# genome length
## number of segments with a clonal aberration
## length of genome with a clonal aberration
## number of segments with a subclonal aberration
## length of genome with a subclonal aberration
## number of segments with a clonal LOH
## length of genome with a clonal LOH
## number of segments with a subclonal LOH
## length of genome with a subclonal LOH

## number of segments clonally 1+1
# genome clonally 1+1
##number of segments subclonally 1+1
# genome subclonally 1+1
##number of segments clonally 2+2
# genome clonally 2+2
##number of segments subclonally 2+2
# genome subclonally 2+2

# number of segments that is clonal where either major or minor copy number is odd
# length of genome that has frac1A=1, ie clonal, where either the major or the minor copy number is odd

# number of segments at 50pc (between 0.95/2 ± 0.05)
# length of genome at 50pc

# number mutations on 1 chrom copy
# number mutations on 2 chrom copy
# number mutations on >1 chrom copy

# number of segments clonal homdels
# total length bp clonal homdels
# longest clonal homdel bp

# number of segments subclonal homdels
# total length bp subclonal homdels
# longest subclonal homdel bp

# number of mutations in regions that are clonally homdel

# sampledir = paste0(resultsdir,"tumo",data[i,3],"_norm",data[i,4])
# participantid = data[i,2]
# tumourplatekey = data[i,3]
# normalplatekey = data[i,4]
# run = 4
# lowerclonal = 0.9
# upperclonal = 1.1

AssessBBDPCRun_afterpeakcalling <- function(sampledir,participantid,tumourplatekey,normalplatekey,run,lowerclonal,upperclonal) {

  tumournormal = paste0('tumo',tumourplatekey,'_norm',normalplatekey)
  sampletumnorm = paste0(participantid,"_",tumournormal)

  print(participantid)
  print(tumournormal)
  print(paste("DPClust run",run))

  # get dpc run optima info
  if (run==1) {
    dpcdir = "/P-DPC1/DPClust/"
    dpcprepdir = "/P-DPC1/DPPrep/"
    subdir = "/O-Postprocessing/"
    puritydir = "/M-CallSubclones/"
    assessmentdir = "/Q-AssessBB1DPC1/"
    dir.create(paste0(sampledir,"/Q-AssessBB1DPC1/"))
  } else if (run==2) {
    dpcdir = "/S-DPC2/DPClust/"
    dpcprepdir = "/S-DPC2/DPPrep/"
    subdir = "/R-BB2/O-Postprocessing/"
    puritydir = "/R-BB2/M-CallSubclones/"
    assessmentdir = "/T-AssessBB2DPC2/"
    dir.create(paste0(sampledir,"/T-AssessBB2DPC2/"))
  } else if (run==3) {
    dpcdir = "/V-DPC3/DPClust/"
    dpcprepdir = "/V-DPC3/DPPrep/"
    subdir = "/U-BB3/O-Postprocessing/"
    puritydir = "/U-BB3/M-CallSubclones/"
    assessmentdir = "/W-AssessBB3DPC3/"
    dir.create(paste0(sampledir,"/W-AssessBB3DPC3/"))
  } else if (run==4) {
    dpcdir = "/Y-DPC4/DPClust/"
    dpcprepdir = "/Y-DPC4/DPPrep/"
    subdir = "/X-BB4/O-Postprocessing/"
    puritydir = "/X-BB4/M-CallSubclones/"
    assessmentdir = "/Z-AssessBB4DPC4/"
    dir.create(paste0(sampledir,"/Z-AssessBB4DPC4/"))
  }

  dpcfilename=paste0(sampledir,dpcdir,tumournormal,"_optimaInfo.txt")

  if (file.exists(dpcfilename)) {

    optima = read.table(dpcfilename,header=T,stringsAsFactors = F)

    # some output
    print(paste("Loaded dpclust optima info file",tumourplatekey))

    # collect dpc metrics across all runs that were done for this sample
    dpcmetrics = matrix(NA,nrow=0,ncol=18)

    # beginning & total number of mutations
    number.mutations=sum(as.numeric(optima[,3]))
    # beginning & total number of clusters
    number.clusters=nrow(optima)

    # remove any cluster with <-1% mutations in
    toremove=which(optima[,3]<=number.mutations/100)

    if(length(toremove)>0) {optima=optima[-toremove,]}

    # total mutations after removal of smaller clusters
    number.mutations.post1pcremoval=sum(as.numeric(optima[,3]))
    number.clusters.post1pcremoval=nrow(optima)

    # position clonal cluster
    dists=abs(optima[,2]-1) #(upperclonal+lowerclonal)/2)
    position.clonalcluster=optima[which(dists==min(dists)),2]

    # number mutations clonal cluster
    numbermutations.clonalcluster=optima[which(dists==min(dists)),3]

    # position top ccf cluster
    position.topccfcluster=optima[nrow(optima),2]

    # number mutations top ccf cluster
    numbermutations.topccfcluster=optima[nrow(optima),3]

    # position cluster most mutations
    # a v small number of samples with v few mutations have clusters with the same number
    # of mutations in them. choose the top CCF cluster in this situation.
    most.muts=which(optima[,3]==max(optima[,3]))
    if (length(most.muts)>1) {
      position.clustermostmutations=optima[most.muts[length(most.muts)],2]
    } else {
      position.clustermostmutations=optima[which(optima[,3]==max(optima[,3])),2]
    }
    # number mutations cluster with most mutations
    if (length(most.muts)>1) {
      numbermutations.clustermostmutations=optima[most.muts[length(most.muts)],3]
    } else {
      numbermutations.clustermostmutations=optima[which(optima[,3]==max(optima[,3])),3]
    }
    # position super clonal cluster
    # number mutations superclonal cluster
    superclonal=which(optima[,2]>upperclonal)

    if(length(superclonal)>0) {
      position.superclonalcluster=optima[max(superclonal),2]
      numbermutations.superclonalcluster=optima[max(superclonal),3]
    } else {
      position.superclonalcluster=0
      numbermutations.superclonalcluster=0
    }

    # position biggest subclonal cluster
    # number mutations biggest subclonal cluster
    subclonal=which(optima[,2]<lowerclonal)
    if (length(subclonal)==0) {
      position.biggestsubclonalcluster=0
      numbermutations.biggestsubclonalcluster=0
    } else if (length(subclonal)==1) {
      position.biggestsubclonalcluster=optima[which(optima[subclonal,3]==max(optima[subclonal,3])),2]
      numbermutations.biggestsubclonalcluster=optima[which(optima[subclonal,3]==max(optima[subclonal,3])),3]
    }  else if (length(subclonal)>1) {
      subclonal.clusters=which(optima[subclonal,3]==max(optima[subclonal,3]))
      position.biggestsubclonalcluster=optima[subclonal.clusters[length(subclonal.clusters)],2]
      numbermutations.biggestsubclonalcluster=optima[subclonal.clusters[length(subclonal.clusters)],3]
    }

    # ask if the cluster closest to clonal has the largest number of mutations or has the highest CCF.
    # if ((optima[which(dists==min(dists)),3]==max(optima[,3]) | optima[which(dists==min(dists)),3]==optima[nrow(optima),3]) & (position.clonalcluster>=lowerclonal & position.clonalcluster<=upperclonal)) {
    #   passeddpc="Yes"
    # } else {
    #   passeddpc="No"
    # }

    # final test for samples which have gone through a fourth round of copy number calling
    # and been run with input parameters from peak calling
    # for the sample to pass, we ask for there to exist a peak between 0.9 and 1.1 and for no superclonal peak after removal of <1% clusters to exist above 1.1
    if (length(which(optima[,2]>=lowerclonal & optima[,2]<=upperclonal))>0) { #  & length(which(optima[,2]>upperclonal))==0
      passeddpc="Yes"
    } else {
      passeddpc="No"
    }

    # information about number of chromosomes each mutation is present on
    # comes from dpinput file
    dpcprepfilename = paste0(sampledir,dpcprepdir,tumournormal,"_DPinput.txt")
    dpcprep = read.table(dpcprepfilename,header=T,stringsAsFactors = F)

    # number of mutations present on one chr
    no.muts.chr1 = length(which(dpcprep[,16]==1))
    # number of mutations present on two chrs
    no.muts.chr2 = length(which(dpcprep[,16]==2))
    # number of mutations present on more than one chrs
    no.muts.chrmore1 =  length(which(dpcprep[,16]>1))

    # is there a cluster around 0.5 of the clonal cluster
    upperbound = (position.clonalcluster/2)+0.05
    lowerbound = (position.clonalcluster/2)-0.05
    if (length(which(optima$location<upperbound & optima$location>lowerbound))==0) {
      cluster.at.50 = FALSE
      position.cluster.at.50 = FALSE
      number.muts.cluster.at.50 = FALSE
    } else {
      cluster.at.50 = TRUE
      position.cluster.at.50 = optima[which(optima$location<upperbound & optima$location>lowerbound),2]
      number.muts.cluster.at.50 = optima[which(optima$location<upperbound & optima$location>lowerbound),3]
    }

    dpcmetrics = cbind(number.mutations,number.mutations.post1pcremoval,number.clusters,number.clusters.post1pcremoval,
                       position.clonalcluster,numbermutations.clonalcluster,position.topccfcluster,numbermutations.topccfcluster,
                       position.clustermostmutations,numbermutations.clustermostmutations,position.superclonalcluster,
                       numbermutations.superclonalcluster,position.biggestsubclonalcluster,numbermutations.biggestsubclonalcluster,
                       no.muts.chr1,no.muts.chr2,no.muts.chrmore1,cluster.at.50,position.cluster.at.50,number.muts.cluster.at.50,
                       passeddpc)

    print(paste("Collected DPClust metrics"))

  } else {

    dpcmetricscolnames = c("number.mutations","number.mutations.post1pcremoval","number.clusters","number.clusters.post1pcremoval",
                           "position.clonalcluster","numbermutations.clonalcluster","position.topccfcluster","numbermutations.topccfcluster",
                           "position.clustermostmutations","numbermutations.clustermostmutations","position.superclonalcluster",
                           "numbermutations.superclonalcluster","position.biggestsubclonalcluster","numbermutations.biggestsubclonalcluster",
                           "no.muts.chr1","no.muts.chr2","no.muts.chrmore1","cluster.at.50","position.cluster.at.50","number.muts.cluster.at.50",
                           "passeddpc")
    dpcmetrics =matrix(NA,1,21)
    colnames(dpcmetrics)=dpcmetricscolnames

    print(paste("No DPClust runs completed"))

  }

  # load subclones and purity files from BB
  subfilename = paste0(sampledir,subdir,tumourplatekey,"_subclones_hg38.txt")
  purityfilename = paste0(sampledir,puritydir,tumourplatekey,"_cellularity_ploidy.txt")

  if (file.exists(subfilename)) {

    sub = read.table(subfilename,header=T,stringsAsFactors = F)
    # some output
    print(paste("Loaded subclones hg38 file",tumourplatekey))

    purityploidy = read.table(purityfilename,header=T,stringsAsFactors = F)
    # some output
    print(paste("Loaded cellularity file",tumourplatekey))

    #  purity and ploidy
    purity = purityploidy[1,1]
    ploidy = purityploidy[1,3]

    # number of segments
    segments = nrow(sub)

    # genome length
    genome.length = as.numeric(sub[,3])-as.numeric(sub[,2])
    genome.length = sum(as.numeric(genome.length))

    # number of segments clonally diploid
    clonallydiploid = sub[which(sub[,8]==1 & sub[,9]==1 & sub[,10]==1),]
    segments.clonallydiploid = nrow(clonallydiploid)
    # genome length clonally diploid
    length.clonallydiploid = sum(as.numeric(clonallydiploid[,3])-as.numeric(clonallydiploid[,2]))

    # number of segments subclonally diploid
    subclonallydiploid = sub[which((sub[,8]==1 & sub[,9]==1 & sub[,10]!=1) | sub[,11]==1 & sub[,12]==1),]
    segments.subclonallydiploid = nrow(subclonallydiploid)
    # genome length subclonally diploid
    length.subclonallydiploid = sum(as.numeric(subclonallydiploid[,3])-as.numeric(subclonallydiploid[,2]))

    # number of segments clonally tetraploid
    clonallytetraploid = sub[which(sub[,8]==2 & sub[,9]==2 & sub[,10]==1),]
    segments.clonallytetraploid = nrow(clonallytetraploid)
    # genome length clonally tetraploid
    length.clonallytetraploid = sum(as.numeric(clonallytetraploid[,3])-as.numeric(clonallytetraploid[,2]))

    # number of segments subclonally tetraploid
    subclonallytetraploid = sub[which((sub[,8]==2 & sub[,9]==2 & sub[,10]!=1) | sub[,11]==2 & sub[,12]==2),]
    segments.subclonallytetraploid = nrow(subclonallytetraploid)
    # genome length subclonally tetraploid
    length.subclonallytetraploid = sum(as.numeric(subclonallytetraploid[,3])-as.numeric(subclonallytetraploid[,2]))

    # number of segments with a clonal aberration
    clonalaberration = sub[which(sub$frac1_A==1 & (sub$nMaj1_A!=1 | sub$nMin1_A!=1)),]
    segments.clonalaberration = nrow(clonalaberration)
    # length of genome with a clonal aberration
    length.clonalaberration = sum(as.numeric(clonalaberration[,3]-clonalaberration[,2]))

    # number of segments with a subclonal aberration (this just means any subclonality)
    subclonalaberration = sub[which(sub$frac1_A!=1 & !is.na(sub$frac1_A)),]
    segments.subclonalaberration = nrow(subclonalaberration)
    # length of genome with a subclonal aberration
    length.subclonalaberration = sum(as.numeric(subclonalaberration[,3]-subclonalaberration[,2]))

    # number of segments with a clonal LOH
    clonalLOH = sub[which(sub$frac1_A==1 & (sub$nMaj1_A==0 | sub$nMin1_A==0)),]
    segments.clonalLOH = nrow(clonalLOH)
    # length of genome with a clonal LOH
    length.clonalLOH = sum(as.numeric(clonalLOH[,3]-clonalLOH[,2]))

    # number of segments with a subclonal LOH
    subclonalLOH = sub[which(sub$frac1_A!=1 & (sub$nMaj1_A==0 | sub$nMin1_A==0 | sub$nMaj2_A==0 | sub$nMin2_A==0)),]
    segments.subclonalLOH = nrow(subclonalLOH)
    # length of genome with a subclonal LOH
    length.subclonalLOH = sum(as.numeric(subclonalLOH[,3]-subclonalLOH[,2]))

    # number of clonal segments where either major or minor copy number is odd
    clonalodd = sub[which(sub$frac1_A==1 & (sub$nMaj1_A%%2!=0 | sub$nMin1_A%%2!=0)),]
    segments.clonalodd = nrow(clonalodd)
    # length of clonal segments where either major or minor copy number is odd
    length.clonalodd = sum(as.numeric(clonalodd[,3]-clonalodd[,2]))

    # number of segments with 2 50% subclones (we are using 0.95 as our 1 for this dataset)
    midpoint = (upperclonal+lowerclonal)/2/2
    fiftylower = midpoint-0.025
    fiftyupper = midpoint+0.025
    fiftypcsubclones = sub[which(sub$frac1_A<fiftyupper & sub$frac1_A>fiftylower),]
    segments.fiftypcsubclones = nrow(fiftypcsubclones)
    # length of segments with 2 50% subclones (we are using 0.95 as our 1 for this dataset)
    length.fiftypcsubclones = sum(as.numeric(fiftypcsubclones[,3]-fiftypcsubclones[,2]))

    # number of segments that are clonally homdel
    clonalhomdels = sub[which(sub$frac1_A==1 & sub$nMaj1_A==0 & sub$nMin1_A==0),]
    segments.clonalhomdels = nrow(clonalhomdels)
    # length clonal homdels
    length.clonalhomdels = sum(as.numeric(clonalhomdels[,3]-clonalhomdels[,2]))
    # max length clonal homdel
    if(segments.clonalhomdels>0){
      maxlength.clonalhomdels = max(clonalhomdels[,3]-clonalhomdels[,2])
    } else {
      maxlength.clonalhomdels =0
    }

    # number of segments that are subclonally homdel
    subclonalhomdels = sub[which(sub$frac1_A!=1 & ((sub$nMaj1_A==0 & sub$nMin1_A==0) | (sub$nMaj2_A==0 & sub$nMin2_A==0))),]
    segments.subclonalhomdels = nrow(subclonalhomdels)
    # total length of subclonal homdels
    length.subclonalhomdels = sum(as.numeric(subclonalhomdels[,3]-subclonalhomdels[,2]))
    # max length subclonal homdel
    if(segments.subclonalhomdels>0){
      maxlength.subclonalhomdels = max(subclonalhomdels[,3]-subclonalhomdels[,2])
    } else {
      maxlength.subclonalhomdels = 0
    }

    bbmetrics = cbind(purity,ploidy,segments,genome.length,segments.clonallydiploid,length.clonallydiploid,
                      segments.subclonallydiploid,length.subclonallydiploid,
                      segments.clonallytetraploid,length.clonallytetraploid,
                      segments.subclonallytetraploid,length.subclonallytetraploid,
                      segments.clonalaberration,length.clonalaberration,
                      segments.subclonalaberration,length.subclonalaberration,
                      segments.clonalLOH,length.clonalLOH,
                      segments.subclonalLOH,length.subclonalLOH,
                      segments.clonalodd,length.clonalodd,
                      segments.fiftypcsubclones,length.fiftypcsubclones,
                      segments.clonalhomdels,length.clonalhomdels,maxlength.clonalhomdels,
                      segments.subclonalhomdels,length.subclonalhomdels,maxlength.subclonalhomdels)

  } else {

    bbmetricscolnames = c("purity","ploidy","segments","genome.length","segments.clonallydiploid","length.clonallydiploid",
                          "segments.subclonallydiploid","length.subclonallydiploid",
                          "segments.clonallytetraploid","length.clonallytetraploid",
                          "segments.subclonallytetraploid","length.subclonallytetraploid",
                          "segments.clonalaberration","length.clonalaberration",
                          "segments.subclonalaberration","length.subclonalaberration",
                          "segments.clonalLOH","length.clonalLOH",
                          "segments.subclonalLOH","length.subclonalLOH",
                          "segments.clonalodd","length.clonalodd",
                          "segments.fiftypcsubclones","length.fiftypcsubclones",
                          "segments.clonalhomdels","length.clonalhomdels","maxlength.clonalhomdels",
                          "segments.subclonalhomdels","length.subclonalhomdels","maxlength.subclonalhomdels")
    bbmetrics =matrix(NA,1,30)
    colnames(bbmetrics)=bbmetricscolnames
    print(paste("No BB runs completed"))

  }

  # combine dpc and bb metrics for this sample
  samplemetrics = as.data.frame(cbind(participantid,tumourplatekey,normalplatekey,bbmetrics,dpcmetrics))

  # write metrics table for this call
  write.csv(samplemetrics,paste0(sampledir,assessmentdir,tumourplatekey,"_metrics_run",run,"_",lowerclonal,"_",upperclonal),
            row.names = F,
            quote = F)

  # write a PASS or FAIL file dependent on DPClust output
  if (samplemetrics$passeddpc=="Yes") {
    write.csv("PASS",paste0(sampledir,assessmentdir,"PASS_",lowerclonal,"_",upperclonal),quote=F,row.names = F)
  } else if (samplemetrics$passeddpc=="No") {
    write.csv("FAIL",paste0(sampledir,assessmentdir,"FAIL_",lowerclonal,"_",upperclonal),quote=F,row.names = F)
  }

  print(paste("Written final metrics and indicator pass/fail file for",tumourplatekey))

  return(samplemetrics)

}




