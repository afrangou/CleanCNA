#----------------------------------------------------#
# Create WGD script and directories
#----------------------------------------------------#
#
# data=read.csv("/re_gecip/cancer_colorectal/analysisResults/2.chromosomeAberrations/I.Battenberg/archive/data.csv",
#               stringsAsFactors = F)
# data=data[,2:ncol(data)]
#
# resultsdir=data$resultsdir[i]
# participantid=data$participantid[i]
# tumourplatekey=data$tumourplatekey[i]
# normalplatekey=data$normalplatekey[i]
# tumourbam=data$tumourbam[i]
# normalbam=data$normalbam[i]
# gender=data$gender[i]
# project.code="re_gecip_cancer_colorectal"
# bb1done=T
# ccubeexists=NA
# vcfdir=NA
# ccubedir = NA
# ccubepurity = data$ccubepurity[i]
# vcffilepath = data$vcffilepath[i]
# ploidy_old = data$ploidy[i]
# purity_old = data$purity[i]
# shouldbediploid = data$shouldbedip[i]
# shouldbetetraploid = data$shouldbetetra[i]

CreateWGDrun <- function(resultsdir,
                         participantid,
                         tumourplatekey,
                         normalplatekey,
                         tumourbam,
                         normalbam,
                         gender,
                         project.code,
                         bb1done,
                         ccubeexists,
                         vcfdir,
                         ccubedir,
                         ccubepurity,
                         vcffilepath,
                         ploidy_old,
                         purity_old,
                         shouldbediploid,
                         shouldbetetraploid) {

  # sample directory
  sampledir = paste0(resultsdir,"tumo",tumourplatekey,"_norm",normalplatekey)

  # gender
  if (gender=="Female") {gendertrue = "FALSE"}
  if (gender=="Male") {gendertrue = "TRUE"}

  # create dirs and script for WGD run if a ploidy is suspected to be wrong (so Setup, FCN, CS, Convert, DPPrep, DPClust, Assess)
  dir.create(paste0(sampledir,"/ZZ-BB_WGD"))
  dir.create(paste0(sampledir,"/ZZ-BB_WGD/L-FitCopyNumber"))
  dir.create(paste0(sampledir,"/ZZ-BB_WGD/M-CallSubclones"))
  dir.create(paste0(sampledir,"/ZZ-BB_WGD/O-Postprocessing"))
  dir.create(paste0(sampledir,"/ZZ-DPC_WGD"))
  dir.create(paste0(sampledir,"/ZZ-DPC_WGD/DPPrep"))
  dir.create(paste0(sampledir,"/ZZ-DPC_WGD/DPClust"))
  dir.create(paste0(sampledir,"/ZZ-AssessBBDPC_WGD"))

  # create dirs and script for WGD run if a ploidy is suspected to be wrong (so Setup, FCN, CS, Convert, DPPrep, DPClust, Assess)
  header=matrix(nrow=26,ncol=1)
  header[1,1]="#!/usr/bin/env bash"
  header[2,1]="CONFIG=$1"
  header[3,1]=""
  header[4,1]=""
  header[5,1]=""
  # header[4,1]="# Set up run WGD parameters"
  # header[5,1]=paste0('bsub ',
  #                    ' -J"SetUpRunWGD_',tumourplatekey,'"',
  #                    # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
  #                    ' -q short -P ',project.code,
  #                    ' -o ',sampledir,"/logs/SetUpRunWGD.%J.out",
  #                    ' -e ',sampledir,"/logs/SetUpRunWGD.%J.err",
  #                    ' /home/AFrangou/battenberg_lsf_pipeline/steps/SetUpRunWGD.sh ',
  #                    sampledir,'/',tumourplatekey,'_configfile.txt')
  header[6,1]=""
  header[7,1]="# Run BB FitCopyNumber WGD"
  header[8,1]=paste0('bsub -J"FitCopyNumberWGD_',tumourplatekey,'"',
                     # ' -R"select[mem>15900] rusage[mem=15900]" -M15900',
                     ' -q short -P ',project.code,
                     ' -o ',sampledir,"/logs/FitCopyNumber_WGD.%J.out",
                     ' -e ',sampledir,"/logs/FitCopyNumber_WGD.%J.err",
                     ' /home/AFrangou/battenberg_lsf_pipeline/steps/FitCopyNumber_WGD.sh ',
                     sampledir,'/',tumourplatekey,'_configfile.txt')
  header[9,1]=""
  header[10,1]="# Run BB CallSubclones WGD"
  header[11,1]=paste0('bsub -w"done(FitCopyNumberWGD_',tumourplatekey,')"',
                      ' -J"CallSubclones_WGD_',tumourplatekey,'"',
                      # ' -R"select[mem>15900] rusage[mem=15900]" -M15900',
                      ' -q short -P ',project.code,
                      ' -o ',sampledir,"/logs/CallSubclones_WGD.%J.out",
                      ' -e ',sampledir,"/logs/CallSubclones_WGD.%J.err",
                      ' /home/AFrangou/battenberg_lsf_pipeline/steps/CallSubclones_WGD.sh ',
                      sampledir,'/',tumourplatekey,'_configfile.txt')
  header[12,1]=""
  header[13,1]="# Convert subclones from CallSubclones WGD"
  header[14,1]=paste0('bsub -w"done(CallSubclones_WGD_',tumourplatekey,')"',
                      ' -J"ConvertSubclonesWGD_',tumourplatekey,'"',
                      # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                      ' -q short -P ',project.code,
                      ' -o ',sampledir,"/logs/ConvertSubclonesWGD.%J.out",
                      ' -e ',sampledir,"/logs/ConvertSubclonesWGD.%J.err",
                      ' /home/AFrangou/battenberg_lsf_pipeline/steps/ConvertSubclonesWGD.sh ',
                      sampledir,'/',tumourplatekey,'_configfile.txt')
  header[15,1]=""
  header[16,1]="# Run DPPrep WGD"
  header[17,1]=paste0('bsub -w"done(ConvertSubclonesWGD_',tumourplatekey,')"',
                      ' -J"DPPrepWGD_',tumourplatekey,'"',
                      # ' -R"select[mem>8000] rusage[mem=8000]" -M8000',
                      ' -q short -P ',project.code,
                      ' -o ',sampledir,"/logs/DPPrepWGD.%J.out",
                      ' -e ',sampledir,"/logs/DPPrepWGD.%J.err",
                      ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPPrepWGD.sh ',
                      sampledir,'/',tumourplatekey,'_configfile.txt')
  header[18,1]=""
  header[19,1]="# Run DPClust WGD"
  header[20,1]=paste0('bsub -w"done(DPPrepWGD_',tumourplatekey,')"',
                      ' -J"DPClustWGD_',tumourplatekey,'"',
                      # ' -R"select[mem>16000] rusage[mem=16000]" -M16000',
                      ' -q short -P ',project.code,
                      ' -o ',sampledir,"/logs/DPClustWGD.%J.out",
                      ' -e ',sampledir,"/logs/DPClustWGD.%J.err",
                      ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPClustWGD.sh ',
                      sampledir,'/',tumourplatekey,'_configfile.txt')
  header[21,1]=""
  header[22,1]="# Assess PASS/FAIL of DPClust WGD"
  header[23,1]=paste0('bsub -w"done(DPClustWGD_',tumourplatekey,')"',
                      ' -J"AssessBBDPCWGD_',tumourplatekey,'"',
                      # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                      ' -q short -P ',project.code,
                      ' -o ',sampledir,"/logs/AssessBBDPCWGD.%J.out",
                      ' -e ',sampledir,"/logs/AssessBBDPCWGD.%J.err",
                      ' /home/AFrangou/battenberg_lsf_pipeline/steps/AssessBBDPCWGD.sh ',
                      sampledir,'/',tumourplatekey,'_configfile.txt')
  header[24,1]=""
  header[25,1]="# Assess PASS/FAIL status based on outcome of peak calling"
  header[26,1]=paste0('bsub -w"done(AssessBBDPCWGD_',tumourplatekey,')"',
                      ' -J"AssessBBDPC_PeaksWGD_',tumourplatekey,'"',
                      # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                      ' -q short -P ',project.code,
                      ' -o ',sampledir,"/logs/AssessBBDPC_PeaksWGD.%J.out",
                      ' -e ',sampledir,"/logs/AssessBBDPC_PeaksWGD.%J.err",
                      ' /home/AFrangou/battenberg_lsf_pipeline/steps/ComputePeaksVAFWGD.sh ',
                      sampledir,'/',tumourplatekey,'_configfile.txt')

  write.table(header,paste0(sampledir,'/RunWGD.sh'),quote=F,col.names=F,row.names=F)
  system(paste0('chmod 750 ',sampledir,'/RunWGD.sh'))

  # take new PLOIDY from old ploidy (either halved or doubled), calculate newpurity
  if(!is.na(shouldbetetraploid)) {
    # from diploid to tetraploid
    newpurity = purity_old/(2-purity_old)
    newploidy = ploidy_old*2
  } else if (!is.na(shouldbediploid)) {
    # from tetraploid to diploid
    newpurity = (2*purity_old)/(purity_old+1)
    newploidy = ploidy_old/2
  }


  # newploidy = ((purity_old * ploidy_old) + 2*(newpurity - purity_old))/newpurity

  print(paste("New purity is",newpurity))
  print(paste("New ploidy is",newploidy))

  system(paste0("echo RHO_WGD=",newpurity,
                " >> ",
                sampledir,"/",
                tumourplatekey,"_configfile.txt"))
  print("Added RHO_WGD to config file")
  system(paste0("echo PSI_WGD=",newploidy,
                " >> ",
                sampledir,"/",
                tumourplatekey,"_configfile.txt"))
  print("Added PSI_WGD to config file")

}



