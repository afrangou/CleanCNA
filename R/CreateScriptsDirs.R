#----------------------------------------------------#
# Create scripts and directories for runs or reruns
#----------------------------------------------------#

# RUNS FROM THE CANCER_ANALYSIS TABLE PLUS A COUPLE OF PARAMETERS:
# WHETHER THE FIRST RUN OF BB HAS ALREADY BEEN DONE
# WHETHER WE HAVE A CCUBE PURITY CALL

# samples = read.csv("~/re_gecip/cancer_sarcoma/0.sarcomaLandscape/2.0.chromosomeAberrations/sampleList.2019-07-02_09-05-00.modified.tsv",
#                    stringsAsFactors=F,
#                    sep="\t")
#


CreateScriptsDirs <- function(resultsdir,
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
                              vcffilepath) {

        # define gender true false
        if (gender=="Female") {gendertrue = "FALSE"}
        if (gender=="Male") {gendertrue = "TRUE"}

        # Make config file for peaks run (make this more manipulable later, just using standard for now)
        system(paste0("cp ",resultsdir,"peakconfigfile.R ",
                      resultsdir,"tumo",tumourplatekey,"_norm",normalplatekey,"/",
                      tumourplatekey,"_peakconfigfile.R"))

        # if don't yet have a first BB run
        if (bb1done==F) {
                # create overall sample dir
                sampledir = paste0(resultsdir,"tumo",tumourplatekey,"_norm",normalplatekey)
                dir.create(sampledir)

                # create subdirs & subsubdirs
                dir.create(paste0(sampledir,"/A-GetAlleleCounts"))
                dir.create(paste0(sampledir,"/B-RunBAFLogR"))
                dir.create(paste0(sampledir,"/C-RunGCcorrect"))
                dir.create(paste0(sampledir,"/D-GenerateImputeInputFromAlleleFrequencies"))
                dir.create(paste0(sampledir,"/E-RunImpute"))
                dir.create(paste0(sampledir,"/F-CombineImputeOutputs"))
                dir.create(paste0(sampledir,"/G-GetHaplotypedBAFs"))
                dir.create(paste0(sampledir,"/H-CleanUp"))
                dir.create(paste0(sampledir,"/I-PlotHaplotypedData"))
                dir.create(paste0(sampledir,"/J-CombineBAFfiles"))
                dir.create(paste0(sampledir,"/K-SegmentBAFphased"))
                dir.create(paste0(sampledir,"/L-FitCopyNumber"))
                dir.create(paste0(sampledir,"/M-CallSubclones"))
                dir.create(paste0(sampledir,"/N-wrapup"))
                dir.create(paste0(sampledir,"/O-Postprocessing"))
                dir.create(paste0(sampledir,"/P-DPC1"))
                dir.create(paste0(sampledir,"/P-DPC1/DPPrep"))
                dir.create(paste0(sampledir,"/P-DPC1/DPClust"))
                dir.create(paste0(sampledir,"/Q-AssessBB1DPC1"))
                dir.create(paste0(sampledir,"/R-BB2"))
                dir.create(paste0(sampledir,"/R-BB2/L-FitCopyNumber"))
                dir.create(paste0(sampledir,"/R-BB2/M-CallSubclones"))
                dir.create(paste0(sampledir,"/R-BB2/O-Postprocessing"))
                dir.create(paste0(sampledir,"/S-DPC2"))
                dir.create(paste0(sampledir,"/S-DPC2/DPPrep"))
                dir.create(paste0(sampledir,"/S-DPC2/DPClust"))
                dir.create(paste0(sampledir,"/T-AssessBB2DPC2"))
                dir.create(paste0(sampledir,"/U-BB3"))
                dir.create(paste0(sampledir,"/U-BB3/L-FitCopyNumber"))
                dir.create(paste0(sampledir,"/U-BB3/M-CallSubclones"))
                dir.create(paste0(sampledir,"/U-BB3/O-Postprocessing"))
                dir.create(paste0(sampledir,"/V-DPC3"))
                dir.create(paste0(sampledir,"/V-DPC3/DPPrep"))
                dir.create(paste0(sampledir,"/V-DPC3/DPClust"))
                dir.create(paste0(sampledir,"/W-AssessBB3DPC3"))
                dir.create(paste0(sampledir,"/X-BB4"))
                dir.create(paste0(sampledir,"/X-BB4/L-FitCopyNumber"))
                dir.create(paste0(sampledir,"/X-BB4/M-CallSubclones"))
                dir.create(paste0(sampledir,"/X-BB4/O-Postprocessing"))
                dir.create(paste0(sampledir,"/Y-DPC4"))
                dir.create(paste0(sampledir,"/Y-DPC4/DPPrep"))
                dir.create(paste0(sampledir,"/Y-DPC4/DPClust"))
                dir.create(paste0(sampledir,"/Z-AssessBB4DPC4"))
                dir.create(paste0(sampledir,"/ZZ-BB_WGD"))
                dir.create(paste0(sampledir,"/ZZ-BB_WGD/L-FitCopyNumber"))
                dir.create(paste0(sampledir,"/ZZ-BB_WGD/M-CallSubclones"))
                dir.create(paste0(sampledir,"/ZZ-BB_WGD/O-Postprocessing"))
                dir.create(paste0(sampledir,"/ZZ-DPC_WGD"))
                dir.create(paste0(sampledir,"/ZZ-DPC_WGD/DPPrep"))
                dir.create(paste0(sampledir,"/ZZ-DPC_WGD/DPClust"))
                dir.create(paste0(sampledir,"/ZZ-AssessBBDPC_WGD"))
                dir.create(paste0(sampledir,"/logs"))

                # create basic config file (parameters will be added later if reruns are required)
                config=matrix(nrow=40,ncol=1)
                config[1,1]=paste0("RUN_DIR=",resultsdir)
                config[2,1]=paste0("LOG_DIR=",sampledir,"/logs/")
                config[3,1]=paste0("OUTPUT_DIR=",sampledir)
                config[4,1]=paste0("SAMPLE_PATH=",sampledir)
                config[5,1]="PIPELINE_DIR=/home/AFrangou/battenberg-lsf-pipeline"
                config[6,1]=paste0("TUMOURNAME=",tumourplatekey)
                config[7,1]=paste0("IS_MALE=",gendertrue)
                config[8,1]=paste0("NORMALNAME=",normalplatekey)
                config[9,1]=paste0("TUMOURBAM=",tumourbam)
                config[10,1]=paste0("NORMALBAM=",normalbam)
                config[11,1]="PLATFORM_GAMMA=1"
                config[12,1]="PHASING_GAMMA=1"
                config[13,1]="SEGMENTATION_GAMMA=10"
                config[14,1]="CLONALITY_DIST_METRIC=0"
                config[15,1]="ASCAT_DIST_METRIC=1"
                config[16,1]="MIN_PLOIDY=1.6"
                config[17,1]="MAX_PLOIDY=4.8"
                config[18,1]="MIN_RHO=0.13"
                config[19,1]="MAX_RHO=1.02"
                config[20,1]="MIN_GOODNESS_OF_FIT=0.63"
                config[21,1]="BALANCED_THRESHOLD=0.51"
                config[22,1]="IMPUTEINFOFILE=/home/AFrangou/ALL_100G_phase1integrated_v3_impute/impute_info.txt"
                config[23,1]="IMPUTE_EXE=/home/AFrangou/impute_v2.3.2_x86_64_static/impute2"
                config[24,1]=paste0("SEED=",as.integer(Sys.time()))
                config[25,1]="MAX_CN_STATE=250"
                config[26,1]="SV_BREAKPOINTS_FILE=NA"
                config[27,1]="PROBLEMLOCI=/home/AFrangou/probloci_270415.txt"
                config[28,1]="G1000_PREFIX_HG38=/home/AFrangou/battenberg_1000genomesloci2012_v3_hg38/hg38_alleleFiles_feb22_chr"
                config[29,1]="G1000_PREFIX_AC_HG38=/home/AFrangou/battenberg_1000genomesloci2012_v3_hg38/hg38_lociFiles_feb22_chr"
                config[30,1]="G1000_PREFIX=/home/AFrangou/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr"
                config[31,1]="G1000_PREFIX_AC=/home/AFrangou/runallelecounter/1000genomesloci2012/1000genomesloci2012_chr"
                config[32,1]="GCCORRECTPREFIX=/home/AFrangou/battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_filled_chr_"
                config[33,1]="MIN_NORMAL_DEPTH=10"
                config[34,1]="ALLELECOUNTER=/home/AFrangou/bin/alleleCounter"
                config[35,1]=paste0("PARTICIPANTID=",participantid)
                config[36,1]=paste0("VCFFILEPATH=",vcffilepath)
                config[37,1]=paste0("CCUBEDIR=",ccubedir)
                config[38,1]=paste0("CCUBEPURITY=",ccubepurity)
                config[39,1]="SHAPEIT_EXE=/home/AFrangou/bin/shapeit"
                config[40,1]="SHAPEITINFOFILE=/home/AFrangou/ALL_100G_phase1integrated_v3_impute/impute_info_shapeit.txt"

                write.table(config,paste0(sampledir,"/",tumourplatekey,"_configfile.txt"),
                            quote=F,
                            col.names=F,
                            row.names=F)

                # create script for first Battenberg run
                header=matrix(nrow=51,ncol=1)
                header[1,1]="#!/usr/bin/env bash"
                header[2,1]="CONFIG=$1"

                header[4,1]="# Get allele frequencies"
                header[5,1]="# Tumour"
                header[6,1]=paste0('bsub -J"GetAlleleFrequenciesTumour_',
                                  tumourplatekey,'[1-23]" -q medium -P ',project.code,
                                  ' -o ',sampledir,
                                  '/logs/GetAlleleFrequenciesTumour.%J.%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/GetAlleleCounts_tumour.sh ',
                                  sampledir,'/',tumourplatekey,'_configfile.txt')

                header[7,1]="# Normal"
                header[8,1]=paste0('bsub -J"GetAlleleFrequenciesNormal_',
                                  normalplatekey,'[1-23]" -q medium -P ',project.code,
                                  ' -o ',sampledir,
                                  '/logs/GetAlleleFrequenciesNormal.%J.%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/GetAlleleCounts_normal.sh ',
                                  sampledir,'/',tumourplatekey,'_configfile.txt')

                header[10,1]="# Remove 'chr' from alleleCounter output, array"
                header[11,1]=paste0('bsub -w"done(GetAlleleFrequenciesTumour_',
                                    tumourplatekey,'[1-23]) && done(GetAlleleFrequenciesNormal_',
                                    normalplatekey,'[1-23])" -J"RunRemoveCHR_',tumourplatekey,
                                    '" -q short -P ',project.code,
                                    ' -o ',sampledir,'/logs/RemoveCHR.%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/RemoveCHR.sh ',
                                    sampledir,'/',tumourplatekey,'_configfile.txt')

                header[13,1]=paste0("# Convert alleleCounter output back to hg37")
                header[14,1]=paste0('bsub -w"done(RunRemoveCHR_',tumourplatekey,')" -J"RunSwitchback_hg38_to_hg37_',
                                    tumourplatekey,'[1-23]" -q short -P ',project.code,
                                   ' -o ',sampledir,'/logs/Switchback_hg38_to_hg37.%J.%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/Switchback_hg38_to_hg37.sh ',
                                   sampledir,'/',tumourplatekey,'_configfile.txt')

                header[16,1]="# Get BAF and logR"
                header[17,1]=paste0('bsub -w"done(RunSwitchback_hg38_to_hg37_',
                                   tumourplatekey,'[1-23])" -R"select[mem>28000] rusage[mem=28000]" -M28000 -J"RunBAFLogR_',
                                   tumourplatekey,'" -q short -P ',project.code,
                                   ' -o ',sampledir,'/logs/RunBAFLogR.%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/RunBAFLogR.sh ',
                                   sampledir,'/',tumourplatekey,'_configfile.txt')

                header[19,1]='# RunGCCorrect_wgs'
                header[20,1]=paste0('bsub -w "done(RunBAFLogR_',tumourplatekey,
                                    ')" -R"select[mem>35000] rusage[mem=35000]" -M35000 -J"runGCcorrect_',
                                    tumourplatekey,'" -q medium -P ',project.code,
                                    ' -o ',sampledir,'/logs/runGCcorrect',tumourplatekey,
                                    '.%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/RunGCcorrect_wgs.sh ',
                                    sampledir,'/',tumourplatekey,'_configfile.txt')

                header[22,1]='# Perform phasing'
                header[23,1]='# GenerateImputeInputFromAlleleFrequencies'
                header[24,1]=paste0('bsub -w "done(RunSwitchback_hg38_to_hg37_',
                                    tumourplatekey,'[1-23])" -R"select[mem>4000] rusage[mem=4000]" -M4000 -J"GenerateImputeInputs_',
                                    tumourplatekey,'[1-23]" -q short -P ',project.code,
                                    ' -o ',sampledir,
                                    '/logs/GenerateImputeInputs.%J.%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/GenerateImputeInputFromAlleleFrequencies.sh ',
                                    sampledir,'/',tumourplatekey,'_configfile.txt')

                header[26,1]='# RunImpute'
                header[27,1]=paste0('bsub -w "done(GenerateImputeInputs_',
                                   tumourplatekey,'[*])" -R"select[mem>8000] rusage[mem=8000]" -M8000 -J"Impute_',
                                   tumourplatekey,'[1-23]" -q short -P ',project.code,
                                   ' -o ',sampledir,
                                   '/logs/Impute.%J.%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/RunImpute.sh ',
                                   sampledir,'/',tumourplatekey,'_configfile.txt')

                header[29,1]='# CombineImputeOutputs'
                header[30,1]=paste0('bsub -w"done(Impute_',tumourplatekey,
                                    '[*])" -R"select[mem>4000] rusage[mem=4000]" -M4000 -J"CombineImputeOutputs_',
                                    tumourplatekey,'[1-23]" -q short -P ',project.code,' -o ',
                                    sampledir,'/logs/CombineImputeOutputs.%J.%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/CombineImputeOutputs.sh ',
                                    sampledir,'/',tumourplatekey,'_configfile.txt')

                header[32,1]='# GetHaplotypeBAFs'
                header[33,1]=paste0('bsub -w"done(CombineImputeOutputs_',
                                   tumourplatekey,'[*])" -R"select[mem>2000] rusage[mem=2000]" -M2000 -J"GetHaplotypedBAFs_',
                                   tumourplatekey,'[1-23]" -q short -P ',project.code,
                                   ' -o ',sampledir,
                                   '/logs/GetHaplotypedBAFs.%J.%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/GetHaplotypedBAFs.sh ',
                                   sampledir,'/',tumourplatekey,'_configfile.txt')

                header[35,1]='# CleanUp'
                header[36,1]=paste0('bsub -w"done(GetHaplotypedBAFs_',
                                   tumourplatekey,'[*])" -J"CleanUp_',
                                   tumourplatekey,'[1-23]" -q short -P ',project.code,
                                   ' -o ',sampledir,
                                   '/logs/CleanUp.%J.%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/CleanUp.sh ',
                                   sampledir,'/',tumourplatekey,'_configfile.txt')

                header[38,1]='# PlotHaplotypedData'
                header[39,1]=paste0('bsub -w"done(CleanUp_',
                                    tumourplatekey,'[*])" -R"select[mem>4000] rusage[mem=4000]" -M4000 -J"PlotHaplotypedData_',
                                    tumourplatekey,'[1-23]" -q short -P ',project.code,
                                    ' -o ',sampledir,
                                    '/logs/PlotHaplotypedData.%J.%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/PlotHaplotypedData.sh ',
                                    sampledir,'/',tumourplatekey,'_configfile.txt')

                header[41,1]='# CombineBAFfiles'
                header[42,1]=paste0('bsub -w"done(PlotHaplotypedData_',
                                   tumourplatekey,')" -R"select[mem>4000] rusage[mem=4000]" -M4000 -J"CombineBAFfiles_',
                                   tumourplatekey,'" -q short -P ',project.code,
                                   ' -o ',sampledir,'/logs/CombineBAFfiles.%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/CombineBAFfiles.sh ',
                                   sampledir,'/',tumourplatekey,'_configfile.txt')

                header[44,1]='# Segmentation and copy number calling'
                header[45,1]=paste0('bsub -w"done(CombineBAFfiles_',
                                    tumourplatekey,')" -R"select[mem>4000] rusage[mem=4000]" -M4000 -J"SegmentBAFphased_',
                                    tumourplatekey,'" -q short -P ',project.code,
                                    ' -o ',sampledir,'/logs/SegmentBAFphased.%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/SegmentBAFphased.sh ',
                                    sampledir,'/',tumourplatekey,'_configfile.txt')

                header[47,1]=paste0('bsub -w"done(SegmentBAFphased_',
                                   tumourplatekey,') && done(runGCcorrect_',
                                   tumourplatekey,')" -R"select[mem>15900] rusage[mem=15900]" -M15900 -J"FitCopyNumber_',
                                   tumourplatekey,'" -q short -P ',project.code,
                                   ' -o ',sampledir,
                                   '/logs/FitCopyNumber.%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/FitCopyNumber.sh ',
                                   sampledir,'/',tumourplatekey,'_configfile.txt')

                header[49,1]=paste0('bsub -w"done(FitCopyNumber_',
                                    tumourplatekey,')" -R"select[mem>15900] rusage[mem=15900]" -M15900 -J"CallSubclones_',
                                    tumourplatekey,'" -q short -P ',project.code,
                                    ' -o ',sampledir,
                                    '/logs/CallSubclones.%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/CallSubclones.sh ',
                                    sampledir,'/',tumourplatekey,'_configfile.txt')

                header[51,1]=paste0('bsub -w"done(CallSubclones_',
                                    tumourplatekey,')" -R"select[mem>500] rusage[mem=500]" -M500 -J"wrapup_',
                                    tumourplatekey,'" -q short -P ',project.code,
                                    ' -o ',sampledir,'/logs/wrapup.%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/wrapup.sh ',
                                    sampledir,'/',tumourplatekey,'_configfile.txt')

                header[which(is.na(header[,1])),]=""
                write.table(header,paste0(sampledir,'/',tumourplatekey,"_submission.sh"),quote=F,col.names=F,row.names=F)

        } else if (bb1done==T) {

                sampledir = paste0(resultsdir,"tumo",tumourplatekey,"_norm",normalplatekey)

                # make additional directories for reruns
                dir.create(paste0(sampledir,"/O-Postprocessing"))
                dir.create(paste0(sampledir,"/P-DPC1"))
                dir.create(paste0(sampledir,"/P-DPC1/DPPrep"))
                dir.create(paste0(sampledir,"/P-DPC1/DPClust"))
                dir.create(paste0(sampledir,"/Q-AssessBB1DPC1"))
                dir.create(paste0(sampledir,"/R-BB2"))
                dir.create(paste0(sampledir,"/R-BB2/L-FitCopyNumber"))
                dir.create(paste0(sampledir,"/R-BB2/M-CallSubclones"))
                dir.create(paste0(sampledir,"/R-BB2/O-Postprocessing"))
                dir.create(paste0(sampledir,"/S-DPC2"))
                dir.create(paste0(sampledir,"/S-DPC2/DPPrep"))
                dir.create(paste0(sampledir,"/S-DPC2/DPClust"))
                dir.create(paste0(sampledir,"/T-AssessBB2DPC2"))
                dir.create(paste0(sampledir,"/U-BB3"))
                dir.create(paste0(sampledir,"/U-BB3/L-FitCopyNumber"))
                dir.create(paste0(sampledir,"/U-BB3/M-CallSubclones"))
                dir.create(paste0(sampledir,"/U-BB3/O-Postprocessing"))
                dir.create(paste0(sampledir,"/V-DPC3"))
                dir.create(paste0(sampledir,"/V-DPC3/DPPrep"))
                dir.create(paste0(sampledir,"/V-DPC3/DPClust"))
                dir.create(paste0(sampledir,"/W-AssessBB3DPC3"))
                dir.create(paste0(sampledir,"/X-BB4"))
                dir.create(paste0(sampledir,"/X-BB4/L-FitCopyNumber"))
                dir.create(paste0(sampledir,"/X-BB4/M-CallSubclones"))
                dir.create(paste0(sampledir,"/X-BB4/O-Postprocessing"))
                dir.create(paste0(sampledir,"/Y-DPC4"))
                dir.create(paste0(sampledir,"/Y-DPC4/DPPrep"))
                dir.create(paste0(sampledir,"/Y-DPC4/DPClust"))
                dir.create(paste0(sampledir,"/Z-AssessBB4DPC4"))

                # add participant ID to the config file
                system(paste0("echo PARTICIPANTID=",participantid,
                              " >> ",
                              sampledir,"/",
                              tumourplatekey,"_configfile.txt"))
                # add vcffilepath to the config file
                system(paste0("echo VCFFILEPATH=",vcffilepath,
                              " >> ",
                              sampledir,"/",
                              tumourplatekey,"_configfile.txt"))
                # add ccubefilepath to the config file
                system(paste0("echo CCUBEDIR=",ccubedir,
                              " >> ",
                              sampledir,"/",
                              tumourplatekey,"_configfile.txt"))
                # add ccube purity to config file
                system(paste0("echo CCUBEPURITY=",if(!is.na(ccubepurity)){ccubepurity/100}else{ccubepurity},
                              " >> ",
                              sampledir,"/",
                              tumourplatekey,"_configfile.txt"))

        }

        # create script for first DPClust run (does subclones conversion, DPPrep, DPClust)
        # this happens for all samples, whether first BB has run or not
        header=matrix(nrow=17,ncol=1)
        header[1,1]="#!/usr/bin/env bash"
        header[2,1]="CONFIG=$1"
        header[3.1]=""
        header[4,1]="# Convert subclones file from BB1 to hg38"
        if (bb1done==T) {
                header[5,1]=paste0('bsub -J"ConvertSubclones1_',
                                   tumourplatekey,'"',
                                   # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                                   ' -q short -P ',project.code,
                                   ' -o ',sampledir,"/logs/ConvertSubclones1.%J.out",
                                   ' -e ',sampledir,"/logs/ConvertSubclones1.%J.err",
                                   ' /home/AFrangou/battenberg_lsf_pipeline/steps/ConvertSubclones.sh ',
                                   sampledir,'/',tumourplatekey,'_configfile.txt')
        } else if (bb1done==F) {
                header[5,1]=paste0('bsub -w"done(wrapup_',tumourplatekey,')"' ,
                                   ' -J"ConvertSubclones1_',
                                   tumourplatekey,'"',
                                   # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                                   ' -q short -P ',project.code,
                                   ' -o ',sampledir,"/logs/ConvertSubclones1.%J.out",
                                   ' -e ',sampledir,"/logs/ConvertSubclones1.%J.err",
                                   ' /home/AFrangou/battenberg_lsf_pipeline/steps/ConvertSubclones.sh ',
                                   sampledir,'/',tumourplatekey,'_configfile.txt')
        }
        header[6,1]=""
        header[7,1]="# RunDPPrep1"
        header[8,1]=paste0('bsub -w"done(ConvertSubclones1_',tumourplatekey,')"',
                           ' -J"DPPrep_',tumourplatekey,'"',
                           # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                           ' -q short -P ',project.code,
                           ' -o ',sampledir,"/logs/DPPrep.%J.out",
                           ' -e ',sampledir,"/logs/DPPrep.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPPrep.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[9,1]=""
        header[10,1]="# RunDPClust1"
        header[11,1]=paste0('bsub -w"done(DPPrep_',tumourplatekey,')"',
                           ' -J"DPClust_',tumourplatekey,'"',
                           # ' -R"select[mem>16000] rusage[mem=16000]" -M16000',
                           ' -q medium -P ',project.code,
                           ' -o ',sampledir,"/logs/DPClust.%J.out",
                           ' -e ',sampledir,"/logs/DPClust.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPClust.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[12,1]=""
        header[13,1]="# Assess PASS/FAIL status based on outcome of DPClust"
        header[14,1]=paste0('bsub -w"done(DPClust_',tumourplatekey,')"',
                           ' -J"AssessBBDPC_',tumourplatekey,'"',
                           # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                           ' -q short -P ',project.code,
                           ' -o ',sampledir,"/logs/AssessBBDPC.%J.out",
                           ' -e ',sampledir,"/logs/AssessBBDPC.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/AssessBBDPC1.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[15,1]=""
        header[16,1]="# Assess PASS/FAIL status based on outcome of peak calling"
        header[17,1]=paste0('bsub -w"done(AssessBBDPC_',tumourplatekey,')"',
                            ' -J"AssessBBDPC_Peaks_',tumourplatekey,'"',
                            # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/AssessBBDPC_Peaks.%J.out",
                            ' -e ',sampledir,"/logs/AssessBBDPC_Peaks.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/ComputePeaksVAF.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')

        write.table(header,paste0(sampledir,'/Run1.sh'),quote=F,col.names=F,row.names=F)


        # create script for second run if a rerun is needed (so Setup, FCN, CS, Convert, DPPrep, DPClust, Assess)
        header=matrix(nrow=26,ncol=1)
        header[1,1]="#!/usr/bin/env bash"
        header[2,1]="CONFIG=$1"
        header[3.1]=""
        header[4,1]="# Set up run 2 (rerun 1) parameters"
        header[5,1]=paste0('bsub -w"done(AssessBBDPC_Peaks_',tumourplatekey,')" -J"SetUpRun2_',tumourplatekey,'"',
                           # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                           ' -q short -P ',project.code,
                           ' -o ',sampledir,"/logs/SetUpRun2.%J.out",
                           ' -e ',sampledir,"/logs/SetUpRun2.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/SetUpRun2.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[6,1]=""
        header[7,1]="# Run BB FitCopyNumber 2"
        header[8,1]=paste0('bsub -w"done(SetUpRun2_',tumourplatekey,')"',
                           ' -J"FitCopyNumber2_',tumourplatekey,'"',
                           # ' -R"select[mem>15900] rusage[mem=15900]" -M15900',
                           ' -q short -P ',project.code,
                           ' -o ',sampledir,"/logs/FitCopyNumber_run2.%J.out",
                           ' -e ',sampledir,"/logs/FitCopyNumber_run2.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/FitCopyNumber_run2.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[9,1]=""
        header[10,1]="# Run BB CallSubclones 2"
        header[11,1]=paste0('bsub -w"done(FitCopyNumber2_',tumourplatekey,')"',
                           ' -J"CallSubclones_run2_',tumourplatekey,'"',
                           # ' -R"select[mem>15900] rusage[mem=15900]" -M15900',
                           ' -q short -P ',project.code,
                           ' -o ',sampledir,"/logs/CallSubclones_run2.%J.out",
                           ' -e ',sampledir,"/logs/CallSubclones_run2.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/CallSubclones_run2.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[12,1]=""
        header[13,1]="# Convert subclones from CallSubclones 2"
        header[14,1]=paste0('bsub -w"done(CallSubclones_run2_',tumourplatekey,')"',
                            ' -J"ConvertSubclones2_',tumourplatekey,'"',
                            # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/ConvertSubclones2.%J.out",
                            ' -e ',sampledir,"/logs/ConvertSubclones2.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/ConvertSubclones2.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[15,1]=""
        header[16,1]="# Run DPPrep 2"
        header[17,1]=paste0('bsub -w"done(ConvertSubclones2_',tumourplatekey,')"',
                            ' -J"DPPrep2_',tumourplatekey,'"',
                            # ' -R"select[mem>8000] rusage[mem=8000]" -M8000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/DPPrep2.%J.out",
                            ' -e ',sampledir,"/logs/DPPrep2.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPPrep2.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[18,1]=""
        header[19,1]="# Run DPClust 2"
        header[20,1]=paste0('bsub -w"done(DPPrep2_',tumourplatekey,')"',
                            ' -J"DPClust2_',tumourplatekey,'"',
                            # ' -R"select[mem>16000] rusage[mem=16000]" -M16000',
                            ' -q medium -P ',project.code,
                            ' -o ',sampledir,"/logs/DPClust2.%J.out",
                            ' -e ',sampledir,"/logs/DPClust2.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPClust2.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[21,1]=""
        header[22,1]="# Assess PASS/FAIL of DPClust 2"
        header[23,1]=paste0('bsub -w"done(DPClust2_',tumourplatekey,')"',
                            ' -J"AssessBBDPC2_',tumourplatekey,'"',
                            # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/AssessBBDPC2.%J.out",
                            ' -e ',sampledir,"/logs/AssessBBDPC2.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/AssessBBDPC2.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[24,1]=""
        header[25,1]="# Assess PASS/FAIL status based on outcome of peak calling"
        header[26,1]=paste0('bsub -w"done(AssessBBDPC2_',tumourplatekey,')"',
                            ' -J"AssessBBDPC_Peaks2_',tumourplatekey,'"',
                            # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/AssessBBDPC_Peaks2.%J.out",
                            ' -e ',sampledir,"/logs/AssessBBDPC_Peaks2.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/ComputePeaksVAF2.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')

        write.table(header,paste0(sampledir,'/Run2.sh'),quote=F,col.names=F,row.names=F)

        # create script for third run if a rerun is needed (so Setup, FCN, CS, Convert, DPPrep, DPClust, Assess)
        header=matrix(nrow=26,ncol=1)
        header[1,1]="#!/usr/bin/env bash"
        header[2,1]="CONFIG=$1"
        header[3.1]=""
        header[4,1]="# Set up run 3 (rerun 2) parameters"
        header[5,1]=paste0('bsub -w"done(AssessBBDPC_Peaks2_',tumourplatekey,')"',
                           ' -J"SetUpRun3_',tumourplatekey,'"',
                           # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                           ' -q short -P ',project.code,
                           ' -o ',sampledir,"/logs/SetUpRun3.%J.out",
                           ' -e ',sampledir,"/logs/SetUpRun3.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/SetUpRun3.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[6,1]=""
        header[7,1]="# Run BB FitCopyNumber 3"
        header[8,1]=paste0('bsub -w"done(SetUpRun3_',tumourplatekey,')"',
                           ' -J"FitCopyNumber3_',tumourplatekey,'"',
                           # ' -R"select[mem>15900] rusage[mem=15900]" -M15900',
                           ' -q short -P ',project.code,
                           ' -o ',sampledir,"/logs/FitCopyNumber_run3.%J.out",
                           ' -e ',sampledir,"/logs/FitCopyNumber_run3.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/FitCopyNumber_run3.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[9,1]=""
        header[10,1]="# Run BB CallSubclones 3"
        header[11,1]=paste0('bsub -w"done(FitCopyNumber3_',tumourplatekey,')"',
                            ' -J"CallSubclones_run3_',tumourplatekey,'"',
                            # ' -R"select[mem>15900] rusage[mem=15900]" -M15900',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/CallSubclones_run3.%J.out",
                            ' -e ',sampledir,"/logs/CallSubclones_run3.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/CallSubclones_run3.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[12,1]=""
        header[13,1]="# Convert subclones from CallSubclones 3"
        header[14,1]=paste0('bsub -w"done(CallSubclones_run3_',tumourplatekey,')"',
                            ' -J"ConvertSubclones3_',tumourplatekey,'"',
                            # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/ConvertSubclones3.%J.out",
                            ' -e ',sampledir,"/logs/ConvertSubclones3.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/ConvertSubclones3.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[15,1]=""
        header[16,1]="# Run DPPrep 3"
        header[17,1]=paste0('bsub -w"done(ConvertSubclones3_',tumourplatekey,')"',
                            ' -J"DPPrep3_',tumourplatekey,'"',
                            # ' -R"select[mem>8000] rusage[mem=8000]" -M8000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/DPPrep3.%J.out",
                            ' -e ',sampledir,"/logs/DPPrep3.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPPrep3.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[18,1]=""
        header[19,1]="# Run DPClust 3"
        header[20,1]=paste0('bsub -w"done(DPPrep3_',tumourplatekey,')"',
                            ' -J"DPClust3_',tumourplatekey,'"',
                            # ' -R"select[mem>16000] rusage[mem=16000]" -M16000',
                            ' -q medium -P ',project.code,
                            ' -o ',sampledir,"/logs/DPClust3.%J.out",
                            ' -e ',sampledir,"/logs/DPClust3.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPClust3.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[21,1]=""
        header[22,1]="# Assess PASS/FAIL of DPClust 3"
        header[23,1]=paste0('bsub -w"done(DPClust3_',tumourplatekey,')"',
                            ' -J"AssessBBDPC3_',tumourplatekey,'"',
                            # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/AssessBBDPC3.%J.out",
                            ' -e ',sampledir,"/logs/AssessBBDPC3.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/AssessBBDPC3.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[24,1]=""
        header[25,1]="# Assess PASS/FAIL status based on outcome of peak calling"
        header[26,1]=paste0('bsub -w"done(AssessBBDPC3_',tumourplatekey,')"',
                            ' -J"AssessBBDPC_Peaks3_',tumourplatekey,'"',
                            # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/AssessBBDPC_Peaks3.%J.out",
                            ' -e ',sampledir,"/logs/AssessBBDPC_Peaks3.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/ComputePeaksVAF3.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')

        write.table(header,paste0(sampledir,'/Run3.sh'),quote=F,col.names=F,row.names=F)

        # create script for fourth run if a rerun is needed (so Setup, FCN, CS, Convert, DPPrep, DPClust, Assess)
        header=matrix(nrow=26,ncol=1)
        header[1,1]="#!/usr/bin/env bash"
        header[2,1]="CONFIG=$1"
        header[3.1]=""
        header[4,1]="# Set up run 4 (rerun 3) parameters"
        header[5,1]=paste0('bsub -w"done(AssessBBDPC_Peaks3_',tumourplatekey,')"',
                           ' -J"SetUpRun4_',tumourplatekey,'"',
                           # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                           ' -q short -P ',project.code,
                           ' -o ',sampledir,"/logs/SetUpRun4.%J.out",
                           ' -e ',sampledir,"/logs/SetUpRun4.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/SetUpRun4.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[6,1]=""
        header[7,1]="# Run BB FitCopyNumber 4"
        header[8,1]=paste0('bsub -w"done(SetUpRun4_',tumourplatekey,')"',
                           ' -J"FitCopyNumber4_',tumourplatekey,'"',
                           # ' -R"select[mem>15900] rusage[mem=15900]" -M15900',
                           ' -q short -P ',project.code,
                           ' -o ',sampledir,"/logs/FitCopyNumber_run4.%J.out",
                           ' -e ',sampledir,"/logs/FitCopyNumber_run4.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/FitCopyNumber_run4.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[9,1]=""
        header[10,1]="# Run BB CallSubclones 4"
        header[11,1]=paste0('bsub -w"done(FitCopyNumber4_',tumourplatekey,')"',
                            ' -J"CallSubclones_run4_',tumourplatekey,'"',
                            # ' -R"select[mem>15900] rusage[mem=15900]" -M15900',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/CallSubclones_run4.%J.out",
                            ' -e ',sampledir,"/logs/CallSubclones_run4.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/CallSubclones_run4.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[12,1]=""
        header[13,1]="# Convert subclones from CallSubclones 4"
        header[14,1]=paste0('bsub -w"done(CallSubclones_run4_',tumourplatekey,')"',
                            ' -J"ConvertSubclones4_',tumourplatekey,'"',
                            # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/ConvertSubclones4.%J.out",
                            ' -e ',sampledir,"/logs/ConvertSubclones4.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/ConvertSubclones4.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[15,1]=""
        header[16,1]="# Run DPPrep 4"
        header[17,1]=paste0('bsub -w"done(ConvertSubclones4_',tumourplatekey,')"',
                            ' -J"DPPrep4_',tumourplatekey,'"',
                            # ' -R"select[mem>8000] rusage[mem=8000]" -M8000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/DPPrep4.%J.out",
                            ' -e ',sampledir,"/logs/DPPrep4.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPPrep4.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[18,1]=""
        header[19,1]="# Run DPClust 4"
        header[20,1]=paste0('bsub -w"done(DPPrep4_',tumourplatekey,')"',
                            ' -J"DPClust4_',tumourplatekey,'"',
                            # ' -R"select[mem>16000] rusage[mem=16000]" -M16000',
                            ' -q medium -P ',project.code,
                            ' -o ',sampledir,"/logs/DPClust4.%J.out",
                            ' -e ',sampledir,"/logs/DPClust4.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPClust4.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[21,1]=""
        header[22,1]="# Assess PASS/FAIL of DPClust 4"
        header[23,1]=paste0('bsub -w"done(DPClust4_',tumourplatekey,')"',
                            ' -J"AssessBBDPC4_',tumourplatekey,'"',
                            # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/AssessBBDPC4.%J.out",
                            ' -e ',sampledir,"/logs/AssessBBDPC4.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/AssessBBDPC4.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[24,1]=""
        header[25,1]="# Assess PASS/FAIL status based on outcome of peak calling"
        header[26,1]=paste0('bsub -w"done(AssessBBDPC4_',tumourplatekey,')"',
                            ' -J"AssessBBDPC_Peaks4_',tumourplatekey,'"',
                            # ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q short -P ',project.code,
                            ' -o ',sampledir,"/logs/AssessBBDPC_Peaks4.%J.out",
                            ' -e ',sampledir,"/logs/AssessBBDPC_Peaks4.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/ComputePeaksVAF4.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')

        write.table(header,paste0(sampledir,'/Run4.sh'),quote=F,col.names=F,row.names=F)

        # create overall script to submit whole sample with all (defined) rounds of recall
        # different for those which have first BB already
        if (bb1done==F) {
                header=matrix(nrow=17,ncol=1)
                header[1,1]="#!/usr/bin/env bash"
                header[2,1]="CONFIG=$1"
                header[3.1]=""
                header[4,1]="# First Battenberg run"
                header[5,1]=paste0('bsub -J"Submission_',
                                   tumourplatekey,'"',
                                   ' -q short -P ',project.code,
                                   ' -o ',sampledir,"/logs/Submission.%J.out",
                                   ' -e ',sampledir,"/logs/Submission.%J.err ",
                                   sampledir,'/',tumourplatekey,'_submission.sh')
                header[6,1]=""
                header[7,1]="# First DPPrep and DPClust run"
                header[8,1]=paste0('bsub -w"done(Submission_',tumourplatekey,')" -J"Run1_',
                                   tumourplatekey,'"',
                                   ' -q short -P ',project.code,
                                   ' -o ',sampledir,"/logs/Run1.%J.out",
                                   ' -e ',sampledir,"/logs/Run1.%J.err ",
                                   sampledir,'/Run1.sh')
                header[9,1]=""
                header[10,1]="# Second Battenberg and DPClust run (new DPClust purity)"
                header[11,1]=paste0('bsub -w"done(Run1_',tumourplatekey,')" -J"Run2_',
                                    tumourplatekey,'"',
                                    ' -q short -P ',project.code,
                                    ' -o ',sampledir,"/logs/Run2.%J.out",
                                    ' -e ',sampledir,"/logs/Run2.%J.err ",
                                    sampledir,'/Run2.sh')
                header[12,1]=""
                header[13,1]="# Third Battenberg and DPClust run (Ccube purity)"
                header[14,1]=paste0('bsub -w"done(Run2_',tumourplatekey,')" -J"Run3_',
                                    tumourplatekey,'"',
                                    ' -q short -P ',project.code,
                                    ' -o ',sampledir,"/logs/Run3.%J.out",
                                    ' -e ',sampledir,"/logs/Run3.%J.err ",
                                    sampledir,'/Run3.sh')
                header[15,1]=""
                header[16,1]="# Fourth Battenberg and DPClust run (Peaks purity)"
                header[17,1]=paste0('bsub -w"done(Run3_',tumourplatekey,')" -J"Run4_',
                                    tumourplatekey,'"',
                                    ' -q short -P ',project.code,
                                    ' -o ',sampledir,"/logs/Run4.%J.out",
                                    ' -e ',sampledir,"/logs/Run4.%J.err ",
                                    sampledir,'/Run4.sh')
        } else if (bb1done==T) {
                  header=matrix(nrow=14,ncol=1)
                  header[1,1]="#!/usr/bin/env bash"
                  header[2,1]="CONFIG=$1"
                  header[3.1]=""

                  header[4,1]="# First DPPrep and DPClust run"
                  header[5,1]=paste0('bsub -J"Run1_',tumourplatekey,'"',
                                     ' -q short -P ',project.code,
                                     ' -o ',sampledir,"/logs/Run1.%J.out",
                                     ' -e ',sampledir,"/logs/Run1.%J.err ",
                                     sampledir,'/Run1.sh')
                  header[6,1]=""
                  header[7,1]="# Second Battenberg and DPClust run (new DPClust purity)"
                  header[8,1]=paste0('bsub -w"done(Run1_',tumourplatekey,')" -J"Run2_',
                                      tumourplatekey,'"',
                                      ' -q short -P ',project.code,
                                      ' -o ',sampledir,"/logs/Run2.%J.out",
                                      ' -e ',sampledir,"/logs/Run2.%J.err ",
                                      sampledir,'/Run2.sh')
                  header[9,1]=""
                  header[10,1]="# Third Battenberg and DPClust run (Ccube purity)"
                  header[11,1]=paste0('bsub -w"done(Run2_',tumourplatekey,')" -J"Run3_',
                                      tumourplatekey,'"',
                                      ' -q short -P ',project.code,
                                      ' -o ',sampledir,"/logs/Run3.%J.out",
                                      ' -e ',sampledir,"/logs/Run3.%J.err ",
                                      sampledir,'/Run3.sh')
                  header[12,1]=""
                  header[13,1]="# Fourth Battenberg and DPClust run (Peaks purity)"
                  header[14,1]=paste0('bsub -w"done(Run3_',tumourplatekey,')" -J"Run4_',
                                      tumourplatekey,'"',
                                      ' -q short -P ',project.code,
                                      ' -o ',sampledir,"/logs/Run4.%J.out",
                                      ' -e ',sampledir,"/logs/Run4.%J.err ",
                                      sampledir,'/Run4.sh')
        }
        write.table(header,paste0(sampledir,'/Submit_sample.sh'),quote=F,col.names=F,row.names=F)

}



# for (i in unfinishedthird) {
#
#   CreateScriptsDirs(resultsdir=data$resultsdir[i],
#                     participantid=data$participantid[i],
#                     tumourplatekey=data$tumourplatekey[i],
#                     normalplatekey=data$normalplatekey[i],
#                     tumourbam=data$tumourbam[i],
#                     normalbam=data$normalbam[i],
#                     gender=data$gender[i],
#                     project.code="re_gecip_cancer_colorectal",
#                     bb1done=T,
#                     ccubeexists=NA,
#                     vcfdir=NA,
#                     ccubedir = NA,
#                     ccubepurity = data$ccubepurity[i],
#                     vcffilepath = data$vcffilepath[i])
#
#
# }
#
# script = matrix(nrow=2500,ncol=1)
# script[1,1]="#!/usr/bin/env bash"
# script[2,1]=""
# script[3,1]=""
# for (i in unfinishedthird) {
#   script[i+3,1]=paste0("sh ","/home/AFrangou/colorectal/battenberg_results_fixvaf/tumo",
#                        data[i,3],"_norm",data[i,4],
#                        "/Submit_sample.sh")
# }
# script[which(is.na(script[,1]))]=""
# write.table(script,"~/colorectal/battenberg_results_fixvaf/Lastfew.sh",row.names=F,col.names=F,quote=F)
# system(paste0("chmod 750 ~/colorectal/battenberg_results_fixvaf/Lastfew.sh"))
#
#
# # check ccube purities, alex changed from 43 to 0.43 so my /100 made purities tiny. check which
# # cases this happened and fix
# for (i in 1:nrow(data)) {
#
#
#
#
#
# }



# # test CreateScriptsDirsFixedParams on  (after running below few lines)
# i=4
# CreateScriptsDirsFixedParams(resultsdir=data[i,1],
#                              participantid=data[i,2],
#                              tumourplatekey=data[i,3],
#                              normalplatekey=data[i,4],
#                              tumourbam=data[i,5],
#                              normalbam=data[i,6],
#                              gender=data[i,7],
#                              project.code=data[i,8],
#                              fixed_rho = 0.5,
#                              fixed_psi = 2)
#
# # colorectal
# samples=read.csv("/re_gecip/cancer_colorectal/V6_sample_sheet.csv",stringsAsFactors=F,sep="\t")
#
# resultsdir = "/home/AFrangou/colorectal/battenberg_results/"
# participantid = samples$participant_id_use
# tumourplatekey = samples$tumour_sample_platekey
# normalplatekey = samples$germline_sample_platekey
# tumourbam = samples$tumour_bam
# normalbam = samples$germline_bam
# gender = samples$sex
# project.code = "re_gecip_cancer_colorectal"
# bb1done = T
# ccubeexists = T
# vcfdir = "/re_gecip/cancer_colorectal/analysisResults/0.variantCalls/A.StrelkaV6/"
# data = cbind(resultsdir,
#              participantid,
#              tumourplatekey,
#              normalplatekey,
#              tumourbam,
#              normalbam,
#              gender,
#              project.code,
#              bb1done,
#              ccubeexists,
#              vcfdir)
# data[which(data[,7]=="FEMALE"),7]="Female"
# data[which(data[,7]=="MALE"),7]="Male"
#
# for (i in 1:10) {
#
#   CreateScriptsDirs(resultsdir=data[i,1],
#                     participantid=data[i,2],
#                     tumourplatekey=data[i,3],
#                     normalplatekey=data[i,4],
#                     tumourbam=data[i,5],
#                     normalbam=data[i,6],
#                     gender=data[i,7],
#                     project.code=data[i,8],
#                     bb1done=data[i,9],
#                     ccubeexists=data[i,10],
#                     vcf=data[i,11])
#
# }
#
#
# numbersamples=10
# script = matrix(nrow=numbersamples+3,ncol=1)
# script[1,1]="#!/usr/bin/env bash"
# script[2,1]=""
# script[3,1]=""
# for (i in 1:numbersamples) {
#   script[i+3,1]=paste0("sh ","/home/AFrangou/colorectal/battenberg_results/tumo",
#                        data[i,3],"_norm",data[i,4],
#                        "/Submit_sample.sh")
# }
# write.table(script,"/home/AFrangou/colorectal/automated_reruns_10.sh",row.names=F,col.names=F,quote=F)
#
# # colorectal v7 from beginning (ie first BB too)
# # sort out first BB runs for v7 CrC
# # get v7 samples to run BB first run on through David's account
# v6=read.csv("/re_gecip/cancer_colorectal/V6_sample_sheet.csv",stringsAsFactors=F,sep="\t") # 1678
# v7=read.csv("/re_gecip/cancer_colorectal/analysisResults/2.chromosomeAberrations/I.Battenberg/A.labkeytables/cancer_analysis_2019-07-30_11-19-10.csv",
#             stringsAsFactors = F)
# v7=v7[which(v7$Library.Type=="PCR-Free" & v7$Preparation.Method=="FF"),]
# inv6=match(v7$Tumour.Sample.Platekey,v6$tumour_sample_platekey)
# # those samples new to v7
# newv7=v7[which(is.na(inv6)),]
# # those samples already in v6
# inv6=v7[which(!is.na(inv6)),]
#
# # check none of the new v7 samples haven't already been run
# # 26 have already been run - copy these over.
# subclones=c()
# subclonesheadnode = c()
# for (i in 1:nrow(newv7)) {
#   # checkwhichran[i]=dir.exists(paste0("~/re_gecip/cancer_colorectal/analysisResults/2.chromosomeAberrations/I.Battenberg/01_Battenberg_original/",
#   #                         newv7$Participant.Id[i],"/tumo",newv7$Tumour.Sample.Platekey[i],"_norm",newv7$Germline.Sample.Platekey[i],"/"))
#   # if(checkwhichran[i]==T) {
#   subclones[i]=file.exists(paste0("/re_gecip/cancer_colorectal/analysisResults/2.chromosomeAberrations/I.Battenberg/01_Battenberg_original/",
#                                   newv7$Participant.Id[i],"/tumo",newv7$Tumour.Sample.Platekey[i],"_norm",newv7$Germline.Sample.Platekey[i],
#                                   "/M-CallSubclones/",newv7$Tumour.Sample.Platekey[i],"_subclones.txt"))
#   subclonesheadnode[i]=file.exists(paste0("/home/AFrangou/colorectal/battenberg_results/tumo",newv7$Tumour.Sample.Platekey[i],"_norm",newv7$Germline.Sample.Platekey[i],
#                                           "/M-CallSubclones/",newv7$Tumour.Sample.Platekey[i],"_subclones.txt"))
# }
#
# newv7=newv7[-which(subclonesheadnode==T),]
# write.csv(newv7,"~/colorectal/newv7samples_torun.csv",row.names=F)
#
# samples=read.csv("~/colorectal/newv7samples_torun.csv",stringsAsFactors=F)
#
# resultsdir = "/home/AFrangou/colorectal/battenberg_results/"
# participantid = samples$Participant.Id
# tumourplatekey = samples$Tumour.Sample.Platekey
# normalplatekey = samples$Germline.Sample.Platekey
# tumourbam = samples$Tumour.Bam
# normalbam = samples$Germline.Bam
# gender = samples$Sex
# project.code = "re_gecip_cancer_colorectal"
# bb1done = F
# ccubeexists = F
# vcfdir = "/re_gecip/cancer_colorectal/analysisResults/0.variantCalls/A.StrelkaV6/"
# data = cbind(resultsdir,
#              participantid,
#              tumourplatekey,
#              normalplatekey,
#              tumourbam,
#              normalbam,
#              gender,
#              project.code,
#              bb1done,
#              ccubeexists,
#              vcfdir)
# data[which(data[,7]=="FEMALE"),7]="Female"
# data[which(data[,7]=="MALE"),7]="Male"
#
# for (i in 1:nrow(data)) {
#
#   CreateScriptsDirs(resultsdir=data[i,1],
#                     participantid=data[i,2],
#                     tumourplatekey=data[i,3],
#                     normalplatekey=data[i,4],
#                     tumourbam=data[i,5],
#                     normalbam=data[i,6],
#                     gender=data[i,7],
#                     project.code=data[i,8],
#                     bb1done=data[i,9],
#                     ccubeexists=data[i,10],
#                     vcf=data[i,11])
#
# }
# numbersamples=nrow(data)
# script = matrix(nrow=numbersamples+3,ncol=1)
# script[1,1]="#!/usr/bin/env bash"
# script[2,1]=""
# script[3,1]=""
# for (i in 1:165) {
#   script[i+3,1]=paste0("sh ","/home/AFrangou/colorectal/battenberg_results/tumo",
#                        data[i,3],"_norm",data[i,4],
#                        "/",data[i,3],"_submission.sh")
# }
# write.table(script,"/home/AFrangou/colorectal/automated_firstbb_runs_v7_DCW.sh",row.names=F,col.names=F,quote=F)
# script = matrix(nrow=numbersamples+3,ncol=1)
# script[1,1]="#!/usr/bin/env bash"
# script[2,1]=""
# script[3,1]=""
# for (i in 166:329) {
#   script[i+3,1]=paste0("sh ","/home/AFrangou/colorectal/battenberg_results/tumo",
#                        data[i,3],"_norm",data[i,4],
#                        "/",data[i,3],"_submission.sh")
# }
# write.table(script,"/home/AFrangou/colorectal/automated_firstbb_runs_v7_NAP.sh",row.names=F,col.names=F,quote=F)
# # change write permissions for these samples so DCW and NAP can write
# for (i in 1:nrow(data)) {
#   system(paste0("chmod -R 775 ~/colorectal/battenberg_results/tumo",data[i,3],"_norm",data[i,4]))
# }

#
#
# #
# # testicular
# # samples=read.csv("/re_gecip/cancer_testicular/0.Testicular_Landscape/0.AnalysisResults/01_Labkey_data/00_labkey_data_tables_v7/Clare_new_samples_Aug2019.csv",
#                  # stringsAsFactors=F)
# samples=read.csv("/re_gecip/cancer_testicular/0.Testicular_Landscape/0.AnalysisResults/01_Labkey_data/00_labkey_data_tables_v7/labkey_testicular_29Jul2019.csv",
#                                   stringsAsFactors=F)
# resultsdir = "/home/AFrangou/testicular_new/battenberg_results/"
# participantid = samples$Participant.Id
# tumourplatekey = samples$Tumour.Sample.Platekey
# normalplatekey = samples$Germline.Sample.Platekey
# tumourbam = samples$Tumour.Bam
# normalbam = samples$Germline.Bam
# gender = samples$Sex
# project.code = "re_gecip_cancer_testicular"
# bb1done = T
# ccubeexists = T
# vcfdir = "/re_gecip/cancer_testicular/0.Testicular_Landscape/0.AnalysisResults/02_Variant_calls/02_Filtered_variants/"
# ccubedir = "/re_gecip/cancer_testicular/0.Testicular_Landscape/0.AnalysisResults/04_Copy_number_calls/03_ccube/results/"
# ccubepurity = NA
# vcffilepath = paste0(vcfdir,"/",participantid,"/",normalplatekey,"_",tumourplatekey,".somatic.ILLUMINA_PASS.split.normalize.uniq.GNOMAD.APR_SOM.OCT_SOM.APR_GL.MAP.trf.closest_GL.closest_Gnomad.SEGDUP.closest_matched_gl_indel.spon.VEP.FILTERED.vcf.gz")
# data = cbind(resultsdir,
#              participantid,
#              tumourplatekey,
#              normalplatekey,
#              tumourbam,
#              normalbam,
#              gender,
#              project.code,
#              bb1done,
#              ccubeexists,
#              vcfdir,
#              ccubedir,
#              ccubepurity,
#              vcffilepath)
# data[which(data[,7]=="FEMALE"),7]="Female"
# data[which(data[,7]=="MALE"),7]="Male"
# for (i in 1:nrow(data)) {
#
#   CreateScriptsDirs(resultsdir=data[i,1],
#                     participantid=data[i,2],
#                     tumourplatekey=data[i,3],
#                     normalplatekey=data[i,4],
#                     tumourbam=data[i,5],
#                     normalbam=data[i,6],
#                     gender=data[i,7],
#                     project.code=data[i,8],
#                     bb1done=data[i,9],
#                     ccubeexists=data[i,10],
#                     vcfdir=data[i,11],
#                     ccubedir = data[i,12],
#                     ccubepurity = data[i,13],
#                     vcffilepath = data[i,14])
#
# }
# #
# numbersamples=nrow(data)
# script = matrix(nrow=numbersamples+3,ncol=1)
# script[1,1]="#!/usr/bin/env bash"
# script[2,1]=""
# script[3,1]=""
# for (i in 1:numbersamples) {
#   script[i+3,1]=paste0("sh ","/home/AFrangou/testicular_new/battenberg_results/tumo",
#                        data[i,3],"_norm",data[i,4],
#                        "/Submit_sample.sh")
# }
# write.table(script,"/home/AFrangou/testicular_new/automated_runs_v7.sh",row.names=F,col.names=F,quote=F)

# testicular reruns
# samples=read.csv("/re_gecip/cancer_testicular/0.Testicular_Landscape/0.AnalysisResults/01_Labkey_data/00_labkey_data_tables_v7/labkey_testicular_29Jul2019.csv",
#                  stringsAsFactors=F)
# resultsdir = "/home/AFrangou/testicular_new/battenberg_results/"
# participantid = samples$Participant.Id
# tumourplatekey = samples$Tumour.Sample.Platekey
# normalplatekey = samples$Germline.Sample.Platekey
# tumourbam = samples$Tumour.Bam
# normalbam = samples$Germline.Bam
# gender = samples$Sex
# project.code = "re_gecip_cancer_testicular"
# bb1done = T
# ccubeexists = F
# vcfdir = "/re_gecip/cancer_testicular/0.Testicular_Landscape/0.AnalysisResults/02_Variant_calls/02_Filtered_variants/"
# data = cbind(resultsdir,
#              participantid,
#              tumourplatekey,
#              normalplatekey,
#              tumourbam,
#              normalbam,
#              gender,
#              project.code,
#              bb1done,
#              ccubeexists,
#              vcfdir)
# data[which(data[,7]=="FEMALE"),7]="Female"
# data[which(data[,7]=="MALE"),7]="Male"
# for (i in 1:nrow(data)) {
#
#   CreateScriptsDirs(resultsdir=data[i,1],
#                     participantid=data[i,2],
#                     tumourplatekey=data[i,3],
#                     normalplatekey=data[i,4],
#                     tumourbam=data[i,5],
#                     normalbam=data[i,6],
#                     gender=data[i,7],
#                     project.code=data[i,8],
#                     bb1done=data[i,9],
#                     ccubeexists=data[i,10],
#                     vcf=data[i,11])
#
# }
# numbersamples=nrow(data)
# script = matrix(nrow=numbersamples+3,ncol=1)
# script[1,1]="#!/usr/bin/env bash"
# script[2,1]=""
# script[3,1]=""
# for (i in 1:numbersamples) {
#   script[i+3,1]=paste0("sh ","/home/AFrangou/testicular_new/battenberg_results/tumo",
#                        data[i,3],"_norm",data[i,4],
#                        "/Submit_sample.sh")
# }
# write.table(script,"/home/AFrangou/testicular_new/automated_runs_v7_Run1Run2.sh",row.names=F,col.names=F,quote=F)

# # testicular reruns
# samples=read.csv("/re_gecip/cancer_testicular/0.Testicular_Landscape/0.AnalysisResults/01_Labkey_data/00_labkey_data_tables_v7/labkey_testicular_29Jul2019.csv",
#                  stringsAsFactors=F)
# resultsdir = "/home/AFrangou/testicular_new/battenberg_results/"
# participantid = samples$Participant.Id
# tumourplatekey = samples$Tumour.Sample.Platekey
# normalplatekey = samples$Germline.Sample.Platekey
# tumourbam = samples$Tumour.Bam
# normalbam = samples$Germline.Bam
# gender = samples$Sex
# project.code = "re_gecip_cancer_testicular"
# bb1done = T
# ccubeexists = T
# vcfdir = "/re_gecip/cancer_testicular/0.Testicular_Landscape/0.AnalysisResults/02_Variant_calls/02_Filtered_variants/"
# ccubedir=paste0("/re_gecip/cancer_testicular/0.Testicular_Landscape/0.AnalysisResults/04_Copy_number_calls/03_ccube/results/",samples$Participant.Id,".tumo",samples$Tumour.Sample.Platekey,"_norm",samples$Germline.Sample.Platekey,".purity.tsv")
# data = cbind(resultsdir,
#              participantid,
#              tumourplatekey,
#              normalplatekey,
#              tumourbam,
#              normalbam,
#              gender,
#              project.code,
#              bb1done,
#              ccubeexists,
#              vcfdir,
#              ccubedir)
# data[which(data[,7]=="FEMALE"),7]="Female"
# data[which(data[,7]=="MALE"),7]="Male"
# for (i in 1:nrow(data)) {
#
#   CreateScriptsDirs(resultsdir=data[i,1],
#                     participantid=data[i,2],
#                     tumourplatekey=data[i,3],
#                     normalplatekey=data[i,4],
#                     tumourbam=data[i,5],
#                     normalbam=data[i,6],
#                     gender=data[i,7],
#                     project.code=data[i,8],
#                     bb1done=data[i,9],
#                     ccubeexists=data[i,10],
#                     vcf=data[i,11],
#                     ccubedir=data[i,12])
#
# }
#


#
#
# # sarcoma
# samples =read.csv("~/sarcoma/cancer_analysis_11_Aug_19.csv",stringsAsFactors = F)
# priorities = read.csv("~/sarcoma/2.6.1.priority.samples.simplified.csv",stringsAsFactors = F)
# touse =c()
# for (i in 1:nrow(priorities)) {
#
#   rows = which(samples[,1]==priorities[i,1])
#   touse = c(touse,rows)
#
# }
# samples = samples[touse,]
# resultsdir = "/home/AFrangou/sarcoma/2.6.battenberg/battenberg_results/"
# participantid = samples$Participant.Id
# tumourplatekey = samples$Tumour.Sample.Platekey
# normalplatekey = samples$Germline.Sample.Platekey
# tumourbam = samples$Tumour.Bam
# normalbam = samples$Germline.Bam
# gender = samples$Sex
# project.code = "re_gecip_cancer_sarcoma"
# bb1done = F
# ccubeexists = F
# vcfdir = "/re_gecip/cancer_sarcoma/0.sarcomaLandscape/0.3.variantCalling/0.3.1.Strelka"
# data = cbind(resultsdir,
#              participantid,
#              tumourplatekey,
#              normalplatekey,
#              tumourbam,
#              normalbam,
#              gender,
#              project.code,
#              bb1done,
#              ccubeexists,
#              vcfdir)
# data[which(data[,7]=="FEMALE"),7]="Female"
# data[which(data[,7]=="MALE"),7]="Male"
# for (i in 1:nrow(data)) {
#
#   CreateScriptsDirs(resultsdir=data[i,1],
#                     participantid=data[i,2],
#                     tumourplatekey=data[i,3],
#                     normalplatekey=data[i,4],
#                     tumourbam=data[i,5],
#                     normalbam=data[i,6],
#                     gender=data[i,7],
#                     project.code=data[i,8],
#                     bb1done=data[i,9],
#                     ccubeexists=data[i,10],
#                     vcf=data[i,11])
#
# }

# numbersamples=nrow(data)
# script = matrix(nrow=numbersamples+3,ncol=1)
# script[1,1]="#!/usr/bin/env bash"
# script[2,1]=""
# script[3,1]=""
# for (i in 1:numbersamples) {
#   script[i+3,1]=paste0("sh ","/home/AFrangou/sarcoma/2.6.battenberg/battenberg_results/tumo",
#                        data[i,3],"_norm",data[i,4],
#                        "/",data[i,3],"_submission.sh")
# }
# write.table(script,"/home/AFrangou/sarcoma/2.6.battenberg/automated_runs_prioritysamples_BB1.sh",row.names=F,col.names=F,quote=F)

#
#
# # # alona's
# samples = read.csv("~/haem/battenberg/cancer_analysis_2019-07-29_17-38-58.csv",stringsAsFactors=F)
# platekeys = read.table("~/haem/battenberg/samples_for_alona",stringsAsFactors=F)
# samples = samples[match(platekeys[,1],samples$Tumour.Sample.Platekey),]
#
# resultsdir = "/home/AFrangou/haem/battenberg/battenberg_results/"
# participantid = samples$Participant.Id
# tumourplatekey = samples$Tumour.Sample.Platekey
# normalplatekey = samples$Germline.Sample.Platekey
# tumourbam = samples$Tumour.Bam
# normalbam = samples$Germline.Bam
# gender = samples$Sex
# project.code = "re_gecip_cancer_haem"
# bb1done = T
# ccubeexists = T
# vcfdir = NA
# ccubedir = NA
# ccubepurity = samples$Tumour.Purity
# vcffilepath = samples$Tumour.Sv.Vcf
# vcffilepath = gsub("SV.","",vcffilepath)
# data = cbind(resultsdir,
#              participantid,
#              tumourplatekey,
#              normalplatekey,
#              tumourbam,
#              normalbam,
#              gender,
#              project.code,
#              bb1done,
#              ccubeexists,
#              vcfdir,
#              ccubedir,
#              ccubepurity,
#              vcffilepath)
# data[which(data[,7]=="FEMALE"),7]="Female"
# data[which(data[,7]=="MALE"),7]="Male"
# for (i in 1:nrow(data)) {
#
#   CreateScriptsDirs(resultsdir=data[i,1],
#                     participantid=data[i,2],
#                     tumourplatekey=data[i,3],
#                     normalplatekey=data[i,4],
#                     tumourbam=data[i,5],
#                     normalbam=data[i,6],
#                     gender=data[i,7],
#                     project.code=data[i,8],
#                     bb1done=data[i,9],
#                     ccubeexists=data[i,10],
#                     vcfdir=data[i,11],
#                     ccubedir = data[i,12],
#                     ccubepurity = as.numeric(data[i,13]),
#                     vcffilepath = data[i,14])
#
# }

# numbersamples=nrow(data)
# script = matrix(nrow=numbersamples+3,ncol=1)
# script[1,1]="#!/usr/bin/env bash"
# script[2,1]=""
# script[3,1]=""
# for (i in 1:numbersamples) {
#   script[i+3,1]=paste0("sh ","/home/AFrangou/haem/battenberg/battenberg_results/tumo",
#                        data[i,3],"_norm",data[i,4],
#                        "/Submit_sample.sh")
# }
# write.table(script,"/home/AFrangou/haem/battenberg/battenberg_results/alonas_rest.sh",row.names=F,col.names=F,quote=F)

# # endometrial v7 and others not yet run
# samples = read.csv("~/endometrial/endo_cancer_analysis_2019-09-16_17-15-21.csv",stringsAsFactors=F)
# samples = samples[which(samples$Preparation.Method=="FF"),]
# resultsdir = "/home/AFrangou/endometrial/battenberg_results/"
# participantid = samples$Participant.Id
# tumourplatekey = samples$Tumour.Sample.Platekey
# normalplatekey = samples$Germline.Sample.Platekey
# tumourbam = samples$Tumour.Bam
# normalbam = samples$Germline.Bam
# gender = samples$Sex
# project.code = "re_gecip_cancer_colorectal"
# bb1done = F
# ccubeexists = F
# vcfdir = "/re_gecip/cancer_ovarian/endometrial_current/analyses_v7/create_filtered_vcf_files_DC/output/"
# ccubedir = NA
# data = cbind(resultsdir,
#              participantid,
#              tumourplatekey,
#              normalplatekey,
#              tumourbam,
#              normalbam,
#              gender,
#              project.code,
#              bb1done,
#              ccubeexists,
#              vcfdir,
#              ccubedir)
# data[which(data[,7]=="FEMALE"),7]="Female"
# data[which(data[,7]=="MALE"),7]="Male"
# # remove those already run
# done = system(paste0("ls ~/endometrial/battenberg_results/tumo*/M-CallSubclones/*subclones.txt"),
#               intern=T)
# doneplatekeys = sapply(done,function(x){strsplit(x,"Subclones/")[[1]][2]})
# doneplatekeys = sapply(doneplatekeys,function(x){strsplit(x,"_subclones")[[1]][1]})
# toremove = match(doneplatekeys,data[,3])
# toremove = toremove[-which(is.na(toremove))]
# data = data[-toremove,]
# data = data[-c(3,5),] # couple of previously listed now probably removed samples
#
# for (i in 1:nrow(data)) {
#
#   CreateScriptsDirs(resultsdir=data[i,1],
#                     participantid=data[i,2],
#                     tumourplatekey=data[i,3],
#                     normalplatekey=data[i,4],
#                     tumourbam=data[i,5],
#                     normalbam=data[i,6],
#                     gender=data[i,7],
#                     project.code=data[i,8],
#                     bb1done=data[i,9],
#                     ccubeexists=data[i,10],
#                     vcfdir = data[i,11],
#                     ccubedir = data[i,12],
#                     ccubepur = NA,
#                     vcffilepath = NA)
#
# }
#
# numbersamples=nrow(data)
# script = matrix(nrow=numbersamples+3,ncol=1)
# script[1,1]="#!/usr/bin/env bash"
# script[2,1]=""
# script[3,1]=""
# for (i in 1:numbersamples) {
#   script[i+3,1]=paste0("sh ","/home/AFrangou/endometrial/battenberg_results/tumo",
#                        data[i,3],"_norm",data[i,4],
#                        "/",data[i,3],"_submission.sh")
# }
# write.table(script,"/home/AFrangou/endometrial/battenberg_results/final_110.sh",row.names=F,col.names=F,quote=F)
#




#


#  kill some jobs
# x=read.csv("~/testicular_new/testictokill5",stringsAsFactors=F)
# test=lapply(x,function(y){strsplit(y,"AFr")})
# test=unlist(test[[1]])
# test=test[seq(1,length(test),by=2)]
# samples=paste0("bkill ",test)
# script = matrix(nrow=length(samples)+3,ncol=1)
# script[1,1]="#!/usr/bin/env bash"
# script[2,1]="module load cluster"
# script[3,1]=""
# for (i in 1:length(samples)) {
#   script[i+3,1]=samples[i]
# }
# write.table(script,"/home/AFrangou/testicular_new/kill_hanging_testicular_jobs.sh",row.names=F,col.names=F,quote=F)

# list of jobs to kill if the config file has an NA in the FIXED_PARAMS_RHO field
# bjobs | grep "Aug 17" | grep PEND  > ~/colorectal/jobstokill
# x=read.csv("~/colorectal/fixed_params_reruns.sh",stringsAsFactors=F,sep=" ")
# dirs = sapply(x[,2],function(y){strsplit(y,"_norm")[[1]][1]})
# tumourstokill = sapply(dirs,function(x){strsplit(x,"tumo")[[1]][2]})
# for (i in 1:nrow(tumourstokill)) {
#   grep
# }

# better way
# samples = dir("~/colorectal/battenberg_results/")
# samples = samples[grep("_norm",samples)]
# dirs = paste0("~/colorectal/battenberg_results/",samples,"/")
# tumour = sapply(dirs,function(x){strsplit(x,"tumo")[[1]][2]})
# tumours = sapply(tumour,function(x){strsplit(x,"_norm")[[1]][1]})
# tokill =
# for (i in 1:length(dirs)) {
#   config = read.table(paste0(dirs[i],tumours[i],"_configfile.txt"),stringsAsFactors=F)
#   fixedrholine = config[grep("FIXED_RHO=",config[,1]),]
#   fixedrho = strsplit(fixedrholine,"FIXED_RHO=")[[1]][2]
#   if (is.na(fixedrho)) {
#
#   }
# }


