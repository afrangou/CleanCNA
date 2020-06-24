#--------------------------------------------------------------------------#
# Create scripts and directories for runs with provided purity and ploidy
#--------------------------------------------------------------------------#

CreateScriptsDirsFixedParams <- function(resultsdir,participantid,tumourplatekey,normalplatekey,tumourbam,normalbam,gender,project.code,fixed_rho,fixed_psi,vcfdir) {

        sampledir = paste0(resultsdir,"tumo",tumourplatekey,"_norm",normalplatekey)
        dir.create(sampledir)

        # create subdirs & subsubdirs
        dir.create(paste0(sampledir,"/X-BB4"))
        dir.create(paste0(sampledir,"/X-BB4/L-FitCopyNumber"))
        dir.create(paste0(sampledir,"/X-BB4/M-CallSubclones"))
        dir.create(paste0(sampledir,"/X-BB4/O-Postprocessing"))
        dir.create(paste0(sampledir,"/Y-DPC4"))
        dir.create(paste0(sampledir,"/Y-DPC4/DPPrep"))
        dir.create(paste0(sampledir,"/Y-DPC4/DPClust"))
        dir.create(paste0(sampledir,"/Z-AssessBB4DPC4"))

        # add fixed rho to the config file
        system(paste0("echo FIXED_RHO=",fixed_rho,
                      " >> ",
                      sampledir,"/",
                      tumourplatekey,"_configfile.txt"))
        # add fixed psi to the config file
        system(paste0("echo FIXED_PSI=",fixed_psi,
                      " >> ",
                      sampledir,"/",
                      tumourplatekey,"_configfile.txt"))
        # add vcffilepath to the config file
        system(paste0("echo VCFFILEPATH=",vcfdir,participantid,"/",normalplatekey,"_",tumourplatekey,
                      ".somatic.ILLUMINA_PASS.split.normalize.uniq.GNOMAD.APR_SOM.OCT_SOM.APR_GL.MAP.trf.closest_GL.closest_Gnomad.SEGDUP.closest_matched_gl_indel.spon.VEP.FILTERED.vcf.gz",
                      " >> ",
                      sampledir,"/",
                      tumourplatekey,"_configfile.txt"))

        # make script for this independent run
        # create script for third run if a rerun is needed (so Setup, FCN, CS, Convert, DPPrep, DPClust, Assess)
        header=matrix(nrow=21,ncol=1)
        header[1,1]="#!/usr/bin/env bash"
        header[2,1]="CONFIG=$1"
        header[3,1]=""
        header[4,1]=""
        header[5,1]="# Run BB FitCopyNumberFixedParam"
        header[6,1]=paste0('bsub -J"FitCopyNumberFixedParam_',tumourplatekey,'"',
                           ' -R"select[mem>15900] rusage[mem=15900]" -M15900',
                           ' -q gecip -P ',project.code,
                           ' -o ',sampledir,"/logs/FitCopyNumberFixedParam.%J.out",
                           ' -e ',sampledir,"/logs/FitCopyNumberFixedParam.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/FitCopyNumberFixedParam.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[7,1]=""
        header[8,1]="# Run BB CallSubclonesFixedParam"
        header[9,1]=paste0('bsub -w"done(FitCopyNumberFixedParam_',tumourplatekey,')"',
                           ' -J"CallSubclonesFixedParam_',tumourplatekey,'"',
                           ' -R"select[mem>15900] rusage[mem=15900]" -M15900',
                           ' -q gecip -P ',project.code,
                           ' -o ',sampledir,"/logs/CallSubclonesFixedParam.%J.out",
                           ' -e ',sampledir,"/logs/CallSubclonesFixedParam.%J.err",
                           ' /home/AFrangou/battenberg_lsf_pipeline/steps/CallSubclonesFixedParam.sh ',
                           sampledir,'/',tumourplatekey,'_configfile.txt')
        header[10,1]=""
        header[11,1]="# Convert subclones from CallSubclonesFixedParam"
        header[12,1]=paste0('bsub -w"done(CallSubclonesFixedParam_',tumourplatekey,')"',
                            ' -J"ConvertSubclonesFixedParam_',tumourplatekey,'"',
                            ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q gecip -P ',project.code,
                            ' -o ',sampledir,"/logs/ConvertSubclonesFixedParam.%J.out",
                            ' -e ',sampledir,"/logs/ConvertSubclonesFixedParam.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/ConvertSubclonesFixedParam.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[13,1]=""
        header[14,1]="# Run DPPrepFixedParam"
        header[15,1]=paste0('bsub -w"done(ConvertSubclonesFixedParam_',tumourplatekey,')"',
                            ' -J"DPPrepFixedParam_',tumourplatekey,'"',
                            ' -R"select[mem>8000] rusage[mem=8000]" -M8000',
                            ' -q gecip -P ',project.code,
                            ' -o ',sampledir,"/logs/DPPrepFixedParam.%J.out",
                            ' -e ',sampledir,"/logs/DPPrepFixedParam.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPPrepFixedParam.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[16,1]=""
        header[17,1]="# Run DPClustFixedParam"
        header[18,1]=paste0('bsub -w"done(DPPrepFixedParam_',tumourplatekey,')"',
                            ' -J"DPClustFixedParam_',tumourplatekey,'"',
                            ' -R"select[mem>16000] rusage[mem=16000]" -M16000',
                            ' -q gecip -P ',project.code,
                            ' -o ',sampledir,"/logs/DPClustFixedParam.%J.out",
                            ' -e ',sampledir,"/logs/DPClustFixedParam.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/DPClustFixedParam.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')
        header[19,1]=""
        header[20,1]="# Assess PASS/FAIL of DPClustFixedParam"
        header[21,1]=paste0('bsub -w"done(DPClustFixedParam_',tumourplatekey,')"',
                            ' -J"AssessBBDPCFixedParam_',tumourplatekey,'"',
                            ' -R"select[mem>4000] rusage[mem=4000]" -M4000',
                            ' -q gecip -P ',project.code,
                            ' -o ',sampledir,"/logs/AssessBBDPCFixedParam.%J.out",
                            ' -e ',sampledir,"/logs/AssessBBDPCFixedParam.%J.err",
                            ' /home/AFrangou/battenberg_lsf_pipeline/steps/AssessBBDPCFixedParam.sh ',
                            sampledir,'/',tumourplatekey,'_configfile.txt')

        write.table(header,paste0(sampledir,'/RunFixedParam.sh'),quote=F,col.names=F,row.names=F)


}

