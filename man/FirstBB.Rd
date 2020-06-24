\name{CreateScriptsDirs}
\alias{CreateScriptsDirs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make scripts for all Battenberg and DPClust runs with options for number and type of reruns
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CreateScriptsDirs(resultsdir, participantid, tumourplatekey, normalplatekey, tumourbam, normalbam, gender, project.code)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultsdir}{
%%     ~~Describe \code{resultsdir} here~~
}
  \item{participantid}{
%%     ~~Describe \code{participantid} here~~
}
  \item{tumourplatekey}{
%%     ~~Describe \code{tumourplatekey} here~~
}
  \item{normalplatekey}{
%%     ~~Describe \code{normalplatekey} here~~
}
  \item{tumourbam}{
%%     ~~Describe \code{tumourbam} here~~
}
  \item{normalbam}{
%%     ~~Describe \code{normalbam} here~~
}
  \item{gender}{
%%     ~~Describe \code{gender} here~~
}
  \item{project.code}{
%%     ~~Describe \code{project.code} here~~
}
  \item{bb1done}{
%%     ~~Describe \code{project.code} here~~
}
  \item{ccubeexists}{
%%     ~~Describe \code{project.code} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (resultsdir, participantid, tumourplatekey, normalplatekey,
    tumourbam, normalbam, gender, project.code)
{
    participantdir = paste0(resultsdir, participantid)
    dir.create(participantdir)
    sampledir = paste0(participantdir, "/tumo", tumourplatekey,
        "_norm", normalplatekey)
    dir.create(sampledir)
    dir.create(paste0(sampledir, "/A-GetAlleleCounts"))
    dir.create(paste0(sampledir, "/B-RunBAFLogR"))
    dir.create(paste0(sampledir, "/C-RunGCcorrect"))
    dir.create(paste0(sampledir, "/D-GenerateImputeInputFromAlleleFrequencies"))
    dir.create(paste0(sampledir, "/E-RunImpute"))
    dir.create(paste0(sampledir, "/F-CombineImputeOutputs"))
    dir.create(paste0(sampledir, "/G-GetHaplotypedBAFs"))
    dir.create(paste0(sampledir, "/H-CleanUp"))
    dir.create(paste0(sampledir, "/I-PlotHaplotypedData"))
    dir.create(paste0(sampledir, "/J-CombineBAFfiles"))
    dir.create(paste0(sampledir, "/K-SegmentBAFphased"))
    dir.create(paste0(sampledir, "/L-FitCopyNumber"))
    dir.create(paste0(sampledir, "/M-CallSubclones"))
    dir.create(paste0(sampledir, "/N-wrapup"))
    dir.create(paste0(sampledir, "/logs"))
    config = matrix(nrow = 34, ncol = 1)
    config[1, 1] = paste0("RUN_DIR=", resultsdir)
    config[2, 1] = paste0("LOG_DIR=", sampledir, "/logs/")
    config[3, 1] = paste0("OUTPUT_DIR=", sampledir)
    config[4, 1] = paste0("SAMPLE_PATH=", sampledir)
    config[5, 1] = "PIPELINE_DIR=/home/AFrangou/battenberg-lsf-pipeline"
    config[6, 1] = paste0("TUMOURNAME=", tumourplatekey)
    config[7, 1] = paste0("IS_MALE=", is.male)
    config[8, 1] = paste0("NORMALNAME=", normalplatekey)
    config[9, 1] = paste0("TUMOURBAM=", tumourbam)
    config[10, 1] = paste0("NORMALBAM=", normalbam)
    config[11, 1] = "PLATFORM_GAMMA=1"
    config[12, 1] = "PHASING_GAMMA=1"
    config[13, 1] = "SEGMENTATION_GAMMA=10"
    config[14, 1] = "CLONALITY_DIST_METRIC=0"
    config[15, 1] = "ASCAT_DIST_METRIC=1"
    config[16, 1] = "MIN_PLOIDY=1.6"
    config[17, 1] = "MAX_PLOIDY=4.8"
    config[18, 1] = "MIN_RHO=0.13"
    config[19, 1] = "MAX_RHO=1.02"
    config[20, 1] = "MIN_GOODNESS_OF_FIT=0.63"
    config[21, 1] = "BALANCED_THRESHOLD=0.51"
    config[22, 1] = "IMPUTEINFOFILE=/home/AFrangou/ALL_100G_phase1integrated_v3_impute/impute_info.txt"
    config[23, 1] = "IMPUTE_EXE=/home/AFrangou/impute_v2.3.2_x86_64_static/impute2"
    config[24, 1] = "SEED=123"
    config[25, 1] = "MAX_CN_STATE=250"
    config[26, 1] = "SV_BREAKPOINTS_FILE=NA"
    config[27, 1] = "PROBLEMLOCI=/home/AFrangou/probloci_270415.txt"
    config[28, 1] = "G1000_PREFIX_HG38=/home/AFrangou/battenberg_1000genomesloci2012_v3_hg38/hg38_alleleFiles_feb22_chr"
    config[29, 1] = "G1000_PREFIX_AC_HG38=/home/AFrangou/battenberg_1000genomesloci2012_v3_hg38/hg38_lociFiles_feb22_chr"
    config[30, 1] = "G1000_PREFIX=/home/AFrangou/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr"
    config[31, 1] = "G1000_PREFIX_AC=/home/AFrangou/runallelecounter/1000genomesloci2012/1000genomesloci2012_chr"
    config[32, 1] = "GCCORRECTPREFIX=/home/AFrangou/battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_filled_chr_"
    config[33, 1] = "MIN_NORMAL_DEPTH=10"
    config[34, 1] = "ALLELECOUNTER=/home/AFrangou/bin/alleleCounter"
    write.table(config, paste0(sampledir, "/", tumourplatekey,
        "_configfile.txt"), quote = F, col.names = F, row.names = F)
    header = matrix(nrow = 51, ncol = 1)
    header[1, 1] = "#!/usr/bin/env bash"
    header[2, 1] = "CONFIG=$1"
    header[4, 1] = "# Get allele frequencies"
    header[5, 1] = "# Tumour"
    header[6, 1] = paste0("bsub -J\"GetAlleleFrequenciesTumour_",
        tumourplatekey, "[1-23]\" -q gecip -P ", project.code,
        " -o ", sampledir, "/logs/GetAlleleFrequenciesTumour.\%J.\%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/GetAlleleCounts_tumour.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[7, 1] = "# Normal"
    header[8, 1] = paste0("bsub -J\"GetAlleleFrequenciesNormal_",
        normalplatekey, "[1-23]\" -q gecip -P ", project.code,
        " -o ", sampledir, "/logs/GetAlleleFrequenciesNormal.\%J.\%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/GetAlleleCounts_normal.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[10, 1] = "# Remove 'chr' from alleleCounter output, array"
    header[11, 1] = paste0("bsub -w\"done(GetAlleleFrequenciesTumour_",
        tumourplatekey, "[1-23]) && done(GetAlleleFrequenciesNormal_",
        normalplatekey, "[1-23])\" -J\"RunRemoveCHR_", tumourplatekey,
        "\" -q gecip -P ", project.code, " -o ", sampledir, "/logs/RemoveCHR.\%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/RemoveCHR.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[13, 1] = paste0("# Convert alleleCounter output back to hg37")
    header[14, 1] = paste0("bsub -w\"done(RunRemoveCHR_", tumourplatekey,
        ")\" -J\"RunSwitchback_hg38_to_hg37_", tumourplatekey,
        "[1-23]\" -q gecip -P ", project.code, " -o ", sampledir,
        "/logs/Switchback_hg38_to_hg37.\%J.\%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/Switchback_hg38_to_hg37.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[16, 1] = "# Get BAF and logR"
    header[17, 1] = paste0("bsub -w\"done(RunSwitchback_hg38_to_hg37_",
        tumourplatekey, "[1-23])\" -R\"select[mem>28000] rusage[mem=28000]\" -M28000 -J\"RunBAFLogR_",
        tumourplatekey, "\" -q gecip -P ", project.code, " -o ",
        sampledir, "/logs/RunBAFLogR.\%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/RunBAFLogR.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[19, 1] = "# RunGCCorrect_wgs"
    header[20, 1] = paste0("bsub -w \"done(RunBAFLogR_", tumourplatekey,
        ")\" -R\"select[mem>35000] rusage[mem=35000]\" -M35000 -J\"runGCcorrect_",
        tumourplatekey, "\" -q gecip -P ", project.code, " -o ",
        sampledir, "/logs/runGCcorrect", tumourplatekey, ".\%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/RunGCcorrect_wgs.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[22, 1] = "# Perform phasing"
    header[23, 1] = "# GenerateImputeInputFromAlleleFrequencies"
    header[24, 1] = paste0("bsub -w \"done(RunSwitchback_hg38_to_hg37_",
        tumourplatekey, "[1-23])\" -R\"select[mem>4000] rusage[mem=4000]\" -M4000 -J\"GenerateImputeInputs_",
        tumourplatekey, "[1-23]\" -q gecip -P ", project.code,
        " -o ", sampledir, "/logs/GenerateImputeInputs.\%J.\%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/GenerateImputeInputFromAlleleFrequencies.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[26, 1] = "# RunImpute"
    header[27, 1] = paste0("bsub -w \"done(GenerateImputeInputs_",
        tumourplatekey, "[*])\" -R\"select[mem>8000] rusage[mem=8000]\" -M8000 -J\"Impute_",
        tumourplatekey, "[1-23]\" -q gecip -P ", project.code,
        " -o ", sampledir, "/logs/Impute.\%J.\%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/RunImpute.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[29, 1] = "# CombineImputeOutputs"
    header[30, 1] = paste0("bsub -w\"done(Impute_", tumourplatekey,
        "[*])\" -R\"select[mem>4000] rusage[mem=4000]\" -M4000 -J\"CombineImputeOutputs_",
        tumourplatekey, "[1-23]\" -q gecip -P ", project.code,
        " -o ", sampledir, "/logs/CombineImputeOutputs.\%J.\%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/CombineImputeOutputs.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[32, 1] = "# GetHaplotypeBAFs"
    header[33, 1] = paste("bsub -w\"done(CombineImputeOutputs_",
        tumourplatekey, "[*])\" -R\"select[mem>2000] rusage[mem=2000]\" -M2000 -J\"GetHaplotypedBAFs_",
        tumourplatekey, "[1-23]\" -q gecip -P ", project.code,
        " -o ", sampledir, "/logs/GetHaplotypedBAFs.\%J.\%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/GetHaplotypedBAFs.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[35, 1] = "# CleanUp"
    header[36, 1] = paste0("bsub -w\"done(GetHaplotypedBAFs_",
        tumourplatekey, "[*])\" -J\"CleanUp_", tumourplatekey,
        "[1-23]\" -q gecip -P ", project.code, " -o ", sampledir,
        "/logs/CleanUp.\%J.\%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/CleanUp.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[38, 1] = "# PlotHaplotypedData"
    header[39, 1] = paste0("bsub -w\"done(GetHaplotypedBAFs_",
        tumourplatekey, "[*])\" -R\"select[mem>4000] rusage[mem=4000]\" -M4000 -J\"PlotHaplotypedData_",
        tumourplatekey, "[1-23]\" -q gecip -P ", project.code,
        " -o ", sampledir, "/logs/PlotHaplotypedData.\%J.\%I.out /home/AFrangou/battenberg_lsf_pipeline/steps/PlotHaplotypedData.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[41, 1] = "# CombineBAFfiles"
    header[42, 1] = paste0("bsub -w\"done(GetHaplotypedBAFs_",
        tumourplatekey, ")\" -R\"select[mem>4000] rusage[mem=4000]\" -M4000 -J\"CombineBAFfiles_",
        tumourplatekey, "\" -q gecip -P ", project.code, " -o ",
        sampledir, "/logs/CombineBAFfiles.\%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/CombineBAFfiles.sh ",
        resultsdir, sampledir, "/", tumourplatekey, "_configfile.txt")
    header[44, 1] = "# Segmentation and copy number calling"
    header[45, 1] = paste0("bsub -w\"done(CombineBAFfiles_",
        tumourplatekey, ")\" -R\"select[mem>4000] rusage[mem=4000]\" -M4000 -J\"SegmentBAFphased_",
        tumourplatekey, "\" -q gecip -P ", project.code, " -o ",
        sampledir, "/logs/SegmentBAFphased.\%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/SegmentBAFphased.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[47, 1] = paste0("bsub -w\"done(SegmentBAFphased_",
        tumourplatekey, ") && done(runGCcorrect_", tumourplatekey,
        ")\" -R\"select[mem>15900] rusage[mem=15900]\" -M15900 -J\"FitCopyNumber_",
        tumourplatekey, "\" -q gecip -P ", project.code, " -o ",
        sampledir, "/logs/FitCopyNumber.\%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/FitCopyNumber.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[49, 1] = paste0("bsub -w\"done(FitCopyNumber_", tumourplatekey,
        ")\" -R\"select[mem>15900] rusage[mem=15900]\" -M15900 -J\"CallSubclones_",
        tumourplatekey, "\" -q gecip -P ", project.code, " -o ",
        sampledir, "/logs/CallSubclones.\%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/CallSubclones.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[51, 1] = paste0("bsub -w\"done(CallSubclones_", tumourplatekey,
        ")\" -R\"select[mem>500] rusage[mem=500]\" -M500 -J\"wrapup_",
        tumourplatekey, "\" -q gecip -P ", project.code, " -o ",
        sampledir, "/logs/wrapup.\%J.out /home/AFrangou/battenberg_lsf_pipeline/steps/wrapup.sh ",
        sampledir, "/", tumourplatekey, "_configfile.txt")
    header[which(is.na(header[, 1])), ] = ""
    write.table(header, paste(sampledir, "/", tumourplatekey,
        "_submission.sh", sep = ""), quote = F, col.names = F,
        row.names = F)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
