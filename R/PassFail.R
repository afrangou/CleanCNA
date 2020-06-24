#--------------------------------------------------------#
# Make a PASS or FAIL file to notify submission script
#--------------------------------------------------------#

PassFail <- function(sampledir,participantid,tumourplatekey,normalplatekey,run) {

        # get run directories
        if (run==1) {
          assessmentdir = "Q-AssessBB1DPC1/"
        } else if (run==2) {
          assessmentdir = "Q-AssessBB2DPC2/"
        } else if (run==3) {
          assessmentdir = "Q-AssessBB3DPC3/"
        }

        # read in metrics table from previous run
        metricsfile = paste0(sampledir,assessmentdir,tumourplatekey,"_metrics_run",run)
        metrics = read.csv(paste0(sampledir,assessmentdir,tumourplatekey,"_metrics_run",run))

        if (metrics$passeddpc=="Yes") {
                write.csv("PASS",paste0(sampledir,assessmentdir,"PASS"),quote=F,row.names = F)
        } else if (metrics$passeddpc=="No") {
                write.csv("FAIL",paste0(sampledir,assessmentdir,"FAIL"),quote=F,row.names = F)
        }


}
