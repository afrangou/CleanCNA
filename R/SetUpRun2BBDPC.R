#-----------------------------------------------------------#
# Based on the assessment of BB+DPC, set up a rerun (or not)
#-----------------------------------------------------------#

# This function collects the assessment of the sample based on the metrics produced
# And sets up a reun of FitCopyNumber and CallSubclones and DPPrep2 and DPC2 if a rerun is needed
# This is done by using the identified closest-to-clonal cluster and adjusting the purity and ploidy from this

# sampledir = "/home/AFrangou/colorectal/battenberg_results/tumoX_normY/"
#
# # sampledir = "~/re_gecip/cancer_colorectal/analysisResults/2.chromosomeAberrations/I.Battenberg/01_Battenberg_original/ID/tumoX_normY/"
# participantid = "ID"
# tumourplatekey = "X"
# normalplatekey = "Y"
# run = 2


SetUpRun2BBDPC <- function(sampledir,participantid,tumourplatekey,normalplatekey,run,ccubedir,ccubepurity) {

                # get run directories
                if (run==1) {
                        # First run of DPC
                        assessmentdir = "/Q-AssessBB1DPC1/"
                } else if (run==2) {
                        # Second run with new purity from clonal peak from DPC
                        assessmentdir = "/T-AssessBB2DPC2/"
                } else if (run==3) {
                        # Third run with Ccube purity
                        assessmentdir = "/W-AssessBB3DPC3/"
                } else if (run==4) {
                        # Fourth run with purity from peaks
                        assessmentdir = "/Z-AssessBB4DPC4/"
                }

                # some info
                print(paste("This is checking run",run))
                print(paste("Checking status of run with",sampledir,assessmentdir,tumourplatekey,"_metrics_run",run))

                # read in bb/dpc metrics table from previous run
                metricsfile = paste0(sampledir,assessmentdir,tumourplatekey,"_metrics_run",run)
                metrics = read.csv(paste0(sampledir,assessmentdir,tumourplatekey,"_metrics_run",run))

                print(metrics)

                # old purity and ploidy
                oldpurity = metrics$purity
                oldploidy = metrics$ploidy
                print(paste("Old purity was",oldpurity))
                print(paste("Old ploidy was",oldploidy))

                if (run==1) {

                        # Take DPClust assessment and peaks assessment

                        # Add new rho and psi from the DPClust assessment to the config file
                        # calculate new purity and ploidy using (old purity) * (position cluster with highest CCF)
                        newpurity = oldpurity * metrics$position.topccfcluster
                        newploidy = ((oldpurity * oldploidy) + 2*(newpurity - oldpurity))/newpurity

                        print(paste("New purity is",newpurity))
                        print(paste("New ploidy is",newploidy))

                        # add these new parameters to the config file
                        system(paste0("echo RHO_NEW=",newpurity,
                                      " >> ",
                                      sampledir,"/",
                                      tumourplatekey,"_configfile.txt"))
                        print("Added RHO_NEW to config file for run2")
                        system(paste0("echo PSI_NEW=",newploidy,
                                      " >> ",
                                      sampledir,"/",
                                      tumourplatekey,"_configfile.txt"))
                        print("Added PSI_NEW to config file for run2")

                } else if (run==2) {

                        # take new purity from Ccube and calculate new ploidy accordingly (as before)
                        if (!is.na(ccubepurity)) {
                                newpurity = as.numeric(ccubepurity)
                        } else if (!is.na(ccubedir)) {
                                newpurity = as.numeric(read.csv(paste0(ccubedir,participantid,
                                                                       ".tumo",
                                                                       tumourplatekey,
                                                                       "_norm",
                                                                       normalplatekey,
                                                                       ".purity.tsv")))
                        }

                        newploidy = ((oldpurity * oldploidy) + 2*(newpurity - oldpurity))/newpurity

                        print(paste("New purity is",newpurity))
                        print(paste("New ploidy is",newploidy))

                        system(paste0("echo RHO_NEW_NEW=",newpurity,
                                      " >> ",
                                      sampledir,"/",
                                      tumourplatekey,"_configfile.txt"))
                        print("Added RHO_NEW_NEW to config file for run3")
                        system(paste0("echo PSI_NEW_NEW=",newploidy,
                                      " >> ",
                                      sampledir,"/",
                                      tumourplatekey,"_configfile.txt"))
                        print("Added PSI_NEW_NEW to config file for run3")


                } else if (run==3) {

                        # take new purity from peak calling and calculate new ploidy
                        # loads Rdata list called 'peaks'
                        load(paste0(sampledir,
                                    "/W-AssessBB3DPC3/",
                                    tumourplatekey,"_peaks_output.Rdata"))

                        # old ploidy comes from ccube run
                        oldploidy = oldploidy
                        # old purity comes from peaks 'peaks' Rdata table - should match purity from output of ccube run
                        # (if they used the last ccube run every time for the peaks assessment)
                        # but looks like they used the first run
                        oldpurity = peaks$summary$purity.old
                        newpurity = peaks$summary$purity.new
                        newploidy = ((oldpurity * oldploidy) + 2*(newpurity - oldpurity))/newpurity

                        print(paste("New purity is",newpurity))
                        print(paste("New ploidy is",newploidy))

                        system(paste0("echo RHO_PEAKS=",newpurity,
                                      " >> ",
                                      sampledir,"/",
                                      tumourplatekey,"_configfile.txt"))
                        print("Added RHO_PEAKS to config file")
                        system(paste0("echo PSI_PEAKS=",newploidy,
                                      " >> ",
                                      sampledir,"/",
                                      tumourplatekey,"_configfile.txt"))
                        print("Added PSI_PEAKS to config file")

                }

}
