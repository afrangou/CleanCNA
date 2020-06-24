#------------------------------------------------------------------------------------------------#
# Function to generate the DPinput file from Battenberg output and sample vcfs which is then used in DPClust
#------------------------------------------------------------------------------------------------#

GenerateDPinput <- function(sampledir,participantid,tumourplatekey,normalplatekey,gender,vcffilepath,rhoandpsifilepath,subcloneshg38filepath,run) {

        # definitions
        nucleotides = c("A","C","G","T")

        # print samplename to logfile
        samplename=paste0("tumo",tumourplatekey,"_norm",normalplatekey)
        print(paste("Sample name:",samplename))
        print(run)
        print(subcloneshg38filepath)

        # get run directories
        if (run==1) {
                dpcoveralldir = "/P-DPC1/"
                dpcdir = "/P-DPC1/DPClust/"
                dpcprepdir = "/P-DPC1/DPPrep/"
                subdir = "/O-Postprocessing/"
                fitcopydir = "/L-FitCopyNumber/"
                puritydir = "/M-CallSubclones/"
                assessmentdir = "/Q-AssessBB1DPC1/"
        } else if (run==2) {
                dpcoveralldir = "/S-DPC2/"
                dpcdir = "/S-DPC2/DPClust/"
                dpcprepdir = "/S-DPC2/DPPrep/"
                subdir = "/R-BB2/O-Postprocessing/"
                fitcopydir = "/R-BB2/L-FitCopyNumber/"
                puritydir = "/R-BB2/M-CallSubclones/"
                assessmentdir = "/T-AssessBB2DPC2/"
        } else if (run==3) {
                dpcoveralldir = "/V-DPC3/"
                dpcdir = "/V-DPC3/DPClust/"
                dpcprepdir = "/V-DPC3/DPPrep/"
                subdir = "/U-BB3/O-Postprocessing/"
                fitcopydir = "/U-BB3/L-FitCopyNumber/"
                puritydir = "/U-BB3/M-CallSubclones/"
                assessmentdir = "/W-AssessBB3DPC3/"
        } else if (run==4) {
                dpcoveralldir = "/Y-DPC4/"
                dpcdir = "/Y-DPC4/DPClust/"
                dpcprepdir = "/Y-DPC4/DPPrep/"
                subdir = "/X-BB4/O-Postprocessing/"
                fitcopydir = "/X-BB4/L-FitCopyNumber/"
                puritydir = "/X-BB4/M-CallSubclones/"
                assessmentdir = "/Z-AssessBB4DPC4/"
        } else if (run=="WGD") {
                if (file.exists(paste0(sampledir,"/",samplename,"/Z-AssessBB4DPC4/",tumourplatekey,"_peaks_output.Rdata"))) {
                        dpcoveralldir = "/Y-DPC4/"
                        dpcdir = "/Y-DPC4/DPClust/"
                        dpcprepdir = "/Y-DPC4/DPPrep/"
                        subdir = "/X-BB4/O-Postprocessing/"
                        fitcopydir = "/X-BB4/L-FitCopyNumber/"
                        puritydir = "/X-BB4/M-CallSubclones/"
                        assessmentdir = "/Z-AssessBB4DPC4/"
                } else if (file.exists(paste0(sampledir,"/",samplename,"/W-AssessBB3DPC3/",tumourplatekey,"_peaks_output.Rdata"))) {
                        dpcoveralldir = "/V-DPC3/"
                        dpcdir = "/V-DPC3/DPClust/"
                        dpcprepdir = "/V-DPC3/DPPrep/"
                        subdir = "/U-BB3/O-Postprocessing/"
                        fitcopydir = "/U-BB3/L-FitCopyNumber/"
                        puritydir = "/U-BB3/M-CallSubclones/"
                        assessmentdir = "/W-AssessBB3DPC3/"
                } else if (file.exists(paste0(sampledir,"/",samplename,"/T-AssessBB2DPC2/",tumourplatekey,"_peaks_output.Rdata"))) {
                        dpcoveralldir = "/S-DPC2/"
                        dpcdir = "/S-DPC2/DPClust/"
                        dpcprepdir = "/S-DPC2/DPPrep/"
                        subdir = "/R-BB2/O-Postprocessing/"
                        fitcopydir = "/R-BB2/L-FitCopyNumber/"
                        puritydir = "/R-BB2/M-CallSubclones/"
                        assessmentdir = "/T-AssessBB2DPC2/"
                } else if (file.exists(paste0(sampledir,"/",samplename,"/Q-AssessBB1DPC1/",tumourplatekey,"_peaks_output.Rdata"))) {
                        dpcoveralldir = "/P-DPC1/"
                        dpcdir = "/P-DPC1/DPClust/"
                        dpcprepdir = "/P-DPC1/DPPrep/"
                        subdir = "/O-Postprocessing/"
                        fitcopydir = "/L-FitCopyNumber/"
                        puritydir = "/M-CallSubclones/"
                        assessmentdir = "/Q-AssessBB1DPC1/"
                }
        }

        # check vcf exists and turn it into format for DPClust, else print that it doesn't exist and the sample is exiting
        if (file.exists(vcffilepath)) {

                # load vcf
                snvs=read.table(vcffilepath,stringsAsFactors=F)
                print(vcffilepath)
                print("vcfloaded")

                # select only those variants which have passed all filters and which are SNVs (not indels)
                snvs=snvs[which(snvs[,7]=="PASS" & nchar(snvs[,4])==1 & nchar(snvs[,5])==1),]

                # make new snvs table for dpinput
                newsnvs=snvs[,c(1,2,2,4,5,10,11)]
                newsnvs[,2]=newsnvs[,3]-1

                # remove stuff not assigned to a chr
                if ("chrX" %in% snvs[,1]) {
                  newsnvs=newsnvs[1:max(which(newsnvs[,1]=="chrX")),]
                }
                # replace chr1 with 1 etc
                newsnvs[,1]=gsub("chr","",newsnvs[,1])

                # turn counts from vcf into required format
                # this is for Dan's vcf format which includes both normal and tumour counts (at some point was using Alona's almost-fully-filtered
                # files which only had counts for the tumour)

                # normal depth
                normdepth=sapply(newsnvs[,6],function(x){strsplit(x,":")[[1]][1]})
                # counts in normal of A,C,G,T
                normA=sapply(newsnvs[,6],function(x){strsplit(x,":")[[1]][5]})
                normA=sapply(normA,function(x){strsplit(x,",")[[1]][1]})
                normC=sapply(newsnvs[,6],function(x){strsplit(x,":")[[1]][6]})
                normC=sapply(normC,function(x){strsplit(x,",")[[1]][1]})
                normG=sapply(newsnvs[,6],function(x){strsplit(x,":")[[1]][7]})
                normG=sapply(normG,function(x){strsplit(x,",")[[1]][1]})
                normT=sapply(newsnvs[,6],function(x){strsplit(x,":")[[1]][8]})
                normT=sapply(normT,function(x){strsplit(x,",")[[1]][1]})
                # tumour depth
                tumdepth=sapply(newsnvs[,7],function(x){strsplit(x,":")[[1]][9]})
                # counts in tumour of A,C,G,T
                tumA=sapply(newsnvs[,7],function(x){strsplit(x,":")[[1]][5]})
                tumA=sapply(tumA,function(x){strsplit(x,",")[[1]][1]})
                tumC=sapply(newsnvs[,7],function(x){strsplit(x,":")[[1]][6]})
                tumC=sapply(tumC,function(x){strsplit(x,",")[[1]][1]})
                tumG=sapply(newsnvs[,7],function(x){strsplit(x,":")[[1]][7]})
                tumG=sapply(tumG,function(x){strsplit(x,",")[[1]][1]})
                tumT=sapply(newsnvs[,7],function(x){strsplit(x,":")[[1]][8]})
                tumT=sapply(tumT,function(x){strsplit(x,",")[[1]][1]})
                # create file type for normal
                normpos=cbind(newsnvs[,c(1,3,4,5),],0,0,0,0,0)
                colnames(normpos)=c("Chromosome","Position","REF","ALT","count_A","count_C","count_G","count_T","total_depth")
                normpos[,5]=normA
                normpos[,6]=normC
                normpos[,7]=normG
                normpos[,8]=normT
                normpos[,9]=normdepth
                # create file type for tumour
                tumpos=cbind(newsnvs[,c(1,3,4,5),],0,0,0,0,0)
                colnames(tumpos)=c("Chromosome","Position","REF","ALT","count_A","count_C","count_G","count_T","total_depth")
                tumpos[,5]=tumA
                tumpos[,6]=tumC
                tumpos[,7]=tumG
                tumpos[,8]=tumT
                tumpos[,9]=tumdepth

                # file containing all SNV positions and counts of all alleles at these positions
                normalpositioncounts=normpos[,c(1,2,5:9)]
                tumourpositioncounts=tumpos[,c(1,2,5:9)]
                # make directory for prep for DPClust
                dir.create(paste0(sampledir,dpcoveralldir))
                dir.create(paste0(sampledir,dpcprepdir))
                dir.create(paste0(sampledir,dpcdir))

                # write count info
                write.table(normalpositioncounts,
                            paste0(sampledir,dpcprepdir,normalplatekey,"_snv_counts_normal.txt"),
                            row.names=F,
                            quote=F,
                            sep="\t")
                write.table(tumourpositioncounts,
                            paste0(sampledir,dpcprepdir,tumourplatekey,"_snv_counts_tumour.txt"),
                            row.names=F,
                            quote=F,
                            sep="\t")

                # file containing just the loci for input into the function below
                tumourloci=tumpos[,c(1:4)]
                write.table(tumourloci,paste0(sampledir,
                                              dpcprepdir,
                                              tumourplatekey,"_loci.txt"),
                                              row.names=F,
                                              col.names=F,
                                              quote=F,
                                              sep="\t")

                # define filenames for function
                loci_file=paste0(sampledir,dpcprepdir,tumourplatekey,"_loci.txt")
                allele_frequencies_file=paste0(sampledir,dpcprepdir,tumourplatekey,"_snv_counts_tumour.txt")
                cellularity_file=paste0(sampledir,fitcopydir,tumourplatekey,"_rho_and_psi.txt")
                subclone_file=paste0(sampledir,subdir,tumourplatekey,"_subclones_hg38.txt")
                output_file=paste0(sampledir,dpcprepdir,samplename,"_DPinput.txt")
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
                                           output_file)


        } else {

                print(paste("No vcf available - exiting"))

        }

}


