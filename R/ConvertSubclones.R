#-------------------------------------------------------------------------------#
# Function to convert the subclones file from hg37 back to hg38 (single sample)
#-------------------------------------------------------------------------------#

# sampledir = "/home/AFrangou/colorectal/battenberg_results/tumoX_normY"
# participantid = "ID"
# tumourplatekey = "tumoX"
# normalplatekey = "normY"
# run = 1
# lowerclonal = 0.9
# upperclonal = 1
# # conversionfilesroot="~/re_gecip/cancer_testicular/04_Copy_number_calls/software/01_Battenberg/battenberg_1000genomesloci2012_v3_hg38/chr"
# conversionfilesroot="/home/AFrangou/battenberg_1000genomesloci2012_v3_hg38/chr"
# conversionfilesend="_hg37-38alleles_full_feb22.txt"

ConvertSubclones <- function(sampledir,participantid,tumourplatekey,normalplatekey,run,conversionfilesroot,conversionfilesend) {

        # loop assigns the conversion files for hg37 to hg38 positions to the name eg 'chr1' for use
        for (i in 1:23) {

                assign(paste0("chr",i),read.table(paste0(conversionfilesroot,i,conversionfilesend),header=T,stringsAsFactors=F))
                print(paste("chr",i,"conversion file assigned to object"))

        }

        # get run directories
        if (run==1) {
                callsubclonesdir = "/"
        } else if (run==2) {
                callsubclonesdir = "/R-BB2/"
        } else if (run==3) {
                callsubclonesdir = "/U-BB3/"
        } else if (run==4) {
                callsubclonesdir = "/X-BB4/"
        } else if (run=="WGD") {
                callsubclonesdir = "/ZZ-BB_WGD/"
        }


        subclonesfile = paste0(sampledir,callsubclonesdir,"M-CallSubclones/",tumourplatekey,"_subclones.txt")

        # original subclones file
        sub=read.table(subclonesfile,header=T,stringsAsFactors=F)
        print(paste("Read subclonesfile",subclonesfile))

        # empty new subclones file
        wholenewsub=as.data.frame(matrix(nrow=0,ncol=ncol(sub)))
        colnames(wholenewsub)=colnames(sub)

        # fill in new subclones file
        for (i in 1:23) {

                pos=get(paste0("chr",i))

                if (i==23) {
                        subsub=sub[which(sub[,1]=="X"),]
                } else {
                        subsub=sub[which(sub[,1]==i),]
                }

                newsub=subsub
                newsub[,2]=pos[match(subsub[,2],pos[,2]),4]
                newsub[,3]=pos[match(subsub[,3],pos[,2]),4]
                wholenewsub=rbind(wholenewsub,newsub)

        }

        # output message
        print(paste0("Converting subclones file for tumo",tumourplatekey,"_norm",normalplatekey," to hg38"))

        #  write subclones conversion
        write.table(wholenewsub,paste0(sampledir,callsubclonesdir,"O-Postprocessing/",tumourplatekey,"_subclones_hg38.txt"),row.names=F,quote=F,sep="\t")

        # let user know
        print(paste0("Sample ",tumourplatekey,"-",normalplatekey," converted to hg38"))

}



