#------------------------------------------------------------------------------------------------#
# Function to generate .sample file per tumour for running ShapeIt2
#------------------------------------------------------------------------------------------------#

# A .sample file is a tsv file containing 4 columns:
# 'ID1' and 'ID2' (which can both be the tumour ID),
# 'missing' (which can just equal 0) and 'sex' (1 for MALE and 2 for FEMALE),
# 'and 3 rows: the header, a row defining what type of data is contained in each column, and the actual sample information.


# gender="MALE"


PrepShapeIt2 <- function(sampledir,participantid,tumourplatekey,normalplatekey,gender) {

        samplefile = matrix(nrow=2,ncol=4)
        colnames(samplefile)=c("ID_1","ID_2","missing","sex")
        samplefile[1:2,]=c(0,0,0,1)
        if (gender=="FEMALE") {samplefile[1:2,4]=2}
        samplefile[2,1:2]=tumourplatekey

        for (chr in 1:23) {
                if (chr==23) {
                        for (filt in 1:2) {
                                write.table(samplefile,paste0(sampledir,"/E-RunShapeIt/",
                                                        tumourplatekey,"_shapeit_input_chr",chr,".filt",filt,".sample"),
                                                        quote=F,
                                                        row.names = F,
                                                        sep="\t")
                        }
                        if (gender=="FEMALE") {
                              write.table(samplefile,paste0(sampledir,"/E-RunShapeIt/",
                                                            tumourplatekey,"_shapeit_input_chr",chr,".filt3.sample"),
                                          quote=F,
                                          row.names = F,
                                          sep="\t")
                        }

                } else {
                        write.table(samplefile,paste0(sampledir,"/E-RunShapeIt/",
                                                        tumourplatekey,"_shapeit_input_chr",chr,".sample"),
                                                        quote=F,
                                                        row.names = F,
                                                        sep="\t")
                }
        }

}

# PrepShapeIt2(sampledir,participantid,tumourplatekey,normalplatekey,gender)
