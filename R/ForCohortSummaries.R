#
# ################ R ################
# library(CleanCNA)
#
# qc = "/filepath/to/QualityControl.tsv" # produced by rerun pipeline
# segfile_dir = "/dir/containing/final/subclones/and/cellularity/files/" # a directory where all the final *subclones.txt and *cellularity_ploidy.txt files have been put
# segfile_name = "REF_NAME_FOR_COHORT_OR_TUMOUR_TYPE" # a name for the cohort (whatever you choose)
# output_dir = "/dir/for/results/to/write/to/" # where the output files will write to
# filestub = segfile_dir
#
# CleanCNA:::CollateSubclones(qc,
#                             segfile_dir,
#                             segfile_name)
#
# CleanCNA:::CodeSegments(segfile_dir,
#                         segfile_name)
#
#
# ################ Command line ################
# load R module
# load bedops module
# load bedtools module
#
# declare -a arr=(homdel loh otherloss nochange gain biggain)
#
# filestub=/dir/containing/final/subclones/and/cellularity/files/ # same as filestub/segfile_dir in R above
# segfile_name=REF_NAME_FOR_COHORT_OR_TUMOUR_TYPE # same as segfile_name in R above
#
# for i in ${arr[@]}
# do
# bedtools genomecov \
# -i ${filestub}${segfile_name}_${i}.bed \
# -g /path/to/hg38_chr_lengths.txt \ # this is at /mnt/lustre/users/afrangou/public/code/CakeTin.nf/data/hg38.chrom.sizes
# -bga > ${filestub}${segfile_name}_${i}_regions.out
# done
#
# ################ R ################
# CleanCNA:::LabelGenomecov(filestub,
#                           segfile_name)
#
#
# ################ Command line ################
# bedops --partition ${filestub}${segfile_name}_homdel_regions.out \
# ${filestub}${segfile_name}_loh_regions.out \
# ${filestub}${segfile_name}_otherloss_regions.out \
# ${filestub}${segfile_name}_nochange_regions.out \
# ${filestub}${segfile_name}_gain_regions.out \
# ${filestub}${segfile_name}_biggain_regions.out > \
# ${filestub}${segfile_name}_bedops
#
# ################ R ################
# CleanCNA:::CleanBedopsOutput(filestub,
#                              segfile_name)
#
# ################ Command line ################
# bedtools intersect \
# -a ${filestub}${segfile_name}_bedops \
# -b ${filestub}${segfile_name}_homdel_regions.out ${filestub}${segfile_name}_loh_regions.out ${filestub}${segfile_name}_otherloss_regions.out ${filestub}${segfile_name}_nochange_regions.out ${filestub}${segfile_name}_gain_regions.out ${filestub}${segfile_name}_biggain_regions.out \
# -wb > ${filestub}${segfile_name}_partitioned_regions_all_CNA_types_overlapped.out
#
# ################ R ################
# CleanCNA:::GenomeWideStackedBar(filestub,
#                                 segfile_name)
#
# bed_file = "/filepath/to/bedfile/" # this file gives chr positions for genes. tab sep. 7 cols, chr,start,end,strand,ens,ens,genename (no colnames), we only use cols 1,2,3,7, others can be empty (but cols must be there)
# driver_file = "/filepath/to/list/of/drivers/called/for/cohort" # this needs columns named $TIER & $SYMBOL, tab sep file, doesn't matter which order columns are in. Tier 1 here refers to those we are def happy to call drivers, Tier 2 maybe a bit less certain - can be defined as you wish essentially.
#
# CleanCNA:::CNAsDriversHeatmapDendrogram(segfile_dir,
#                                         segfile_name,
#                                         bed_file,
#                                         driver_file)
#
# CleanCNA:::CNAsDriversComplexHeatmap(segfile_dir,
#                                      segfile_name,
#                                      plot_name,
#                                      bed_file,
#                                      driver_file,
#                                      tomelt)









