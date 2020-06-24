run.beagle5=function(beaglejar, vcfpath, reffile, outpath, plinkfile) {
  cmd <- paste0("java -Xmx15g -jar ", beaglejar,
                " gt=", vcfpath, " ref=", reffile, " out=", outpath,
                " map=", plinkfile, " nthreads=1 window=40 overlap=4 impute=false")
  system(cmd, wait = T)
}


