\name{RunDPClust_Beagle_nf}
\alias{RunDPClust_Beagle_nf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run DPClust (set up for Nextflow)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RunDPClust_Beagle_nf(tumourplatekey, normalplatekey, rhoandpsifilepath, dpinfofilepath, output_folder)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tumourplatekey}{
%%     ~~Describe \code{tumourplatekey} here~~
}
  \item{normalplatekey}{
%%     ~~Describe \code{normalplatekey} here~~
}
  \item{rhoandpsifilepath}{
%%     ~~Describe \code{rhoandpsifilepath} here~~
}
  \item{dpinfofilepath}{
%%     ~~Describe \code{dpinfofilepath} here~~
}
  \item{output_folder}{
%%     ~~Describe \code{output_folder} here~~
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
function (tumourplatekey, normalplatekey, rhoandpsifilepath,
    dpinfofilepath, output_folder)
{
    nucleotides = c("A", "C", "G", "T")
    samplename = paste0("tumo", tumourplatekey, "_norm", normalplatekey)
    print(paste0(tumourplatekey, " ", normalplatekey, " preparing"))
    if (file.exists(rhoandpsifilepath) & file.exists(dpinfofilepath)) {
        cellularity = read.table(rhoandpsifilepath, header = T,
            stringsAsFactors = F, sep = "\t")[2, 1]
        DP.info = read.table(dpinfofilepath, sep = "\t", header = T,
            stringsAsFactors = F)
        DP.info = DP.info[!is.na(DP.info$subclonal.fraction) &
            DP.info$subclonal.fraction > 0, ]
        mutCount = array(DP.info$mut.count, c(nrow(DP.info),
            1))
        WTCount = array(DP.info$WT.count, c(nrow(DP.info), 1))
        totalCopyNumber = array(DP.info$subclonal.CN, c(nrow(DP.info),
            1))
        copyNumberAdjustment = array(DP.info$no.chrs.bearing.mut,
            c(nrow(DP.info), 1))
        mutation.copy.number = array(DP.info$mutation.copy.number,
            c(nrow(DP.info), 1))
        no.iters = 1000
        output_folder = output_folder
        print(paste("Using cellularity file", rhoandpsifilepath))
        print(paste("Using DPinput file", dpinfofilepath))
        print(paste("Writing output to", output_folder))
        DirichletProcessClustering(mutCount = mutCount, WTCount = WTCount,
            totalCopyNumber = totalCopyNumber, copyNumberAdjustment = copyNumberAdjustment,
            mutation.copy.number = mutation.copy.number, cellularity = cellularity,
            output_folder = output_folder, no.iters = no.iters,
            no.iters.burn.in = ceiling(no.iters/5), subsamplesrun = samplename,
            samplename = samplename, conc_param = 1, cluster_conc = 5,
            mut.assignment.type = 1, most.similar.mut = NA, mutationTypes = NA,
            min.frac.snvs.cluster = NA, max.considered.clusters = 30)
        print(paste0("Clustering has run for ", tumourplatekey,
            "_", normalplatekey))
    }
    else {
        print(paste0(rhoandpsifilepath, " or ", dpinfofilepath,
            " not found, exiting."))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
