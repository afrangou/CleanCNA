\name{RunShapeIt2}
\alias{RunShapeIt2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run ShapeIt2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RunShapeIt2(sampledir, tumourplatekey, imputeinfofile, is.male, chrom, filename.input, filename.output, filename.shapeit2, seed = as.integer(Sys.time()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampledir}{
%%     ~~Describe \code{sampledir} here~~
}
  \item{tumourplatekey}{
%%     ~~Describe \code{tumourplatekey} here~~
}
  \item{imputeinfofile}{
%%     ~~Describe \code{imputeinfofile} here~~
}
  \item{is.male}{
%%     ~~Describe \code{is.male} here~~
}
  \item{chrom}{
%%     ~~Describe \code{chrom} here~~
}
  \item{filename.input}{
%%     ~~Describe \code{filename.input} here~~
}
  \item{filename.output}{
%%     ~~Describe \code{filename.output} here~~
}
  \item{filename.shapeit2}{
%%     ~~Describe \code{filename.shapeit2} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
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
function (sampledir, tumourplatekey, imputeinfofile, is.male,
    chrom, filename.input, filename.output, filename.shapeit2,
    seed = as.integer(Sys.time()))
{
    impute.info = read.table(imputeinfofile, stringsAsFactors = F)
    colnames(impute.info) = c("chrom", "impute_legend", "genetic_map",
        "impute_hap", "impute_sample", "start", "end", "is_par")
    if (is.male) {
        impute.info = impute.info[impute.info$is_par == 1, ]
    }
    chr_names = unique(impute.info$chrom)
    if (!is.na(chrom)) {
        impute.info = impute.info[impute.info$chrom == chr_names[chrom],
            ]
    }
    n = nrow(impute.info)
    for (r in 1:n) {
        if (chrom == 23) {
            filename.input.ok <- paste0(filename.input, ".filt",
                r)
            legend <- read.delim(impute.info$impute_legend[r],
                sep = " ", header = T, stringsAsFactors = F)
            impute <- read.table(paste0(sampledir, "/D-GenerateImputeInputFromAlleleFrequencies/",
                tumourplatekey, "_impute_input_chr23.txt"), stringsAsFactors = F)
            impute <- impute[impute[[3]] \%in\% legend$position,
                ]
            indels.inds <- which(nchar(impute[[4]]) > 1 | nchar(impute[[5]]) >
                1)
            indels <- impute[indels.inds, ]
            gen <- impute[-indels.inds, ]
            write.table(gen, file = paste0(filename.input.ok,
                ".gen"), quote = F, sep = " ", row.names = F,
                col.names = F)
            write.table(indels, file = paste0(filename.input.ok,
                ".gen.missing"), quote = F, sep = " ", row.names = F,
                col.names = F)
        }
        else {
            filename.input.ok <- filename.input
            impute <- read.table(paste0(sampledir, "/D-GenerateImputeInputFromAlleleFrequencies/",
                tumourplatekey, "_impute_input_chr", chrom, ".txt"),
                stringsAsFactors = F)
            indels.inds <- which(nchar(impute[[4]]) > 1 | nchar(impute[[5]]) >
                1)
            indels <- impute[indels.inds, ]
            gen <- impute[-indels.inds, ]
            write.table(gen, file = paste0(filename.input.ok,
                ".gen"), quote = F, sep = " ", row.names = F,
                col.names = F)
            write.table(indels, file = paste0(filename.input.ok,
                ".gen.missing"), quote = F, sep = " ", row.names = F,
                col.names = F)
        }
        cmd = paste(filename.shapeit2, "-G", filename.input.ok,
            "-M", impute.info$genetic_map[r], "-R", impute.info$impute_hap[r],
            impute.info$impute_legend[r], impute.info$impute_sample[r],
            "-O", ifelse(n == 1, filename.output, paste(filename.output,
                r, sep = ".")), ifelse(!impute.info$is_par[r],
                "--chrX", ""), "--thread 1", "--seed", seed,
            "--force")
        system(cmd, wait = T)
    }
    if (n > 1) {
        filenames.haps <- paste(filename.output, 1:n, "haps",
            sep = ".")
        haps <- do.call("rbind", sapply(filenames.haps, read.delim,
            sep = " ", header = F, stringsAsFactors = F, simplify = F))
        haps <- haps[!duplicated(haps[[3]]), ]
        haps <- haps[order(haps[[3]]), ]
        write.table(haps, file = paste0(filename.output, ".haps"),
            quote = F, sep = " ", row.names = F, col.names = F)
        system(paste0("cp ", filename.output, ".1.sample ", filename.output,
            ".sample"))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
