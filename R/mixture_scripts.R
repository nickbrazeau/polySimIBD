#' @title make mixture
#'
#' @description Make mixture from list of samples
#'
#' @param strains named list; sample names with a vector of 0,2 corresponding to haplotypes
#' @param strains.prop numeric vector; proportion of the strains
#' @param CHROMPOS matrix <character, numeric>; genomic coordinates corresponding to chromosome and position
#' @param dp numeric; average read depth per loci
#' @param dpstd numeric; standard deviation around read depth
#' @param samplename character; sample name for resulting mixture
#'
#' @export

make_mixture <- function(strains = list(),
                         strains.prop = list(),
                         CHROMPOS,
                         dp, dpstd,
                         samplename){

  if(length(strains) != length(strains.prop)){
    stop("strains and strain proportions must be of same length")
  }

  dp <- floor( rnorm(n = nrow(CHROMPOS), mean = dp, sd = dpstd) )

  strains.mat <- unname(cbind.data.frame(strains))
  gt <- apply(strains.mat, 1, function(x){
    dpit <- floor( rnorm(1, mean = dp, sd = dpstd) )
    samp <- sample(x, dpit, prob = strains.prop, replace = T)
    refad <- sum(unlist(samp) == 0)
    altad <- sum(unlist(samp) == 2)
    return( paste0("./.:", refad, ",", altad, ":", dpit) ) # for now keep GT blank
  })

  gt <- cbind(FORMAT = "GT:AD:DP", gt)

  # add sample name
  colnames(gt)[2] <- samplename

  fix <- as.matrix(data.frame(CHROM = CHROMPOS[,1], POS = CHROMPOS[,2],
                              ID = ".", REF = ".", ALT = ".", QUAL = ".", FILTER = ".", INFO = "."))

  meta <- c("##fileformat=VCFv4.3",
            "##SimpleSim=This vcf was created using the simpleIBD approach")

  # write out new vcfRobj
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  return(newvcfR)

}

