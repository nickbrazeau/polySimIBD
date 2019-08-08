#' @title make mixture
#'
#' @description Make mixture from list of samples
#'
#' @param strains list;
#' @param strains.prop list;
#' @param pos numeric vector; genomic coordinates
#' @param dp numeric; read depth
#'
#' @export

#TODO add sample name
make_mixture <- function(strains = list(), strains.prop = list(),
                         pos, dp){

  if(length(strains) != length(strains.prop)){
    stop("strains and strain proportions must be of same length")
  }


  strains.mat <- unname(cbind.data.frame(strains))
  gt <- apply(strains.mat, 1, function(x){
    samp <- sample(x, dp, prob = unlist(strains.prop), replace = T)
    refad <- sum(unlist(samp) == 0)
    altad <- sum(unlist(samp) == 2)
    return( paste0("./.:", refad, ",", altad, ":", dp) ) # for now keep GT blank
  })

  gt <- cbind(FORMAT = "GT:AD:DP", gt)

  fix <- as.matrix(data.frame(CHROM = "contig1", POS = pos,
                              ID = ".", REF = ".", ALT = ".", QUAL = ".", FILTER = ".", INFO = "."))

  meta <- append("##fileformat=VCFv4.3",
                 "##SimpleSim=This vcf was created using the simpleIBD approach")

  # write out new vcfRobj
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  return(newvcfR)

}

