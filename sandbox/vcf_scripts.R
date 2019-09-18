#' @title Join two vcfR objects that have the same CHROM/POS Columns
#'
#' @description This assumes that we have two samples that came from the same original VCF and were split
#' for post-processing and are now being re-merged. Will not run if CHROMPOS for both samples are not identical
#'
#' @param vcflist list of vcfR object
#' @importFrom magrittr %>%
#' @export

join_simvcfs <- function(vcfRlist){

  warning("This join function overwrites any information in the
           REF, ALT,
           ID, QUAL, FILTER, INFO
          fields under the assumption that this data isn't meaningul for sims")

  chrompos <- lapply(vcfRlist, function(x){ return( vcfR::getFIX(x)[,1:2] )})
  if(length(unique(chrompos)) != 1 ){
    stop("The CHROMPOS Genomic Coordinates for the VCFs do not match. Cannot join.")
  }

  vcfRlist.gtmat <- lapply(vcfRlist, function(x){
    chrompos <- vcfR::getFIX(x)[,1:2]
    gt <- x@gt
    ret <- cbind.data.frame(chrompos, gt)

    })

  vcfR.CHROMPOSgtmat <- vcfRlist.gtmat %>%
    purrr::reduce(dplyr::left_join, by = c("CHROM", "POS", "FORMAT"))

  # now make in vcfR object

  fix <- as.matrix(data.frame(CHROM = vcfR.CHROMPOSgtmat[,1], POS = vcfR.CHROMPOSgtmat[,2],
                              ID = ".", REF = ".", ALT = ".", QUAL = ".", FILTER = ".", INFO = "."))

  meta <- c("##fileformat=VCFv4.3",
            "##SimpleSim=This vcf was created using the simpleIBD approach",
            "##SimpleSimJoin=This vcf was joined using an internal function from the R package polySimIBD")

  gt <- vcfR.CHROMPOSgtmat %>%
    dplyr::select(-c("CHROM", "POS")) %>%
    as.matrix(.)

  # write out new vcfRobj
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  return(newvcfR)


}




#' @title Converit SimDataPop to vcfR Object
#'
#' @description This converts an object of \{simPopData} to a \{vcfR} object
#' @rdname simPopData
#' @inheritParams make_mixture
#'
#' @param simPopData
#' @importFrom magrittr %>%
#' @export

simPopData2vcfR <- function(simPopData,
                            dp, dpstd){


  # double for loop to make gt matrix
  gtmat <- apply(simPopData$genmat[,3:ncol(simPopData$genmat)], 2,
               function(x){
                 sapply(x, function(y){
                   dpit <- floor( rnorm(1, mean = dp, sd = dpstd) )
                   if( y == "0"){
                     refad = dpit; altad = 0
                   } else if( y == "2"){
                     refad = 0; altad = dpit
                   }
                   return( paste0("./.:", refad, ",", altad, ":", dpit) ) # for now keep GT blank
                 })
               }
  )

  # now coerce to vcfR object
  meta <- c("##fileformat=VCFv4.3",
            "##SimpleSim=This vcf was created using the simPopData approach from the R package polySimIBD")

  fix <- as.matrix(data.frame(CHROM = simPopData$genmat$CHROM,
                              POS = simPopData$genmat$POS,
                              ID = ".", REF = ".", ALT = ".",
                              QUAL = ".", FILTER = ".", INFO = "."))

  gt <- cbind(FORMAT = "GT:AD:DP", gtmat)


  # write out new vcfRobj
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  return(newvcfR)


}








