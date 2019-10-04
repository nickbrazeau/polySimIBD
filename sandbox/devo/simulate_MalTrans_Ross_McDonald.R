
#' @title Simulate IBD-WF
#'
#' @desciption Simulate a population with recombination that approximates the Wright Fisher Process
#' and tracks haplotype identity by descent
#'
#' @param chrompos matrix <factor><numeric>; the genomic coordinates for chromosome and position of the sites
#' @param Nes numeric; The starting effective population size
#' @param rho numeric; expected recombination rate
#' @param k numeric; number of generation to simulate forward
#' @param l.cEIR numeric; lambda that a mosquito bites two different hosts (case EIR)
#' @param l.ooc numeric; lambda for oocyst count
#' @param l.spor numeric; lambda for the number of sporozoites \class{hmp} that will be passed in a cotransmission event
#' @param l.superinfxn numeric; lambda for the superinfxn

#'@details Data is simulated under the following framework:
#'   \enumerate{
#'        \item Start with a set of unique parasite haploid genotypes (n = \code{Nes})
#'        \item Haploid genotypes are repeated based on the population census size (n = \code{Nc}).
#'        \item Encode haplotypes with a \emph{int}
#'        \item Infection the individuals that make up the population census (n = \code{Nc}) with the parasite genomes
#'        \item Simulate forward in time ...
#'   }
#'
#'
#' @export

simulate_MalTrans_WF <- function(chrompos, Nc, Nes, rho, k, l.cEIR, l.ooc, l.spor, l.superinfxn){

  #.....................................................................................................
  # Initial Population
  #.....................................................................................................

  Nes.haplos <- lapply(1:Nes, function(x) return( new("simhaplo") ))
  # loop through effective pop
  for(i in 1:length(Nes.haplos)){
    Nes.haplos[[i]]@haploint <- rep(i, nrow(chrompos))
  }


  popdf <- lapply(1:(Nc*k), function(x) return(list()))
  dim(popdf) <- c(k, Nc)
  popdf[1, ] <- Nes.haplos

  #.....................................................................................................
  # Simulate Forward in Time
  #.....................................................................................................

    # first get what calls are contributing to next generation
    for(i in 2:nrow(popdf)){
      for(j in 1:ncol(popdf)){
        #.......................
        # Mosquito External
        #.......................
        # Do we bite one or multipe hosts
        bites <- 1 + rpois(1, l.cEIR)
        #.......................
        # Mosquito Midgut
        #.......................
        strns <- sample(x = 1:ncol(popdf), size = bites, replace = T)
        strns <- popdf[i-1, strns]
        # in first row, only haplotypes
        # but in second, we can have MOI which is a list of lists, so unnest
        if(i-1 > 1){
         strns <- unlist(strns)
        }


        oocystcount <- 1 + rpois(1, l.ooc)
        sporozoites <- list()

        for(o in 1:oocystcount){
          draws <- sample(strns, 2, replace = T)
          hmps12 <- makecrossover(p1 = draws[[1]], p2 = draws[[2]], rho = rho, chrompos = chrompos)
          hmps34 <- makecrossover(p1 = draws[[1]], p2 = draws[[2]], rho = rho, chrompos = chrompos)
          hmps <- c(hmps12, hmps34)
          sporozoites <- append(sporozoites, hmps)
        }
        #.......................
        # Mosquito Salivary Glands
        #.......................
        nhmps <- 1 + rpois(1, l.spor)
        ret <- sample(sporozoites, size = nhmps, replace = T)

        #.......................
        # Does superinfxn occur
        #.......................
        superinfxns <- rpois(1, l.superinfxn) # here we allow zeros
        if(superinfxns != 0){
          superinfxns <- sample(x = 1:ncol(popdf), size = superinfxns, replace = T)
          superinfxns <- popdf[i-1, superinfxns] # sample from previous generation
          #.......................
          # Second Mosquito Salivary Glands
          #.......................
          nspor.superinfxn <- 1 + rpois(1, l.spor)
          superinfxns <- sample(superinfxns, size = nspor.superinfxn, replace = T)

          ret <- append(ret, superinfxns)
        }

        popdf[i,j] <- list( ret )

    }
  } # end nested for loop

  return(popdf)
}
