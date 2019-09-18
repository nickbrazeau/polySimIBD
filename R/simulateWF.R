
#' @title Simulate IBD-WF
#'
#' @description Simulate a population with recombination that approximates the Wright Fisher Process
#' and tracks haplotype identity by descent
#'
#' @param chrompos matrix <factor><numeric>; the genomic coordinates for chromosome and position of the sites
#' @param Nc numeric; The starting census size
#' @param Nes numeric; The starting effective population size
#' @param rho numeric; expected recombination rate
#' @param k numeric; number of generation to simulate forward
#'
#'@details Data is simulated under the following framework:
#'   \enumerate{
#'        \item Haploid genotypes are drawn from the \code{Nes}
#'        \item Haploid genotypes are repeated based on the census size, \code{Nc}
#'        \item Encode haplotypes with a \emph{int}
#'        \item Simulate forward in time from the census population \code{Nc} by drawing two samples
#'        and allowing sexual recombination to occur.
#'   }
#'
#'   The forward in time simulation is performed by randomly drawing two strains from the previous generation,
#'   where strains can either be homologous or heterologous haplotypes depending on the haplotype frequencies,
#'   and allowing for mating to occur. Within the formed zygote, recombination events occur at a rate proportional
#'   to the recombination rate and distance along the genome: \eqn{1-e^{-rho*d}}. This process is modeled with a
#'   transition matrix and determines whether two samples enter a state of IBD or not. Once the recombination block
#'   has been determined, the four potential haploid meiotic products are created.
#'   TODO We then randomly pick one of the four which may not be ideal
#'
#' @export

simulate_IBD_WF <- function(chrompos, Nc, Nes, rho, k){


  Nes.haplos <- lapply(1:Nes, function(x) return( new("simhaplo") ))
  # loop through effective pop
  for(i in 1:length(Nes.haplos)){
    Nes.haplos[[i]]@haploint <- rep(i, nrow(chrompos))
  }

  # based on Census Size, inflate population
  Nc.draws <- sample(x = 1:Nes, size = Nc, replace = T)
  Nc.haplos <- Nes.haplos[ Nc.draws ]

  #.....................................................................................
  # Simulate Forward in Time
  #.....................................................................................
  nedf <- lapply(1:(Nc*k), function(x) return(new("simhaplo")))
  dim(nedf) <- c(k, Nc)
  nedf[1, ] <- Nc.haplos

  for(i in 2:nrow(nedf)){
    for(j in 1:ncol(nedf)){
      # pick two samples from previous generation
      draws <- nedf[i-1, sample(x = 1:ncol(nedf), size = 2)]
      # make crossover/recombination event
      child <- makecrossover(p1 = draws[[1]], p2 = draws[[2]], rho = rho, chrompos = chrompos)
      nedf[i,j] <- list(child)
    }
  }

  return(nedf)

}












