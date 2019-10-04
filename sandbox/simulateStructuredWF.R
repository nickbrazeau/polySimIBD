
#' @title The Structured Wright Fisher Model for IBD
#'
#' @description Simulate a population with recombination that approximates the Structured Wright Fisher Process
#' and tracks haplotype identity by descent where individuals represent demes, such that within
#' deme Ne approximates the individual-level COI
#'
#' @param chrompos dataframe <factor><numeric>; the genomic coordinates for chromosome and position of the sites
#' @param Nes numeric; The starting effective population size
#' @param Noff numeric; Number of offspring to simulate in each generation, t. Note, this should be a large number ("infinity in a WF sense").
#' @param K numeric; The number of individuals (e.g. demes) to consider
#' @param m numeric; Probability of migration where m represents the probability of moving from \eqn{deme_{origin}} to  \eqn{deme_{new}} by \eqn{m*(1-1/K)}
#' @param coi numeric; The lambda of a right-shifted Poisson process, \eqn{1 + Pos(lambda)} representing the average COI of an individual deme
#' @param rho numeric; expected recombination rate
#' @param G numeric; number of generation to simulate forward
#' @param interference numeric; probability of a double crossover
#'
#'@details Data is simulated under the following framework:
#'   \enumerate{
#'        \item Haploid genotypes are drawn from the \code{Nes}
#'        \item Haploid genotypes are repeated based on the census size, \code{Nc}
#'        \item Encode haplotypes with a \emph{int}
#'        \item ... TODO>..
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

simulate_IBD_SWF <- function(chrompos, Nes, Noff, K, m, coi, rho, G, interference = 0){

  #.........................
  # assertions
  #.........................
  assert_dataframe(chrompos)
  assert_numeric(Nes)
  assert_numeric(Noff)
  assert_numeric(K)
  assert_numeric(m)
  assert_numeric(coi)
  assert_numeric(rho)
  assert_numeric(G)

  #..........................................................
  # Initial Generation
  #..........................................................
  current <- lapply(1:K, function(x) return(list()))
  dim(current) <- c(1, K)
  # initial COIs
  coi.int <- 1 + rpois(n = K, lambda = coi)

  # get haplos
  Nes.haplos <- lapply(1:Nes, function(x) return( new("simhaplo") ))
  # loop through effective pop to make init haplo objects
  for(i in 1:length(Nes.haplos)){
    Nes.haplos[[i]]@haploint <- rep(i, nrow(chrompos))
    Nes.haplos[[i]]@GA <- matrix(NA, nrow = G, ncol = nrow(chrompos))
    Nes.haplos[[i]]@GA[1, ] <- i # self
  }

  # fill demes
  for(k in 1:K){
    current[[k]] <- sample(Nes.haplos, coi.int[k], replace = T)
  }


  #..........................................................
  # Subsequent Generations
  #..........................................................
  for(g in 2:G){
    proposed <- lapply(1:K, function(x) return(list()))
    dim(proposed) <- c(1, K)

    #.........................
    # Perform recombination
    #.........................
    for(k in 1:K){ # for each deme
      for(c in 1:coi.int[[k]]){
        for(n in 1:Noff){ # make this many recombinants
        parents <- sample(current[[k]], size = 2, replace = T)
        child <- makecrossover(p1 = parents[[1]], p2 = parents[[2]], chrompos = chrompos, rho = rho, interference = interference)[[1]]

        # child inherits topology, G1- from parent under assumption parent lineage exchangeable
        # child has their new topology, G, from crossover
        child@GA <- parents[[1]]@GA
        child@GA[g,] <- child@haploint

        # append child to deme
        proposed[[k]] <- append(proposed[[k]], child)
        }
      }
    }

    #.........................
    # iterate through migration
    #.........................
    dem.init.size <- sapply(proposed, length) # note, we track deme size as we do not want to allow an individual to emmigrate twice in one generation. You can imagine this happeneing if an individual from deme 1 emmigrated to deme 5 and then as we were looping through, we gave it a chance to emmigrate again (as K is moving up k++). Instead append it at the end and use the recombo lengths to keep emmigration to (potentially) one per sample one per generation
    for(k in 1:K){ # for each deme
        emigrators <- c() # track who emigrates
        for(p in 1:dem.init.size[k]){ # for each haplotype/chromosome in a given deme
          if(rbinom(1, 1, m)){
            destination <- sample(1:K, size = 1) # find destination
            proposed[[destination]] <- append(proposed[[destination]], proposed[[k]][[p]]) # append haplotype to list
            emigrators <- c(emigrators, p)
          }
        }
        proposed[[k]] <- proposed[[k]][! 1:length(proposed[[k]]) %in% emigrators ] # remove haplotypes that have emigrated
      }; cat("Migration completed for Generation ", g, "\n")


    #.........................
    # UPDATE and store new ancestors
    #.........................
    coi.int <- 1 + rpois(n = K, lambda = coi)

    # partition samples and fill in matrix
    for(k in 1:length(coi.int)){
        current[k] <- list( sample(proposed[[k]], size = coi.int[[k]]) )
    }

  } # end of Generations loop


  return(current)

}

