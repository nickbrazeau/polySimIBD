
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

simulate_IBD_SWF <- function(chrompos, Nes, Noff, K, m, coi, rho, G){

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


  #.........................
  # setup
  #.........................
  # setup haplos
  Nes.haplos <- lapply(1:Nes, function(x) return( new("simhaplo") ))
  # loop through effective pop
  for(i in 1:length(Nes.haplos)){
    Nes.haplos[[i]]@haploint <- rep(i, nrow(chrompos))
  }
  # setup matrix timeframes

  wf.dem <- lapply(1:(G*K), function(x) return(list()))
  dim(wf.dem) <- c(G, K)

  #.........................
  # initialize
  #.........................
  coi.int <- 1 + rpois(n = K, lambda = coi)
  haplos <- Nes.haplos[sample(x = 1:length(Nes.haplos), size = sum(coi.int), replace = T)]

  # partition samples init
  # and fill in first row (parental row) of WF model
  smpl <- c(0, cumsum(coi.int))
  for(i in 1:length(coi.int)){
    wf.dem[1, i] <- list( haplos[smpl[i]:smpl[i+1]] )
  }; cat("First Generation Population Set \n")


  #..........................................................
  # MAIN
  #..........................................................
  for(g in 2:G){

    #.........................
    # Perform recombination
    #.........................
    recomb.dem <- wf.dem[g-1, ]

    for(k in 1:K){ # for each deme
      # need to initialize list
      parents <- sample(recomb.dem[[k]], size = 2, replace = T)
      child <- makecrossover(p1 = parents[[1]], p2 = parents[[2]], chrompos = chrompos, rho = rho)
      recomb.dem[[k]] <- child

      for(n in 2:Noff){ # make this many recombinants
        parents <- sample(recomb.dem[[k]], size = 2, replace = T)
        child <- makecrossover(p1 = parents[[1]], p2 = parents[[2]], chrompos = chrompos, rho = rho)
        recomb.dem[[k]] <- append(recomb.dem[[k]], child)
      }
    }

    #.........................
    # iterate through migration
    #.........................
    migr.dem <- recomb.dem # note, we make a copy because we don't want to allow an individual to emmigrate twice in one generation. You can imagine this happeneing if an individual from deme 1 emmigrated to deme 5 and then as we were looping through, we gave it a chance to emmigrate again (as K is moving up k++). Instead append it at the end and use the recombo lengths to keep emmigration to (potentially) one per sample one per generation
    for(k in 1:K){ # for each deme
        emigrators <- c() # track who emigrates
        for(p in 1:length(recomb.dem[[k]])){ # for each haplotype/chromosome in a given deme
          if(rbinom(1, 1, m)){
            destination <- sample(1:K, size = 1) # find destination
            migr.dem[[destination]] <- append(migr.dem[[destination]], migr.dem[[k]][[p]]) # append haplotype to list
            emigrators <- c(emigrators, p)
          }
        }
        migr.dem[[k]] <- migr.dem[[k]][! 1:length(migr.dem[[k]]) %in% emigrators ] # remove haplotype that has emigrated
      }
    }; cat("Migration completed for Generation ", g, "\n")


    #.........................
    # UPDATE and store
    #.........................
    coi.int <- 1 + rpois(n = K, lambda = coi)

    # partition samples and fill in matrix
    for(k in 1:length(coi.int)){
        wf.dem[g, k] <- list( sample(migr.dem[[k]], size = coi.int[[k]]) )
    }

  return(wf.dem)

}

