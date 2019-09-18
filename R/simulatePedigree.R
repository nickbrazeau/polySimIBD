#' @title Simulate Pedigree IBD
#' @param chrompos matrix <factor><numeric>; the genomic coordinates for chromosome and position of the sites
#' @param PLAF numeric vector; the population-level allele frequencies to simulate genotypes from
#' @param rho numeric; expected recombination rate
#' @param k numeric; number of generation to simulate forward
#' @param inbreeding numeric; The probability of inbreeding for each generation. Backcrosses are made from the n-1 generations. The probability is bounded by 0-1.
#'
#'@details Data is simulated under the following framework:
#'   \enumerate{
#'        \item Parental genotypes are drawn from the \code{PLAF}
#'        \item F1 crosses are made from parental genotypes
#'        \item Outgroup genotype that are unique to the parental and F1 progeny are made from the \code{PLAF}
#'        \item Encode haplotypes with a \emph{bit}
#'        \item Crossover F-progeny with outgroups unless inbreeding is specified. If the probability of inbreeding within each generation is specified,
#'        then either a relative or outgroup is selected following a binomial distribution: \eqn{binom(1, p_inbreeding)} for that generation's mate.
#'   }
#'
#' @export

simulate_IBD_pop_Pedigree <- function(chrompos, rho, k,
                                      inbreeding = NA){

  # assert that inbreeding is not =<0 or 1=<
  #
  #..........................
  # Simulate Parents in Pedigree
  #..........................
  p1 <- new("simhaplo")
  p2 <- new("simhaplo")
  p1@haploint <- rep(1, nrow(chrompos))
  p2@haploint <- rep(2, nrow(chrompos))

  #..........................
  # Simulate F1 progeny
  #..........................
  f1.1 <- makecrossover(p1 = p1, p2 = p2, chrompos = chrompos, rho = rho)[[1]] # just grab first child
  f1.2 <- makecrossover(p1 = p1, p2 = p2, chrompos = chrompos, rho = rho)[[1]] # just grab first child


  #..........................
  # Simulate Mixing of Outgroups only
  # for each F1 progeny
  #..........................
  og <- lapply(1:((k-1)*2), function(x) return(new("simhaplo")) )
  for(i in 1:length(og)){
    og[[i]]@haploint <- rep( i+2 , nrow(chrompos) ) # 1,2 taken by parents
  }

  #..........................
  # Create potential outgroup matings
  #..........................
  f1matings <- og[1:(length(og)/2)]
  f2matings <- og[(length(og)/2 + 1):length(og)]

  f1.1.ks <- list(p1, f1.1)
  f1.2.ks <- list(p1, f1.2)

  #..........................
  # Simualte Matings
  #..........................
  # simulate through lineage 1
  for(i in 1:length(f1matings)){
    if(!is.na(inbreeding)){ # if user specified inbreeding
      if(runif(1) < inbreeding){ # flip weighted coin for whether it is inbred or outgroup
        mate1 <- f1.1.ks[[ sample(1:(length(f1.1.ks)-1), 1) ]] # -1 so no selfings
      } else{
        mate1 <- f1matings[[i]]
      }
    } else { # no user specified inbreeding
      mate1 <- f1matings[[i]]
    }

    f1.1 <- makecrossover(p1 = f1.1, p2 = mate1, chrompos = chrompos, rho = rho)[[1]] # just grab fist child
    f1.1.ks <- append(f1.1.ks, f1.1)
  } #end for loop f1

  # simulate through lineage 2
  for(i in 1:length(f2matings)){
    if(!is.na(inbreeding)){ # if user specified inbreeding
      if(runif(1) < inbreeding){ # flip weighted coin for whether it is inbred or outgroup
        mate2 <- f1.2.ks[[ sample(1:(length(f1.2.ks)-1), 1) ]] # -1 so no selfings
      } else{
        mate2 <- f2matings[[i]]
      }
    } else { # no user specified inbreeding
      mate2 <- f2matings[[i]]
    }

    f1.2 <- makecrossover(p1 = f1.2, p2 = mate2, chrompos = chrompos, rho = rho)[[1]] # just grab first child
    f1.2.ks <- append(f1.2.ks, f1.2)
  } #end for loop f2



  #............................
  # returns
  #............................

  ret <- list(
    k = k,
    f1.1lineage = f1.1.ks,
    f1.2lineage = f1.2.ks,
    kprogeny1 = f1.1,
    kprogeny2 = f1.2
    )

  class(ret) <- "simulatedpedigree"

  return(ret)



}

