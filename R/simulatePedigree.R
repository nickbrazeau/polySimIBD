

#' @title Simulate Pedigree IBD
#' @param pos numeric vector; the genomic coordinates for chromosome and position of the sites
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

simulate_IBD_pop_Pedigree <- function(pos, rho, k,
                                      inbreeding = 0){


  #..........................
  # Assertions
  #..........................
  #assert_matrix(chrompos)
  assert_numeric(pos)
  assert_numeric(rho)
  assert_numeric(k)
  assert_single_bounded(inbreeding, left = 0, right = 1, inclusive_left = T, inclusive_right = T)


  #..........................
  # Simulate Parents in Pedigree
  #..........................
  p1 <- rep(1, length(pos))
  p2 <- rep(2, length(pos))

  #..........................
  # Simulate F1 progeny
  #..........................
  f1.1 <- recombine(p1 = unique(p1), p2 = unique(p2), pos = pos, rho = rho)
  f1.2 <- recombine(p1 = unique(p1), p2 = unique(p2), pos = pos, rho = rho)


  #..........................
  # Simulate Mixing of Outgroups only
  # for each F1 progeny
  #..........................
  og <- 1:((k-1)*2) + 2 # for parent index

  #..........................
  # Create potential outgroup matings
  #..........................
  f1matings <- og[1:(length(og)/2)]
  f2matings <- og[(length(og)/2 + 1):length(og)]

  f1.1.ks <- list(p1, f1.1)
  f1.2.ks <- list(p2, f1.2)

  #..........................
  # Simualte Matings
  #..........................
  # simulate through lineage 1
  for(i in 1:length(f1matings)){

    if(runif(1) < inbreeding){ # flip weighted coin for whether it is inbred or outgroup
      mate1 <- f1.1.ks[[ sample(1:(length(f1.1.ks)-1), 1) ]] # -1 so no selfings
    } else{
      mate1 <- f1matings[[i]]
    }


    f1.1_new <- recombine(p1 = 1, p2 = mate1, pos = pos, rho = rho)
    # get haploint back
    f1.1 <- ifelse(f1.1_new == 1, f1.1, mate1)
    f1.1.ks <- append(f1.1.ks, list(f1.1))
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

    f1.2_new <- recombine(p1 = 2, p2 = mate2, pos = pos, rho = rho)
    # get haploint back
    f1.2 <- ifelse(f1.2_new == 2, f1.2, mate2)
    f1.2.ks <- append(f1.2.ks, list(f1.2))
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

