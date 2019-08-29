
#' @title Simulate IBD-WF
#'
#' @desciption Simulate a population with recombination that approximates the Wright Fisher Process
#' and tracks haplotype identity by descent
#'
#' @param chrompos matrix <factor><numeric>; the genomic coordinates for chromosome and position of the sites
#' @param PLAF numeric vector; the population-level allele frequencies to simulate genotypes from
#' @param Nes numeric; The starting effective population size
#' @param rho numeric; expected recombination rate
#' @param k numeric; number of generation to simulate forward
#'
#'@details Data is simulated under the following framework:
#'   \enumerate{
#'        \item Haploid genotypes are drawn from the \code{PLAF}
#'        \item Encode haplotypes with a \emph{bit}
#'        \item Simulate forward in time from the inital effective population \code{Nes} by drawing two samples
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

simulate_IBD_pop_WF <- function(chrompos, PLAF, Nes, rho, k){


  Nes <- lapply(1:Nes, function(x) return( new("simhaplo") ))
  # loop through first generation based on WF model
  for(i in 1:length(Nes)){
    ret <- c()
    for(j in 1:nrow(chrompos)){
      ret <- append( ret, sample(x=c(0,1), size=1, prob = c(PLAF[j], 1-PLAF[j])) )
    }
    Nes[[i]]@haplogt <- ret
  }

  #.....................................................................................
  # SET UP Initial Haplotype at Generation 1
  # need to account for any first generation duplicated haplotypes
  # by applying the same haplotype bit
  #.....................................................................................
  haplo <- lapply(Nes, function(x){paste0(x@haplogt, collapse = "")})
  if(any(duplicated(haplo))){

    # find offenders
    prws <- as.matrix( t( combn(c(1:length(haplo)), 2) ) )
    prs <- mapply(function(x,y){identical(haplo[[x]], haplo[[y]])}, prws[,1], prws[,2])
    prws <- cbind.data.frame(prws, prs) %>%
      magrittr::set_colnames(c("pair1", "pair2", "ident")) %>%
      dplyr::filter(ident == TRUE) %>%
      dplyr::select(-c("ident"))
    conj <- dplyr::intersect(prws$pair1, prws$pair2) # have conjugancy
    prws <- prws %>%
      dplyr::filter(! pair1 %in% conj)

    duplvls <- length(levels(factor(prws$pair1)))
    prws$haplobit <- factor(prws$pair1, labels = letters[1:duplvls])
    prws <- tidyr::gather(prws, key = "pair", value = "haplonum", 1:2) %>%
      dplyr::select(-c("pair")) %>%
      dplyr::filter(!duplicated(haplonum))

    # now that we found duplicated offenders, find unique haplotypes
    uniqhaps <- data.frame(haplonum = which(! 1:length(Nes) %in% prws$haplonum ),
                           haplobit = NA)
    uniqhaps$haplobit <- letters[(duplvls+1):(nrow(uniqhaps)+(duplvls))]


    # apply haplotypes
    hp <- rbind.data.frame(prws, uniqhaps)

    for(i in 1:length(Nes)){
      Nes[[i]]@haplobit <- rep( as.character( hp$haplobit[ hp$haplonum == i ] ), length(sites) )
    }

  } else { # no duplicated haplotypes in initial draw

    for(i in 1:length(Nes)){
      Nes[[i]]@haplobit <- rep( letters[i], length(sites) )
    }
  }

  #.....................................................................................
  # Simulate Forward in Time
  #.....................................................................................

  # TODO find work around for this nested s4 and this suppress warning issue
  nedf <- suppressWarnings( tibble::as_tibble(matrix(list(), k, length(Nes))) )
  suppressWarnings(nedf[1,] <- Nes)

  for(i in 2:nrow(nedf)){
    for(j in 1:ncol(nedf)){
      # pick two samples from previous generation
      # note, we can pick the same sample, which would represent gen identical
      # if we pick two different samples we have gen distinct
      s1 <- nedf[i-1, sample(x = 1:ncol(nedf), size = 1)][[1]][[1]] # R's nested S4
      s2 <- nedf[i-1, sample(x = 1:ncol(nedf), size = 1)][[1]][[1]] # R's nested S4

      #..............................
      # Make Recombination Block
      # (sporozoite stage)
      #..............................
      recombo.block <- findrecombination(chrompos = chrompos)

      #..............................
      # Recombination
      #..............................
      child <- makecrossover(p1 = s1, p2 = s2, rho = rho, chrompos = chrompos)
      nedf[i,j] <- child


    }} # end nested-for loop


  return(nedf)

}












