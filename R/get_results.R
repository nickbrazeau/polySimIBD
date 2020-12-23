
#------------------------------------------------
#' internal class
#' @noRd
setClass("bvtree",
         slots=list(c="numeric", t="numeric", z="numeric"))


#' @title Extract Effective COI by Loci from SWF Simulation for a Single Host
#' @inheritParams get_arg
#' @details Only accepts a single host 
#' @return vector of effective COI by loci
#' @export
get_realized_coi <- function(swf, host_index = NULL) {
  
  # checks
  assert_custom_class(swf, "swfsim")
  assert_single_int(host_index)
  
  # get ARG from swf and host_index
  arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
  
  # get connections
  conn <- purrr::map(arg, "c")
  # get effective COI over loci 
  effCOI <- purrr::map_dbl(conn, function(x){sum(x == -1)})
  return(effCOI)
}


#' @title Get Within-Host IBD from SWF Simulation
#' @inheritParams get_arg
#' @details Only accepts a single host 
#' @return double of within-host IBD
#' @export
get_within_host_IBD <- function(swf, host_index = NULL) {
  
  # checks
  assert_custom_class(swf, "swfsim")
  assert_single_int(host_index)
  
  # get ARG from swf and host_index
  arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
  
  # get connections
  conn <- purrr::map(arg, "c")
  # get effective IBD over loci 
  numerator <- purrr::map_dbl(conn, function(x){sum(x != -1)})
  # -1 here for the SELF comparison
  denom <-  (swf$coi[[host_index]]-1) * length(conn)
  # out
  wthnIBD <- sum(numerator)/denom
  return(wthnIBD)
}

#' @title  Effective IBD by Loci from SWF Simulation for a Pair of Hosts
#' @inheritParams get_arg
#' @description Assumes that the minimum realized COI between the pairs of host determines
#' the denominator for the between realized IBD
#' @details  Only accepts a pair of hosts. Ignores mutations as interrupting IBD segments. 
#' @return double of pairwise IBD 
#' @export
get_realized_pairwise_ibd <- function(swf, host_index = NULL) {
  # checks
  assert_custom_class(swf, "swfsim")
  assert_int(host_index)
  assert_length(host_index, 2)
  
  # get ARG from swf and host_index
  arg <- polySimIBD::get_arg(swf = swf, host_index = host_index)
  # take max as it return loci effective COI. Absolute max is 
  # how many strains we could observe (even over few loci)
  coi1 <- polySimIBD::get_realized_coi(swf = swf, host_index = host_index[[1]])
  coi1 <- max(coi1)
  coi2 <- polySimIBD::get_realized_coi(swf = swf, host_index = host_index[[2]])
  coi2 <- max(coi2)
  
  # get connections
  conn <- purrr::map(arg, "c")
  # get timing of connections
  tm <- purrr::map(arg, "t")
  
  #......................
  # get pairwise ibd internal function
  #......................
  get_pairwise_ibd <- function(conni, tmi, this_coi) {
    smpl1con <- conni[1:this_coi[1]]
    smpl2con <- conni[(this_coi[1]+1):(cumsum(this_coi)[2])]
    # get IBD
    # connections between 1 and 2
    pwconn <- which(smpl2con %in% 0:(this_coi[1]-1) )
    locimatches <- rep(1, length(pwconn))
    # note we are 0 based in connections
    # note bvtrees always point left
    # catch if there are multiple matches within sample 2 to the pairwise
    # this is a coalescent tree that looks like below if host COI is 2,2
    # c: -1 -1 1 2
    # t: -1 -1 5 1
    if (length(pwconn) != 0) {
      for (i in 1:length(pwconn)) {
        haplotypeindex <- this_coi[1] + pwconn[i] - 1 # -1 for 0-based
        internalconn <- which(smpl2con %in% haplotypeindex )
        if (length(internalconn) != 0) {
          for (i in 1:length(internalconn)) {
            internalhaplotypeplace <- this_coi[1] + internalconn[i] # here 1-based in R
            if (tmi[internalhaplotypeplace] < tmi[this_coi[1] + internalconn[i]]) { # here 1-based in R
              locimatches[i] <- locimatches[i] + 1
            }
          }
        }
      }
    }
    return(sum(locimatches))
  }
  
  #......................
  # get numerator and denominator 
  # NB liftover such that the effective COI determines the max pairwise relat
  # as well as the denom
  #......................
  effCOI <- c(coi1, coi2)
  numerator <- purrr::map2_dbl(.x = conn, .y = tm,
                               .f = get_pairwise_ibd, this_coi = swf$coi[host_index])
  numerator <- ifelse(numerator > min(effCOI), min(effCOI), numerator)
  pairwiseIBD <- sum(numerator)/(min(effCOI) * length(conn)) # min combn * num Loci
  
  # out
  return(pairwiseIBD)
}


#' Extract haplotypes from ARG
#' @param ARG set of bvtrees
#' @return hapmat numeric matrix; a matrix of mutliallelic haplotypes for each parasite considered. Loci are in
#' rows and parasites (haplotypes) are in columns. 
#' @export
get_haplotype_matrix <- function(ARG){
  
  # convert trees into matrix of alleles
  # each column is therefore a haplotype since we consider parasite by parasite
  hap_mat <- t(mapply(function(x) {
    c <- x@c
    ret <- c
    ret[ret == -1] <- 1:sum(ret == -1)
    while (any(c != -1)) {
      w <- which(c == -1)
      c[-w] <- c[c[-w]+1]
      ret[-w] <- ret[ret[-w]+1]
    }
    return(ret)
  }, ARG))
  return(hap_mat)
}








