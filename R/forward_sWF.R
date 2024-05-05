#' @title The Structured Wright Fisher Model for IBD
#' @param pos vector; the genomic coordinates for chromosome and position of the sites
#' @param N integer vector; The number of individuals to consider in each deme
#' @param m numeric numeric; Probability of internal migration where m represents the probability of moving from host_{origin} to host_{new} by m*(1-1/N) of each deme
#' @param mean_coi numeric vector; The lambda of a right-shifted Poisson process, 1 + Pos(lambda) representing the average COI of each deme
#' @param migr_mat numeric matrix; Migrations rates or probabilities between destination and origin. Note, because this is a Wright-Fisher model, we are drawing parents and therefore migration matrix is parameterized towards "where one came from" versus "where one is headed": origin specified as columns and destination in rows Default value of 1 indicates non-spatial model. Note, if probability matrix, rows must sum to 1 (valid marginal probability); otherwise, values will be assumed to be rates and converted to probabilities
#' @param rho numeric; expected recombination rate
#' @param tlim numeric; the maximum number of generations to consider before exiting gracefully if all samples have not coalesced
#' @param verbose boolean
#' @return Returns a list of length six that contains the COI of each individual. A recombination list of length of tlim
#' where each element contains the recombination block -- as a boolean -- of the two parental haplotypes.   
#'  (the number of generations it took for all lineages to coalesce). Finally, there are lists for the parental host 
#'  and parental haplotype assignments for the "paternal" and "maternal" haplotypes (1 and 2), respectively. 
#' @description Simulate a population forwards with recombination that approximates the Structured Wright Fisher Process
#' and tracks haplotype identity by descent where individuals represent demes, such that within
#' a deme individual-level COI is considered.
#' @details Demes are assumed to be ordered throughout (i.e. the order needs to be consistent between N, m, mean_coi, and the rows and columns of the migration matrix).
#' @details Migration matrix is assumed to be a distance matrix that is either a rate or a probability. The program
#'          will coerce the matrix into a probability distribution between origin and destination based on the row-sums. 
#' @details This function is intended to be fed into the [polySimIBD::get_arg] function for interpretability. 
#'
#' @export

sim_swf <- function(pos, N, m, rho, mean_coi, tlim,
                    migr_mat = 1, 
                    verbose = FALSE){
  
  # assertions
  goodegg::assert_vector(pos)
  goodegg::assert_numeric(pos)
  goodegg::assert_increasing(pos)
  if (length(migr_mat) == 1) {
    goodegg::assert_eq(migr_mat, 1, 
              message = "Distance matrix can be set to 1 to indicate a non-spatial model. Otherwise, needs to be a matrix")
    goodegg::assert_length(N, 1,
              message = "If non-spatial model considered, deme size is of length 1") 
    goodegg::assert_length(m, 1,
              message = "If non-spatial model considered, internal migration rates is of length 1") 
    goodegg::assert_length(mean_coi, 1,
              message = "If non-spatial model considered, mean coi is of length 1") 
  } else {
    goodegg::assert_matrix(migr_mat)
    goodegg::assert_eq(ncol(migr_mat), length(N),
              message = "Migration matrix and deme sizes must be of same relative dimensions. In other words, you need to specify a deme size for every deme")
    goodegg::assert_eq(ncol(migr_mat), length(m),
              message = "Migration matrix and internal migration rates must be of same relative dimensions. In other words, you need to specify a deme size for internal migration")
    goodegg::assert_eq(ncol(migr_mat), length(mean_coi),
              message = "Migration matrix and mean COI must be of same relative dimensions. In other words, you need to specify a mean COI for every deme")
  }
  goodegg::assert_pos_int(N, zero_allowed = FALSE)
  goodegg::assert_bounded(m, left = 0, right = 1)
  goodegg::assert_bounded(rho, left = 0, right = 1, inclusive_left = FALSE, inclusive_right = FALSE)
  goodegg::assert_pos(mean_coi, zero_allowed = FALSE)
  goodegg::assert_single_pos_int(tlim, zero_allowed = FALSE)
  
  
  # precalculate probability of an odd number of recombination events between
  # any interval
  odd_prob <- exp(-rho*diff(pos))*sinh(rho*diff(pos))
  if (any(!is.finite(odd_prob))) {  # catch underflow issue
    odd_prob[is.nan(odd_prob)] <- 0.5
  }
  
  # calculate migration rates and probability from distance matrix
  if (is.matrix(migr_mat)) {
    if(any(rowSums(migr_mat) != 1)) { # if not probability, assume that it is rates
      migr_mat <- 1 - exp(-migr_mat)
      migr_mat <- migr_mat/rowSums(migr_mat)
    } 
    
    # split out for cpp import
    migr_mat <- split(migr_mat, 1:nrow(migr_mat))
    migr_mat <- unname(migr_mat)
  } else {
    migr_mat <- 1
  }

  # define argument list
  args <- list(pos = pos,
               maxN = max(N),
               demecnt = length(N),
               mig_mat_prob = migr_mat,
               N = N,
               m = m,
               mean_coi = mean_coi,
               tlim = tlim,
               odd_prob = odd_prob)
  
  # run efficient C++ function
  ret <- sim_swf_cpp(args)
  
  # keep pos around for ease
  pos <- list(pos = pos)
  ret <- append(pos, ret)  
  # return as custom class
  class(ret) <-"swfsim"
  return(ret)
}




#' @title Get ancestral recombination graph from forward simulations
#' @description Given an object \code{swf}, which is the result of forward
#'   simulation using the function \code{sim_swf()}, walks backwards through the
#'   ancestry and calculates the coalescent tree at every locus for the
#'   specified hosts and/or haplotypes.
#' @param swf result of forwards simulation using the function \code{sim_swf()}
#' @param host_index a vector of target hosts. Defaults to all hosts
#' @param haplo_index a list of target haplotypes within the hosts specified by
#'   \code{host_index}. Defaults to all haplotypes within the specified hosts
#' @importFrom methods new
#' @export

get_arg <- function(swf, host_index = NULL, haplo_index = NULL) {
  
  # check inputs and define defaults
  goodegg::assert_class(swf, "swfsim")
  if (is.null(host_index)) {
    host_index <- 1:length(swf$coi)
  }
  if (is.null(haplo_index)) {
    haplo_index <- mapply(function(x) 1:x, swf$coi[host_index], SIMPLIFY = FALSE)
  }
  goodegg::assert_vector(host_index)
  goodegg::assert_pos_int(host_index, zero_allowed = FALSE)
  goodegg::assert_list(haplo_index)
  goodegg::assert_same_length(host_index, haplo_index)
  goodegg::assert_pos_int(unlist(haplo_index), zero_allowed = FALSE)
  
  # define arguments
  args <- c(swf, list(host_index = rep(host_index, mapply(length, haplo_index)) - 1,
                      haplo_index = unlist(haplo_index) - 1))
  
  # pass to efficient C++ function
  output_raw <- get_arg_cpp(args)
  
  # create a list of bvtrees
  L <- length(output_raw$coalesce_target)
  ARG <- lapply(1:L, function(x) return(new("bvtree")))
  for (l in 1:L) {
    ARG[[l]]@c = output_raw$coalesce_target[[l]]
    ARG[[l]]@t = output_raw$coalesce_time[[l]]
    t_temp <- ARG[[l]]@t
    t_temp[which(t_temp == -1)] <- NA
    ARG[[l]]@z = order(t_temp) - 1
  }
  
  # custom class
  class(ARG) <- "argraph"
  
  return(ARG)
}




#' @title Subset an object of class bvtree
#'
#' @description Given a bvtree and a vector of indices \code{s}, creates a new
#'   tree which is a subset of the original tree focusing only on the elements
#'   \code{s}.
#'
#' @param bvtree an object of class "bvtree"
#' @param s a vector specifying which elements in the bvtree to focus on
#'
#' @export

subset_bvtree <- function(bvtree, s) {
  
  # check inputs
  goodegg::assert_class(bvtree, "bvtree")
  goodegg::assert_vector(s)
  goodegg::assert_gr(length(s), 1)
  goodegg::assert_noduplicates(s)
  goodegg::assert_pos_int(s, zero_allowed = FALSE)
  goodegg::assert_leq(s, length(bvtree@c))
  
  # create mask vector
  m <- rep(0, length(bvtree@c))
  m[s] <- 1
  
  # run efficient C++ function
  output_raw <- subset_bvtree_cpp(bvtree@c, bvtree@t, m)
  
  # update tree
  bvtree@c <- output_raw$c[s]
  bvtree@t <- output_raw$t[s]
  
  # update indices
  bvtree@c <- c(-1, (1:length(s))-1)[match(bvtree@c, c(-1,s-1))]
  
  # update z
  tmp <- bvtree@t
  tmp[tmp == -1] <- NA
  bvtree@z <- order(tmp) - 1
  
  return(bvtree)
}






#' @title Get Between-Host Identity by Descent from forward simulations
#'
#' @description Given an object \code{swf}, which is the result of forward
#'   simulation using the function \code{sim_swf()}, walks backwards through the
#'   ancestry and calculates the between host identity by descent. Calculation is 
#'   based on \cite{Verity et. al 2020, Nat Comms, PMC7192906}.    
#'
#' @inheritParams get_arg
#' @importFrom methods new
#' @export

get_bvibd <- function(swf, host_index = NULL, haplo_index = NULL) {
  
  # check inputs and define defaults
  goodegg::assert_class(swf, "swfsim")
  if (is.null(host_index)) {
    host_index <- 1:length(swf$coi)
  }
  if (is.null(haplo_index)) {
    haplo_index <- mapply(function(x) 1:x, swf$coi[host_index], SIMPLIFY = FALSE)
  }
  goodegg::assert_vector(host_index)
  goodegg::assert_pos_int(host_index, zero_allowed = FALSE)
  goodegg::assert_list(haplo_index)
  goodegg::assert_same_length(host_index, haplo_index)
  goodegg::assert_pos_int(unlist(haplo_index), zero_allowed = FALSE)
  
  # define arguments
  args <- c(swf, list(host_index = rep(host_index, mapply(length, haplo_index)) - 1,
                      haplo_index = unlist(haplo_index) - 1))
  
  # pass to efficient C++ function
  output_raw <- get_bvibd_cpp(args)
  numerator <- output_raw$ibd_target[-1]
  
  # under SNP vs PSMC (Li/Durbin model) don't know begin and end, so treat as missing info - ie burn first loci
  wi <- diff(swf$pos)/sum(diff(swf$pos))
  # weighted average (each loci, denom is 1)
  bv_ibd <- sum( numerator*wi ) 
  
  return(bv_ibd)
}

