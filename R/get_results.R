
#------------------------------------------------
#' internal class
#' @noRd
setClass("bvtree",
         slots=list(c="numeric", t="numeric", z="numeric"))

#------------------------------------------------
#' @title Extract Effective COI by Loci from SWF Simulation for a Single Host
#' @inheritParams get_arg
#' @details Only accepts a single host 
#' @description TODO 
#' @return vector of effective COI by loci
#' @export
get_effective_coi <- function(swf, host_index = NULL) {
  
  # checks
  assert_custom_class(swf, "swfsim")
  assert_single_pos_int(host_index)
  
  # get ARG from swf and host_index
  arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
  
  # get connections
  conn <- purrr::map(arg, "c")
  # get effective COI over loci 
  effCOI <- purrr::map_dbl(conn, function(x){sum(x == -1)})
  return(effCOI)
}

#------------------------------------------------
#' @title Calculate Within-Host IBD 
#' @inheritParams get_arg
#' @details Only accepts a single host 
#' @description TODO
#' @return double of within-host IBD
#' @export
get_within_ibd <- function(swf, host_index = NULL) {
  # checks
  assert_custom_class(swf, "swfsim")
  assert_single_pos_int(host_index)
  # get ARG from swf and host_index
  # need to call ARG again to make extraction straightforward (eg don't know what user did to ARG upstream)
  arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
  # get connections
  conn <- purrr::map(arg, "c")
  # get within host IBD per loci 
  numerator <- sapply(conn, function(x){sum(x != -1)})
  # -1 here for self comparison
  denom <-  (swf$coi[[host_index]]-1)
  # under SNP vs PSMC (Li/Durbin model) don't know begin and end, so treat as missing info
  wi <- diff(swf$pos)/sum(diff(swf$pos))
  # out
  wthnIBD <-sum( (numerator[1:(length(numerator) - 1)] / denom) * wi )
  return(wthnIBD)
}



#------------------------------------------------
#' @title Get Connection Intervals 
#' @description Index where unique connections are in the entire genome for proper weighting
#' @param uniqueconn unique bvtree connections from the ARG
#' @param allconn all bvtree connections from the ARG
#' @noRd
#' @noMd

get_conn_intervals <- function(uniqueconn, allconn){
  names(uniqueconn) <- 1:length(uniqueconn)
  intervals <- lapply(uniqueconn,
                      function(uni){
                        return(sapply(allconn, function(x){paste(uni, collapse = "") == paste(x, collapse = "")}))})
  mint <- rep(NA, length(allconn))
  for(i in 1:length(intervals)) {
    mint[intervals[[i]]] <- names(intervals)[i]
  }
  return(as.numeric(mint))
}

#------------------------------------------------
#' @title Calculate Between Host (pairwise) IBD 
#' @inheritParams bvtree
#' @param int internal identification of root indices 
#' @description Internal function: Recursively loop through tree starting with root to find haplo-indices in the "bvtrees c slot" that are connected
#' @noMd
#' @noRd
get_conn_from_root <- function(root, c) {
  # init
  tree <- newroot <- root
  #   
  while(length(newroot) != 0) {
    newroot <- which(c %in% (newroot-1)) # match c value and update new root 
    tree <- c(tree, newroot)
  }
  # out
  return(tree)
}

#------------------------------------------------
#' @title Find relevant connections for within IBD calculation
#' @inheritParams bvtree
#' @description Internal function: subset bvtree "c" slot to haplo-indices that contain a between sample connections and therefore contribute to within-sample IBD in a pairwise comparison 
#' @noMd
#' @noRd

get_withinIBD_bvtree_subset <- function(c, coi1, coi2) {
  # find roots
  roots <- which(c == -1)
  # subset to tree based on roots 
  subset_trees <- lapply(roots, get_bvtree_from_root, c = c)
  # get btwn conn
  btwnconn <- which(c[(coi1+1):(coi2+coi1)] %in% 0:(coi1-1)) + coi1
  
  # drop spmls w/ no btwn 
  subset_trees <- subset_trees[ sapply(subset_trees, function(x, btwn){any(x %in% btwn)}, btwn = btwnconn) ]
  # out
  return(sort(unique(unlist(subset_trees))))
}

#------------------------------------------------
#' @title Calculate Between Host (pairwise) IBD 
#' @description Given an object \code{swf}, ***
#' @inheritParams get_arg
#' @description ***
#' @details  Only accepts a pair of hosts (i.e. pairwise). Ignores mutations as interrupting IBD segments. 
#' @return ***
#' @export
get_pairwise_ibd <- function(swf, host_index = NULL) {
  # check inputs and define defaults
  assert_custom_class(swf, "swfsim")
  assert_vector(host_index)
  assert_noduplicates(host_index)
  assert_pos_int(host_index, zero_allowed = FALSE)
  if(length(host_index) != 2) {
    stop("host_index must be of length 2 for pairwise comparison", call. = FALSE)
  }
  # get ARG from swf and host_index
  # need to call ARG again to make extraction straightforward (eg don't know what user did to ARG upstream)
  arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
  # subset to unique loci for speed 
  conn <- purrr::map(arg, "c")
  uniconn <- unique(conn)
  
  # define arguments for fast cpp function 
  argums <- list(conn = uniconn, 
                 host_haplo_cnt = swf$coi[host_index])
  
  # we define btwn_pairwise_ibd as: 
  # (n_{coal-btwn} + n_{coal-win-btwn}) / (n_{strains1} * n{strains2})
  # where w/in coal share a btwn coal 
  
  # pass to efficient C++ function for quick between tree look up 
  # to determine between host IBD
  output_raw <- calc_between_IBD_cpp(argums)$ibd_numerator
  
  #tidy raw
  # catch if no btwn, no w/in or extra work needed
  if(sum(output_raw == 0)) { return(0)}
  # if within, need to do additional work 
  
  
  # find the locations of unique loci locations for later expansion
  conn_indices <- polySimIBD:::get_conn_intervals(uniqueconn = uniconn, allconn = conn)
  
  # subset to relevant haploindices in bvtrees for correct w/in IBD calculation 
  # i.e. only within IBD that has a pairwise connection, or btwn smpl connection, contributes to 
  # overall calculation of pairwise IBD 
  relconn <- lapply(uniconn, polySimIBD:::get_withinIBD_bvtree_subset, 
                    coi1 = swf$coi[host_index][1],
                    coi2 = swf$coi[host_index][2])
  wunniconn <- mapply(function(x, y){return(x[y])}, x = uniconn, y = relconn, SIMPLIFY = F)
  
  # need to still respect COI, which is [ ] below 
  # within host IBD for smpl 1
  win_smpl1 <- mapply(x = wunniconn, y = relconn, function(x,y){
    z <- x[ y %in% 1:swf$coi[host_index][1] ]
    return(sum(z != -1)) }, SIMPLIFY = T)
  # within host IBD for smpl 1
  win_smpl2 <- mapply(x = wunniconn, y = relconn, function(x,y){
    z <- x[ y %in% swf$coi[host_index][1]:sum(swf$coi[host_index]) ]
    return(sum(z != -1)) }, SIMPLIFY = T)
  
  # expand out unique loci intervals from above
  outputraw <- outputraw[conn_indices]
  win_smpl1 <- win_smpl1[conn_indices]
  win_smplw <- win_smplw[conn_indices]
  
  # numerator of btwn and w/in IBD
  numerator <- output_raw + win_smpl1 + win_smpl2
  
  # under SNP vs PSMC (Li/Durbin model) don't know begin and end, so treat as missing info
  wi <- diff(swf$pos)/sum(diff(swf$pos))
  # weighted average
  return( sum( (numerator / prod(swf$coi[host_index]))*wi ) )
}



#------------------------------------------------
#' Extract haplotypes from ARG
#' @param arg set of bvtrees
#' @return hapmat numeric matrix; a matrix of mutliallelic haplotypes for each parasite considered. Loci are in
#' rows and parasites (haplotypes) are in columns. 
#' @export
get_haplotype_matrix <- function(arg){
  
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
