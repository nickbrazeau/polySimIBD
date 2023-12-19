
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
  goodegg::assert_class(swf, "swfsim")
  goodegg::assert_single_pos_int(host_index)
  
  # get ARG from swf and host_index
  arg <- polySimIBD::get_arg(swf = swf, host_index = host_index)
  
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
  goodegg::assert_class(swf, "swfsim")
  goodegg::assert_single_pos_int(host_index)
  # get ARG from swf and host_index
  # need to call ARG again to make extraction straightforward (eg don't know what user did to ARG upstream)
  arg <- polySimIBD::get_arg(swf = swf, host_index = host_index)
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
#' @title Identify sub-trees recursively from root 
#' @param root internal identification of root indices 
#' @param c connections from bvtree slot
#' @description Internal function: Recursively loop through tree starting with root to find haplo-indices in the "bvtrees c slot" that are connected
#' @details Internal function, not intended for general use
#' @export
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
  # subset to subtree based on roots (ie extract out tree that is based connected to root)
  subset_trees <- lapply(roots, get_conn_from_root, c = c)
  # get btwn conn
  btwnconn <- which(c[(coi1+1):(coi2+coi1)] %in% 0:(coi1-1)) + coi1
  
  # drop spmls w/ no btwn 
  subset_trees <- subset_trees[ sapply(subset_trees, function(x, btwn){any(x %in% btwn)}, btwn = btwnconn) ]
  # out
  return(sort(unique(unlist(subset_trees))))
}

