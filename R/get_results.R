
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
#' @title Calculate Between Host (pairwise) IBD by Interhost Relatedness
#' @inheritParams get_arg
#' @description Between host, or pairwise, IBD as presented in 
#' Verity et al. 2020 (erity realized ibd (PMC7192906). The calculation 
#' ignores intra-host relatedness and does not account for COI. Instead, 
#' if there is any recent coalescence between any of the strains between 
#' host 1 and host 2, the locus is condered to be IBD
#' @return double of between-host IBD
#' @export
get_pairwise_bv_ibd <- function(swf, host_index = NULL) {
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
  # extract connectins 
  conn <- purrr::map(arg, "c")
  # subset to unique loci for speed 
  uniconn <- unique(conn)
  # find the locations of unique loci locations for later expansion
  conn_indices <- polySimIBD:::get_conn_intervals(uniqueconn = uniconn, allconn = conn)
  
  # get between connections per locus for pair 
  # bvtrees point left --  ignring w/in ibd - if any between ibd, then locus is ibd 
  get_loci_pairwise_ibd <- function(conni, this_coi) {
    smpl1con <- conni[1:this_coi[1]]
    smpl2con <- conni[(this_coi[1]+1):(cumsum(this_coi)[2])]
    # get IBD connections between 1 and 2
    pwconn <- which(smpl2con %in% 0:(this_coi[1]-1) )
    # if any, IBD
    pwconn <- sum(pwconn)
    return(pwconn >= 1)
  }
  
  # iterate through/apply function
  lociIBD <- purrr::map_dbl(.x = uniconn,
                            .f = get_loci_pairwise_ibd, this_coi = swf$coi[host_index])
  # expand out unique loci intervals from above
  numerator <- lociIBD[conn_indices]
  
  # under SNP vs PSMC (Li/Durbin model) don't know begin and end, so treat as missing info - ie burn first loci
  numerator <- numerator[-1]
  wi <- diff(swf$pos)/sum(diff(swf$pos))
  # weighted average (each loci, denom is 1)
  return( sum( numerator*wi ) )
  
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
#' @title Identify sub-trees recursively from root 
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
  # subset to subtree based on roots (ie extract out tree that is based connected to root)
  subset_trees <- lapply(roots, polySimIBD:::get_conn_from_root, c = c)
  # get btwn conn
  btwnconn <- which(c[(coi1+1):(coi2+coi1)] %in% 0:(coi1-1)) + coi1
  
  # drop spmls w/ no btwn 
  subset_trees <- subset_trees[ sapply(subset_trees, function(x, btwn){any(x %in% btwn)}, btwn = btwnconn) ]
  # out
  return(sort(unique(unlist(subset_trees))))
}


#------------------------------------------------
#' @title Calculate Between Host (pairwise) IBD accounting for COI
#' @description Given an object \code{swf}, ***
#' @inheritParams get_arg
#' @description ***
#' @details  Only accepts a pair of hosts (i.e. pairwise). Ignores mutations as interrupting IBD segments. 
#' @return ***
#' @export
get_pairwise_coi_ibd <- function(swf, host_index = NULL) {
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
  output_raw <- calc_between_coi_IBD_cpp(argums)$ibd_numerator
  
  #tidy raw
  # catch if no btwn, no w/in or extra work needed
  if(sum(output_raw) == 0) { return(0)}
  # if within, need to do additional work 
  
  
  # find the locations of unique loci locations for later expansion
  conn_indices <- polySimIBD:::get_conn_intervals(uniqueconn = uniconn, allconn = conn)
  
  # subset to relevant haploindices in bvtrees for correct w/in IBD calculation 
  # i.e. only within IBD that has a pairwise connection, or btwn smpl connection, contributes to 
  # overall calculation of pairwise IBD 
  haploind <- lapply(uniconn, polySimIBD:::get_withinIBD_bvtree_subset, 
                     coi1 = swf$coi[host_index][1],
                     coi2 = swf$coi[host_index][2])
  # subset to relevant haploindices for sample 1 based on original COI 
  win_smpl1 <- mapply(function(x, y){
    smpl1_haplo_ind_rel <- which(y %in% 1:swf$coi[host_index][1])
    return(x[smpl1_haplo_ind_rel])
  }, 
  x = uniconn, 
  y = haploind, 
  SIMPLIFY = F)
  # loop through now for w/in calc smpl1 
  win_smpl1 <- sapply(win_smpl1, function(x){
    sum(x %in% 0:(swf$coi[host_index][1]-1))
  })
  
  # subset to relevant haploindices for sample 2 based on original COI 
  win_smpl2 <- mapply(function(x, y){
    smpl2_haplo_ind_rel <- which(y %in% (swf$coi[host_index][1]+1):sum(swf$coi[host_index]))
    return(x[smpl2_haplo_ind_rel])
  }, 
  x = uniconn, 
  y = haploind, 
  SIMPLIFY = F)
  # loop through now for w/in calc smpl2 
  win_smpl2 <- sapply(win_smpl2, function(x){
    sum(x %in% swf$coi[host_index][1]:(sum(swf$coi[host_index])-1))
  })
  
  # expand out unique loci intervals from above
  output_raw <- output_raw[conn_indices]
  win_smpl1 <- win_smpl1[conn_indices]
  win_smpl2 <- win_smpl2[conn_indices]
  
  # numerator of btwn and w/in IBD
  numerator <- output_raw + win_smpl1 + win_smpl2
  
  # under SNP vs PSMC (Li/Durbin model) don't know begin and end, so treat as missing info - ie burn first loci
  numerator <- numerator[-1]
  wi <- diff(swf$pos)/sum(diff(swf$pos))
  # weighted average
  return( sum( (numerator / prod(swf$coi[host_index]))*wi ) )
}



#------------------------------------------------
#' Extract haplotypes from ARG
#' @param arg set of bvtrees
#' @return hapmat numeric matrix; a matrix of multiallelic haplotypes for each parasite considered. Loci are in
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
