test_that("slower R function matches Cpp function for IBD", {
  
  #......................
  # sim framework
  #......................
  
  
  
  #......................
  # old function
  #......................  
  old_get_realized_pairwise_ibd <- function(swf, host_index = NULL) {
    # checks
    assert_custom_class(swf, "swfsim")
    assert_pos_int(host_index)
    assert_length(host_index, 2)
    
    # get ARG from swf and host_index
    arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
    # take max as it return loci effective COI. Absolute max is 
    # how many strains we could observe (even over few loci)
    coi1 <- polySimIBD::get_realized_coi(swf = swf, host_index = host_index[[1]]) # [[1]] here for first sample, accessing vector position from input above
    coi1 <- max(coi1)
    coi2 <- polySimIBD::get_realized_coi(swf = swf, host_index = host_index[[2]]) # per above [[2]]
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
  
  
})
