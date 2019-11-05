library(tidyverse)

simfiles <- list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/polySimIBD/_rslurm_edasim_Forward_sWF_runs/",
                       ".RDS", full.names = T)
simfiles <- simfiles[!grepl("f.RDS|params.RDS", simfiles)]
params <- readRDS("~/Documents/MountPoints/mountedScratchLL/Projects/polySimIBD/_rslurm_edasim_Forward_sWF_runs/params.RDS")

# fix order
simfiles <- simfiles %>%
  as_tibble(.) %>%
  magrittr::set_colnames("path") %>%
  dplyr::mutate(ord = basename(path),
                ord = stringr::str_extract(ord, "[0-9]+"),
                ord = as.numeric(ord)) %>%
  dplyr::arrange(., ord)

#......................................
# Now liftover and get to coal
#......................................

readSim_coaltimes <- function(path){
  ret <- readRDS(path) # read in paths that have multiple runs of ARGsim
  numsims <- length(ret) # get number of sims
  locinum <- length(ret[[1]]$ARG) # get loci count for split
  k <- factor( sort( rep(seq(1:numsims), locinum) ) )

  ret <- purrr::map(ret, "ARG") # pull out arg
  ret <-  unlist(ret, recursive=FALSE) # take one level up
  ret <- lapply(ret, function(x) return(x@t)) # get t

  ret <- split(ret, k) # split
  ret <- lapply(ret, function(x){ # bind times
    x <- as_tibble( t( rbind.data.frame(x) ) ) # tibble doesn't drop to vector when column is 1 like dataframe
    x <- x[, !apply(x, 2, function(z){all(z==-1)})] # drop root
    # sort rows
    if(ncol(x) > 1){
     x <- as_tibble( t(apply(x, 1, function(y){return(sort(y, decreasing = F))})) )
    }

    colnames(x) <- paste0("t", seq(1:ncol(x)))
    rownames(x) <- NULL
    return(x)
    })
  return(ret)
}

simtimes <- lapply(simfiles$path, readSim_coaltimes)
#......................................
# expand for chrom pos
#......................................
params$sims <- unlist(simtimes, recursive=FALSE)


locinum <- nrow(params$sims[[1]])
# expand out and find intervals
maxpos <- max(params$pos[[1]])
pos <- seq(from = 0, to = maxpos)
intvl <- tibble(intvl = findInterval(pos, params$pos[[1]]))

params$simsposexp <- lapply(params$sims, function(x, intvl){
  locinum <- seq(from = 1, to = nrow(x))
  ret <- cbind.data.frame(intvl = locinum, x)
  ret <- dplyr::left_join(intvl, ret, by = "intvl")
  return(ret)
}, intvl = intvl)


#......................................
# save out
#......................................
params <- params %>%
  dplyr::select(c("sample_size", "simsposexp")) %>%
  tidyr::unnest(cols = simsposexp)

saveRDS(params, file = "simulations/msprimesims/simdata/polySimIBD_sims.RDS")




