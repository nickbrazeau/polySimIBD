library(tidyverse)

read_msprime_sims <- function(path){
  smplsize <- stringr::str_extract(path, "smplsize[0-9]")
  smplsize <- stringr::str_extract(smplsize, "[0-9]")
  simnum <- stringr::str_extract(basename(path), "[0-9]+")


  dat <- readr::read_tsv(path, col_names = F)
  dat$simnum <- simnum
  colnames(dat) <- c("treeindex", "interval", "node", "nodeparent", "time2parent", "simnum")

  dat <- data.frame(smplsize = smplsize, dat)
  return(dat)
}


#....................................
# discrete wf from msprime
#....................................
wfsmpl2 <- list.files("simulations/msprimesims/simdata/msprimesims/discreteWF/smplsize2/", full.names = T)
wfsmpl2 <- lapply(wfsmpl2, read_msprime_sims) %>%
  dplyr::bind_rows()

wfsmpl3 <- list.files("simulations/msprimesims/simdata/msprimesims/discreteWF/smplsize3/", full.names = T)
wfsmpl3 <- lapply(wfsmpl3, read_msprime_sims) %>%
  dplyr::bind_rows()

wfsmpl5 <- list.files("simulations/msprimesims/simdata/msprimesims/discreteWF/smplsize5/", full.names = T)
wfsmpl5 <- lapply(wfsmpl5, read_msprime_sims) %>%
  dplyr::bind_rows()


wfsims <- rbind.data.frame(wfsmpl2, wfsmpl3, wfsmpl5) %>%
  dplyr::mutate(model = "discrete wright fisher")


#....................................
# cleanup and process raw msprime style dataframe
#....................................

process_msprime <- function(dat_raw) {

  # fix data types
  dat_raw$simnum <- as.numeric(dat_raw$simnum)
  dat_raw$smplsize <- as.numeric(as.character(dat_raw$smplsize))

  # extract interval start and end
  interval_start_end <- mapply(function(x) {
    as.numeric(strsplit(x, split = "\\(|\\)|\\,")[[1]][2:3])
  }, dat_raw$interval)

  # append interval to raw data
  dat_raw$interval0 <- interval_start_end[1,]
  dat_raw$interval1 <- interval_start_end[2,]

  # return
  return(dat_raw)
}


# extract coalescent times at every locus from msprime style dataframe
msprime_to_coal <- function(dat_raw, pos) {

  # extract sample size (n)
  n <- dat_raw$smplsize[1]

  # extract coalescent times of all sims
  ret <- mapply(function(i) {

    # subset to this sim
    dat_sub_sim <- subset(dat_raw, simnum == i)

    # extract times for every interval in this sim
    t_mat <- t(mapply(function(j) {

      # subset to this tree
      dat_sub_tree <- subset(dat_sub_sim, treeindex == j)

      # extract times
      w <- which(!dat_sub_tree$node %in% 0:(n-1))
      internal <- dat_sub_tree$node[w]
      m <- match(dat_sub_tree$nodeparent, internal)
      t <- sort(rep(dat_sub_tree$time2parent[w], times = tabulate(m) - 1))

      # append interval information and return
      ret <- c(interval0 = dat_sub_tree$interval0[1],
               t = t)
      return(ret)

    }, unique(dat_sub_sim$treeindex)))

    # convert to coalescent time per locus
    ret <- t_mat[findInterval(pos, t_mat[,1]), -1, drop = FALSE]

    return(ret)

  }, unique(dat_raw$simnum), SIMPLIFY = FALSE)

  # return list
  return(ret)
}


# loop through parameter combinations and get coalescent times into nested list
get_output <- function(dat_raw, pos, model_vec, smplsize_vec) {

  # loop through parameter combinations
  out_list <- list()
  for (i in 1:length(model_vec)) {
    out_list[[i]] <- list()
    for (j in 1:length(smplsize_vec)) {

      # subset data
      dat_sub <- subset(dat_raw, model == model_vec[i] & smplsize == smplsize_vec[j])

      # extract coalescent times
      coal_list <- msprime_to_coal(dat_raw = dat_sub,
                                   pos = pos)

      # save to list
      out_list[[i]][[j]] <- coal_list

    }
    names(out_list[[i]]) <- smplsize_vec
  }
  names(out_list) <- model_vec

  # return
  return(out_list)
}

# cleanup and process
msprocessed <- process_msprime(wfsims)
# define parameters
pos <- seq(0, 1000, 1)
model_vec <- c("discrete wright fisher")
smplsize_vec <- c(2,3,5)

# loop through parameter combinations and get coalescent times
msprimeoutput <- get_output(dat_raw = msprocessed,
                            pos = pos,
                            model_vec = model_vec,
                            smplsize_vec = smplsize_vec)



#....................................
# write out
#....................................
saveRDS(msprimeoutput, file = "simulations/msprimesims/simdata/msprimesims.rds")

