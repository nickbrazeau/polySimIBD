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
wfsmpl2 <- list.files("~/Desktop/msprimesims/discreteWF/smplsize2/", full.names = T)
wfsmpl2 <- lapply(wfsmpl2, read_msprime_sims) %>%
  dplyr::bind_rows()

wfsmpl3 <- list.files("~/Desktop/msprimesims/discreteWF/smplsize3/", full.names = T)
wfsmpl3 <- lapply(wfsmpl3, read_msprime_sims) %>%
  dplyr::bind_rows()

wfsmpl5 <- list.files("~/Desktop/msprimesims/discreteWF/smplsize5/", full.names = T)
wfsmpl5 <- lapply(wfsmpl5, read_msprime_sims) %>%
  dplyr::bind_rows()


wfsims <- rbind.data.frame(wfsmpl2, wfsmpl3, wfsmpl5) %>%
  dplyr::mutate(model = "discrete wright fisher")


#....................................
# write out
#....................................


saveRDS(wfsims, file = "~/Desktop/msprimesims.rds")

#....................................
# make plots
#....................................


wfsimsplot <- wfsims %>%
  dplyr::filter(time2parent != 0) %>% # terminal branch points
  dplyr::group_by(simnum, treeindex) %>%
  dplyr::arrange(simnum, treeindex, time2parent) %>%
  dplyr::mutate(ET = paste0("T", row_number()),
                ET = factor(ET, levels = paste0("T", 1:100))) %>%
  ggplot() +
  ggridges::geom_density_ridges(aes(x = time2parent, y = ET, fill = ET, color = ET), alpha = 0.6) +
  xlim(0, 1000) +
  facet_wrap(. ~ smplsize) +
  theme_bw() +
  theme(legend.position = "none")

jpeg("~/Desktop/msprime_dtwf_results.jpg", height = 8, width = 8, res = 250, units = "in")
plot(wfsimsplot)
graphics.off()
