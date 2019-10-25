library(tidyverse)


#....................................
# Read in msprimes and manipulate
#....................................
msprimesims <- readRDS("msprimesims/simdata/msprimesims.rds")
msprimesims <- unlist(msprimesims, recursive = F)

msprimesims$`discrete wright fisher.2` <- # overwrite t1 name for easier bind
  lapply(msprimesims$`discrete wright fisher.2`, function(x){
  colnames(x) <- c("t1")
  return(x)
})

msprimesims.df <- tibble(sample_size = c(2,3,5),
                         simsposexp = msprimesims)
msprimesims.df <- msprimesims.df %>%
  tidyr::unnest(cols = simsposexp)

# this craziness to get R purrr to cooperate
msprimesims.df$simsposexp <- purrr::map(msprimesims.df$simsposexp,
                                        tibble::as_tibble)

msprimesims.df <- msprimesims.df %>%
  tidyr::unnest(cols = simsposexp) %>%
  dplyr::mutate(model = "msprime")

#....................................
# Read in polysimibd
#....................................
polySimIBDsims.df <- readRDS("msprimesims/simdata/polySimIBD_sims.RDS") %>%
  dplyr::select(-c("intvl")) %>%
  dplyr::mutate(model = "polySimIBD")


#....................................
# Comparisons
#....................................
simdf <- rbind.data.frame(msprimesims.df, polySimIBDsims.df)

simdf <- simdf %>%
  dplyr::select(c("sample_size", "model", dplyr::everything())) %>%
  tidyr::gather(., key = "Tn", value = "coaltime", 3:ncol(.)) %>%
  dplyr::filter(!is.na(coaltime))



simdfplotObj <- simdf %>%
  ggplot() +
  geom_density(aes(x = coaltime, y = stat(count),
                   color = factor(model)), alpha = 0.5) +
  facet_grid(coaltime ~ Tn)


jpeg("msprimesims/compare_polysim_mpsrime.jpg", width = 11, height = 8, res = 300, units = "in")
plot(simdfplotObj)
graphics.off()



simdf.summary <- simdf %>%
  dplyr::group_by(model, Tn) %>%
  dplyr::summarise(
    minTime = min(coaltime),
    meanTime = mean(coaltime),
    medTime = median(coaltime),
    maxTime = max(coaltime),
    sdTime = sd(coaltime)
  ) %>%
  dplyr::arrange(., Tn)



