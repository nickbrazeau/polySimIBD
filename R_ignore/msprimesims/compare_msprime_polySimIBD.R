library(tidyverse)


#....................................
# Read in msprimes and manipulate
#....................................
msprimesims <- readRDS("simulations/msprimesims/simdata/msprimesims.rds")
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
msprimesims.df$pos <- 1:1001


#....................................
# Read in polysimibd
#....................................
polySimIBDsims.df <- readRDS("simulations/msprimesims/simdata/polySimIBD_sims.RDS") %>%
  dplyr::mutate(model = "polySimIBD") %>%
  dplyr::select(-c("intvl"))
polySimIBDsims.df$pos <- 1:1001

#################################################################################
###################                Comparisons                ###################
##################################################################################

#....................................
# CONTINUOUS time for msprime
#....................................
simdf <- rbind.data.frame(msprimesims.df, polySimIBDsims.df)

simdf <- simdf %>%
  dplyr::select(c("sample_size", "pos", "model", dplyr::everything())) %>%
  tidyr::gather(., key = "Tn", value = "coaltime", 4:ncol(.)) %>%
  dplyr::filter(!is.na(coaltime))



# drop rows for reasonable plot
rows <- sample(1:nrow(simdf), size = 1e4)
rows <- c( 1:nrow(simdf) %in% sort(rows) )
simdfplotObj <- simdf %>%
  dplyr::filter(rows) %>%
  ggplot() +
  geom_density(aes(x = coaltime, y = stat(count),
                   color = factor(model), fill = factor(model)), alpha = 0.1) +
  facet_wrap(Tn~.) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))


jpeg("simulations/msprimesims/compare_polysim_mpsrime_cont.jpg", width = 11, height = 8, res = 300, units = "in")
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



#....................................
# DISCRETE time for msprime
#....................................

# find discrete rows
discretepos <- seq(1, 1001, by = 100)
discreterows <- which(simdf$pos %in% discretepos)
discreterows <- c(1:nrow(simdf) %in% discreterows)



simdfplotObj <- simdf %>%
  dplyr::filter(discreterows) %>%
  ggplot() +
  geom_density(aes(x = coaltime, y = stat(count),
                   color = factor(model), fill = factor(model)), alpha = 0.1) +
  facet_wrap(Tn~.) +
  theme_classic()


jpeg("simulations/msprimesims/compare_polysim_mpsrime_discrete.jpg", width = 11, height = 8, res = 300, units = "in")
plot(simdfplotObj)
graphics.off()





simdf.summary <- simdf %>%
  dplyr::filter(discreterows) %>%
  dplyr::group_by(model, Tn) %>%
  dplyr::summarise(
    minTime = min(coaltime),
    meanTime = mean(coaltime),
    medTime = median(coaltime),
    maxTime = max(coaltime),
    sdTime = sd(coaltime)
  ) %>%
  dplyr::arrange(., Tn)






