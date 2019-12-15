#--------------------------------------------------------------------
# Purpose of this script is to PROCESS the simulations that will allow us to 
# estimate the bias of the
# MLE inbreeding function from the Verity Nat Comm Submission
# using a discrete loci, discrete time structured wright fisher model
#--------------------------------------------------------------------
library(tidyverse)
library(polySimIBD)

# rsync -av nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/polySimIBD/R_ignore/NatComms_VerityAB2019_Sims/results/_rslurm_verity_nat_comm_sims/ /Users/nickbrazeau/Documents/GitHub/polySimIBD/R_ignore/NatComms_VerityAB2019_Sims/results/ 

plot_theme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
                    axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
                    axis.text.y = element_text(family = "Helvetica", hjust = 0.5, size = 11),
                    axis.text.x = element_text(family = "Helvetica", hjust = 0.5, size = 11, angle = 90),
                    legend.position = "right",
                    legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.85, size = 12),
                    legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
                    axis.line = element_line(color = "#000000", size = 1))


#..............................................................
# read in sims
#..............................................................
paramsdf <- readRDS("R_ignore/NatComms_VerityAB2019_Sims/results/params.RDS")

simfiles <- list.files("R_ignore/NatComms_VerityAB2019_Sims/results/",
                       pattern = ".RDS", full.names = T)
simfiles <- simfiles[!grepl("f.RDS|params.RDS", simfiles)]

simfiles.df <- tibble::tibble(
  path = simfiles ) %>% 
  dplyr::mutate(
  iter = stringr::str_split_fixed(simfiles, "/", 3)[,3],
  iter = stringr::str_extract(iter, "[0-9]+"),
  iter = as.numeric(iter)
  )  %>% 
  dplyr::arrange(iter)

paramsdf$results <- unlist( purrr::map(simfiles.df$path, function(x) readRDS(x)), 
                            recursive = F) # extract results -- use unlist without recursion to match


#..............................................................
# plot results
#..............................................................
plotdf <- paramsdf %>% 
  tidyr::unnest(cols = "results")

# housekeeping
coilvls <- levels(factor(plotdf$mean_coi))

plotObj <- plotdf %>% 
  tidyr::gather(., key = "IBD", value = "IBDest", 11:12) %>% 
  dplyr::group_by(mean_coi, m, N, IBD) %>% 
  dplyr::summarise(
    n = n(),
    meanIBD = mean(IBDest),
    seIBD = sd(IBDest)/sqrt(n),
    IBDLL = meanIBD - 1.96*seIBD,
    IBDUL = meanIBD + 1.96*seIBD
  ) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(mean_coi = factor(mean_coi, 
                                  levels = coilvls,
                                  labels = c("Mean COI: 1", "Mean COI: 2", "Mean COI: 3")),
                m = factor(m, 
                           levels = c("0", "0.25", "0.5", "1"),
                           labels = c("Migr: 0", "Migr: 0.25", "Migr: 0.5", "Migr: 1")),
                IBD = factor(IBD, 
                             levels = c("btwnIBD", "malecotf"),
                             labels = c("Truth", "MLE")),
                logN = log10(N)) %>% 
  ggplot() + 
  geom_pointrange(aes(x = logN, y = meanIBD, ymin = IBDLL, ymax = IBDUL, 
                      group = factor(IBD), color = factor(IBD)), alpha = 0.5) +
  scale_color_manual("IBD Measure", values = c("#54278f", "#e08214")) +
  facet_grid(mean_coi ~ m) + 
  ylab("IBD") + xlab("Effective Population (log10-transformed)") +
  plot_theme 

jpeg("R_ignore/NatComms_VerityAB2019_Sims/MLE_vs_DLWTsWFmod.jpg", 
     width = 11, height = 8, units = "in", res = 300)
plot(plotObj)
graphics.off()





