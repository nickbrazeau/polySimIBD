#--------------------------------------------------------------------
# Purpose of this script is to PROCESS the simulations that will allow us to 
# estimate the bias of the
# MLE inbreeding function from the Verity Nat Comm Submission
# using a discrete loci, discrete time structured wright fisher model
#--------------------------------------------------------------------
library(tidyverse)
library(polySimIBD)


#..............................................................
# Plot Results
#..............................................................
plot_theme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
                    axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
                    axis.text.y = element_text(family = "Helvetica", hjust = 0.5, size = 11),
                    axis.text.x = element_text(family = "Helvetica", hjust = 0.5, size = 11, angle = 90),
                    legend.position = "right",
                    legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.85, size = 12),
                    legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
                    axis.line = element_line(color = "#000000", size = 1))

plotdf <- paramsdf %>% 
  dplyr::mutate(simout = purrr::map(paramsdf$nat_com_sims, "sim_out")) %>% 
  tidyr::unnest(cols = simout)

plotObj <- plotdf %>% 
  tidyr::gather(., key = "IBD", value = "IBDest", 12:16) %>% 
  dplyr::filter(IBD %in% c("btwnIBDprop", "malecotf")) %>% 
  dplyr::group_by(mean_coi, m, N, IBD) %>% 
  dplyr::summarise(
    n = n(),
    meanIBD = mean(IBDest),
    seIBD = sd(IBDest)/sqrt(n),
    IBDLL = meanIBD - 1.96*seIBD,
    IBDUL = meanIBD + 1.96*seIBD
    #meanIBD = mean(IBDest > 0.9)
    
  ) %>% 
  dplyr::mutate(logN = log10(N)) %>% 
  ggplot() + 
  geom_pointrange(aes(x = logN, y = meanIBD, ymin = IBDLL, ymax = IBDUL, 
                      group = factor(IBD), color = factor(IBD)), alpha = 0.5) +
  #geom_point(aes(x = logN, y = meanIBD,
  #               group = factor(IBD), color = factor(IBD)), alpha = 0.5) +                     
  #scale_color_manual("IBD Measure", values = c("#542788", "#f768a1", "#e08214", "#2171b5", "#08519c")) +
  scale_color_manual("IBD Measure", values = c("#006d2c", "#542788", "#e08214")) +
  #scale_color_manual("IBD Measure", values = c("#006d2c", "#e08214", "#ef6548")) +
  facet_grid(mean_coi ~ m) + 
  ylab("IBD") + xlab("Effective Population (log10-transformed)") +
  plot_theme 

jpeg("~/Desktop/temp_polysimibd_btwnIBD_nomut_rho.jpg", width = 11, height = 8, units = "in", res = 300)
plot(plotObj)
graphics.off()





