#' @title Plot Simulated Pedigree
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export

plot.simulatedpedigree <- function(simulatedpedigree){

  # assert that class simulated pedigree

  #...........................
  # Looking at result prog
  #...........................
  progf1k <- data.frame(lineage = "F1.1_K",
                        loci = 1:length(unlist(simulatedpedigree$kprogeny1@haplobit)),
                        haplobit1 = unlist(simulatedpedigree$kprogeny1@haplobit)
  )
  progf2k <- data.frame(lineage = "F1.2_K",
                        loci = 1:length(unlist(simulatedpedigree$kprogeny2@haplobit)),
                        haplobit2 = unlist(simulatedpedigree$kprogeny2@haplobit)
  )

  progfk <- dplyr::left_join(progf1k, progf2k, by = "loci")
  progfk$IBD <- mapply(identical, as.character(progfk$haplobit1),
                            as.character(progfk$haplobit2))

  plotObj1 <- progfk %>%
    ggplot() +
    geom_tile(aes(x=loci, y=lineage.x, fill=factor(IBD))) +
    scale_fill_viridis_d("IBD") +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank())




  #...........................
  # Looking at Full K
  #...........................
  f1line <- data.frame(k = 1:simulatedpedigree$k, lineage = "F1.1")
  f1line$haplobit <- purrr::map(simulatedpedigree$f1.1lineaage, "haplobit")
  f1lines2c <- lapply(f1line$haplobit, function(x){unlist(x)}) %>%
    do.call("rbind", .)
  f1line <- cbind.data.frame(f1line, f1lines2c)

  f2line <- data.frame(k = 1:simulatedpedigree$k, lineage = "F1.2")
  f2line$haplobit <- purrr::map(simulatedpedigree$f1.2lineaage, "haplobit")
  f2lines2c <- lapply(f2line$haplobit, function(x){unlist(x)}) %>%
    do.call("rbind", .)
  f2line <- cbind.data.frame(f2line, f2lines2c)

  plotObj2 <- rbind.data.frame(f1line, f2line) %>%
    dplyr::select(-c("haplobit")) %>%
    tidyr::gather(key = "loci", value = "haplo", 3:ncol(.)) %>%
    dplyr::mutate(haplo = factor(haplo),
                  loci = as.numeric(loci)) %>%
    ggplot() +
    geom_tile(aes(x=loci, y = k, fill = haplo)) +
    facet_grid(~lineage) +
    scale_fill_viridis_d("Haplotype-Bit") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(face = "bold"))

  cowplot::plot_grid(plotObj1, plotObj2, align = "h",
                    nrow=2)


}
