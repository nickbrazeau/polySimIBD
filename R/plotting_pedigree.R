#' @title Plot Between Lineage IBDness from Simulated Pedigree
#' @import ggplot2
#'
#' @importFrom magrittr %>%
#' @export

plot_between_lineages_simulatedpedigree <- function(simulatedpedigree){

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
    scale_fill_viridis_d("Between-Lineages IBD \n at Generation K") +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "top")




  #...........................
  # Looking at Full K
  #...........................
  f1line <- tibble::tibble(k = 0:simulatedpedigree$k, lineage = "F1.1_K")
  f1line$haplobit <- purrr::map(simulatedpedigree$f1.1lineage, "haplobit")
  f1line <- unnest(f1line) %>%
    dplyr::group_by(k) %>%
    dplyr::mutate( loci = 1:length(simulatedpedigree$f1.1lineage[[1]]@haplogt)) # just find number of sites


  f2line <- tibble::tibble(k = 0:simulatedpedigree$k, lineage = "F1.2_K")
  f2line$haplobit <- purrr::map(simulatedpedigree$f1.2lineage, "haplobit")
  f2line <- unnest(f2line) %>%
    dplyr::group_by(k) %>%
    dplyr::mutate( loci = 1:length(simulatedpedigree$f1.2lineage[[1]]@haplogt)) # just find number of sites


  plotObj2 <- rbind.data.frame(f1line, f2line) %>%
    ggplot() +
    geom_tile(aes(x=loci, y = factor(k), fill = factor(haplobit))) +
    facet_grid(~lineage) +
    scale_fill_viridis_d("Haplotype-Bit") +
    scale_y_discrete("K", breaks = 0:max(f1line$k)) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(face = "bold"))

  cowplot::plot_grid(plotObj1, plotObj2, nrow=2)


}




#' @title Plot Within Lineage IBD from Simulated Pedigree
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export


plot_within_lineage_simulatedpedigree <- function(sim, lineage){

  # assert lineage is 1 or 2
  if(lineage == 1){
    ret <- tibble::tibble(k = 0:(length(sim$f1.1lineage)-1))
    ret$lineage <- purrr::map(sim$f1.1lineage, "haplobit")
    ret <- unnest(ret) %>%
      dplyr::group_by(k) %>%
      dplyr::mutate(loci = 1:length(sim$f1.1lineage[[1]]@haplogt))
  } else if(lineage == 2){
    ret <- tibble::tibble(k = 0:(length(sim$f1.2lineage)-1))
    ret$lineage <- purrr::map(sim$f1.2lineage, "haplobit")
    ret <- unnest(ret) %>%
      dplyr::group_by(k) %>%
      dplyr::mutate(loci = 1:length(sim$f1.2lineage[[1]]@haplogt))
  }


  plotObj2 <- ret %>%
    ggplot() +
    geom_tile(aes(x=loci, y=factor(k), fill=factor(lineage))) +
    scale_fill_viridis_d("Haplotype-Bit") +
    scale_y_discrete("K", breaks = 0:max(ret$k)) +
    theme_classic()

  #...................
  # Find IBD Segments
  #....................
  IBDness <- ret %>%
    tidyr::spread(., key = "k", value = "lineage")

  IBDness.mat <- matrix(NA, nrow(IBDness), ncol(IBDness)-2)
  for(i in 3:ncol(IBDness)){
    for(j in 1:nrow(IBDness)){
      IBDness.mat[j, i-2] <- IBDness[j, 2] == IBDness[j, i]

    }
  }

  plotObj1 <- cbind.data.frame(loci = IBDness$loci, IBDness.mat) %>%
    tidyr::gather(., key = "k", value = "IBD", 2:ncol(.)) %>%
    dplyr::mutate(k = factor(k, levels = 1:ncol(IBDness.mat), labels = paste0("F", 1:ncol(IBDness.mat)))) %>%
    ggplot() +
    geom_tile(aes(x=loci, y=k, fill=factor(IBD))) +
    scale_fill_viridis_d("Progeny IBD \ to Parent") +
    theme_classic() +
    theme(axis.title.y = element_blank())

  cowplot::plot_grid(plotObj1, plotObj2, nrow=2)

}









