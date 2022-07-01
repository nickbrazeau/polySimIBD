#' plot the bvtrees
#' @param ARG S4 object; The ARG object that was created from a structured Wright
#'        Fisher Simulation via `sim_structured_WF` and then processed with `get_ARG`.
#' @param loci numeric vector; Loci to plot
#' @import ggplot2
#' @return ggplot object of geom_segments
#' @export

plot_coalescence_trees <- function(ARG, loci = NULL){

  # assertions
  if(!is.null(loci)){
    assert_vector(loci)
    assert_numeric(loci)
    # subset if not null
    ARG <- ARG[[loci]]
  }
  
  if (length(ARG) == 1) {
    assert_custom_class(ARG, "bvtree", message  =  "Elements within the %s must inherit from class '%s'")
  } else {
    assert_custom_class(ARG[[1]], "bvtree", message  =  "Elements within the %s must inherit from class '%s'")
  }


  # coerce this into a tidy format
  if (length(ARG) == 1) {
    ARG.tidy <-  tibble::tibble(
        loci = loci, 
        c = ARG@c,
        t = ARG@t,
        z = ARG@z
      ) 
  } else {
    loci <- 1:length(ARG)
    ARG.tidy <- lapply(loci, function(x){
      ret <- tibble::tibble(
        loci = x, 
        c = ARG[[x]]@c,
        t = ARG[[x]]@t,
        z = ARG[[x]]@z,
      )
      return(ret)
    }) %>% 
      dplyr::bind_rows()
  }

  
  # tidy up inputs
  ARGsimdf <- ARG.tidy %>%
    dplyr::group_by(loci) %>%
    dplyr::mutate(time = ifelse(t == -1, max(t) + 5, t),
                  root = which(t == -1)[1], # pick first root
                  haplotype = seq(from = 0, to = (length(t)-1), by = 1),
                  con = ifelse(t == -1, root, c))

  # take back to mapdf because facet has issue with multiple lines
  ARGsimdf.plotdf <- ARGsimdf %>%
    dplyr::group_by(loci) %>%
    tidyr::nest()

  make_marginal_tree_plot <- function(x){
    plotObj <- ggplot(data = x) +
      geom_segment(aes(x = haplotype, xend = con,
                       y = time, yend = time), # horizontal lines
                   size = 1) +
      geom_segment(aes(x = haplotype, xend = haplotype,
                       y = 0, yend = time, color = factor(haplotype)),
                   size = 1) +  # veritical lines
      xlab("Haplotypes") + ylab("Generations") +
      theme_bw() +
      scale_x_continuous(breaks = 1:max(ARGsimdf$haplotype)) +
      theme(axis.title = element_text(family = "Helvetica", face = "bold", size = 16),
            axis.text.x = element_blank(),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 16),
            legend.position = "none",
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "#bdbdbd", size = 1.1)
      )
    return(plotObj)
  }

  ARGsimdf.plotdf$plots <- purrr::map(ARGsimdf.plotdf$data, make_marginal_tree_plot)

  ret <- cowplot::plot_grid(plotlist = ARGsimdf.plotdf$plots, align = "v")
  return(ret)

}
