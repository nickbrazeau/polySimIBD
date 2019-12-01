#' plot the bvtrees
#' @param ARGsim S4 object; The ARGsim object that was created from a structured Wright
#'        Fisher Simulation via `sim_structured_WF` and then processed with `get_ARG`.
#' @param loci numeric vector; Loci to plot
#' @import ggplot2
#' @return ggplot object of geom_segments
#' @export

plot_coalescence_trees <- function(ARGsim, loci){

  assert_custom_class(x = ARGsim, c = "ARGsim")

  # extract needed objects
  parasites <- ARGsim$parasites
  ARGsim <- ARGsim$ARG[loci]

  # make this into a tidy dataframe for ggplot
  ARGsimdf <- tibble::tibble(loci = loci)
  ARGsimdf$con <- purrr::map(ARGsim, "c")
  ARGsimdf$time <-  purrr::map(ARGsim, "t")

  # add in node information
  assert_same_length(ARGsimdf$time[[1]], parasites)
  ARGsimdf$parasites <- list(parasites)

  ARGsimdf <- ARGsimdf %>%
    tidyr::unnest(cols = c("loci", "parasites", "con", "time")) %>%
    dplyr::group_by(loci) %>%
    dplyr::mutate(time = ifelse(time == -1, max(time) + 5, time),
                  con = ifelse(con == -1, max(parasites), con))

  # take back to mapdf
  ARGsimdf.plotdf <- ARGsimdf %>%
    dplyr::group_by(loci) %>%
    tidyr::nest()

  make_marginal_tree_plot <- function(x){
    plotObj <- ggplot(data = x) +
      geom_segment(aes(x = parasites-0.02, xend = con,
                       y = time, yend = time), # horizontal lines
                   size = 2) +
      geom_segment(aes(x = parasites, xend = parasites,
                       y = 0, yend = time, color = factor(parasites)),
                   size = 2) +  # veritical lines
      xlab("Parasites") + ylab("Generations") +
      theme_bw() +
      scale_x_continuous(breaks = 1:max(ARGsimdf$parasites)) +
      scale_color_viridis_d() +
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
