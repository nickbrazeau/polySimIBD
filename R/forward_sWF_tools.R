#' internal class
#' @noRd
setClass("bvtree",
         slots=list(c="numeric", t="numeric", z="numeric"))


#' Find Coalescence
#' @param swf S4 object;
#' @importFrom magrittr %>%
#' @export
get_ARG <- function(swf){

  # assertions
  assert_custom_class(x = swf, c = "sWFsim")

  # breakdown for easier input into cpp
  L <- length(swf$pos)
  coi <- sum( swf$coi )
  G <- length(swf$anc)

  # need to go from list of matrices to array
  anc <- swf$anc
  # note we can have flucuating population sizes but not loci
  maxn <- max( sapply(anc, nrow) )

  anc.fill <- list()
  for(i in 1:length(anc)){
    expmat <- matrix(NA, nrow = maxn, ncol = L)
    ancdim <- dim(anc[[i]])
    expmat[1:ancdim[1], 1:ancdim[2]] <- anc[[i]]
    expmat[is.na(expmat)] <- -9 # missing is -9 for cpp int
    anc.fill[[i]] <- expmat
  }

  # get pair comparisons
  pairs <- as.data.frame( expand.grid(1:coi, 1:coi) )
  pairs <- pairs[pairs[,1] != pairs[,2], ]
  p1 <- pairs[,1]
  p2 <- pairs[,2]


  # storage objects for cpp output
  ARG <- lapply(1:L, function(x) return(new("bvtree")))
  coaltime.squaremat <- lapply(1:L, function(x) return(matrix(NA, nrow = coi, ncol = coi)))

  # get arguments in list form
  args <- list(L = L, G = G, anc = anc.fill, p1 = p1, p2 = p2, coi = coi)


  # main cpp function
  output_raw <- get_ARG_cpp(args)

  # process output
  message("processing output")

  t <- output_raw$t
  c <- output_raw$c
  z <- output_raw$z
  coal_times <- output_raw$coal_times



  # TODO
  # get coal time square matrix

  coaltime.squaremat[[l]] <-  pairs %>%
    tidyr::spread(., key = "Var2", value = "coaltime") %>%
    dplyr::select(-c("Var1")) %>%
    as.dist(.)

  ARGlist <- list(ARG = ARG,
                  coal_times = coaltime.squaremat)

  class(ARGlist) <- "ARGsim"

  return(ret)
}







#' plot the bvtrees
#' @param ARGsim S4 object;
#' @param loci numeric vector; loci which we want to plot
#' @importFrom magrittr %>%
#' @return ggplot object of geom_segments
#' @export

plot_coalescence_trees <- function(ARGsim, loci){
  assert_custom_class(x = ARGsim, c = "ARGsim")

  ARGsim <- ARGsim$ARG[loci]

  # make this into a tidy dataframe for ggplot
  ARGsimdf <- tibble::tibble(loci = loci)
  ARGsimdf$con <- purrr::map(ARGsim, "c")
  ARGsimdf$time <-  purrr::map(ARGsim, "t")
  # add in node information
  ARGsimdf$nodes <- list( 1:length(ARGsim[[1]]@c) )

  ARGsimdf <- ARGsimdf %>%
    tidyr::unnest(cols = c("loci", "nodes", "con", "time")) %>%
    dplyr::group_by(loci) %>%
    dplyr::mutate(time = ifelse(time == -1, max(time) + 5, time),
                  con = ifelse(con == -1, nodes, con))

  # make vertical lines
  plotObj <- ggplot(data = ARGsimdf) +
    geom_segment(aes(x = nodes-0.02, xend = con,
                     y = time, yend = time), # horizontal lines
                 size = 1.1) +
    geom_segment(aes(x = nodes, xend = nodes,
                     y = 0, yend = time, color = factor(nodes)), size = 1.1) +  # veritical lines
    facet_wrap(facets = loci, scales = "free_y") +
    xlab("Nodes") + ylab("Generations") +
    theme_bw() +
    scale_x_continuous(breaks = 1:max(ARGsimdf$nodes)) +
    scale_color_viridis_d() +
    theme(axis.title = element_text(family = "Helvetica", face = "bold"),
          axis.text = element_text(family = "Helvetica", face = "bold"),
          legend.position = "none",
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "#bdbdbd", size = 0.8)
    )

  return(plotObj)


}

