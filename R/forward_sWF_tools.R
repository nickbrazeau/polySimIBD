#' internal class
#' @noRd
setClass("bvtree",
         slots=list(c="numeric", t="numeric", z="numeric"))


#' Find Coalescence
#' @param swf S4 object; A discrete-loci, discrete-time structured Wright Fisher Simulation
#'            called from `sim_structured_WF`.
#' @param parasites int vector; Terminal nodes, or parasites from the Structured Wright Fisher
#'        Simulation for consideration.
#' @export

get_ARG <- function(swf, parasites = NULL){

  # graceful exits
  if(sum(swf$coi) == 1){
    warning("Your simulation returned one parasite among all hosts (did you run N = 1 and a low mean_moi?). As such, no pairwise comparisons can be made.")
    return("Only one parasite, no inference can be made")
  }

  # assertions
  assert_custom_class(x = swf, c = "sWFsim")
  assert_gr(length(parasites), 1, message = "Parasites must be a vector of terminal nodes for consideration that is longer than 1.")

  # tidy up params
  L <- length(swf$pos)
  anc <- swf$anc
  coi <- sum(swf$coi)

  if(is.null(parasites)){
    parasites <- 1:coi
  }


  ARG <- lapply(1:L, function(x) return(new("bvtree")))
  coaltime.squaremat <- lapply(1:L, function(x) return(matrix(NA, nrow = length(parasites), ncol = length(parasites))))

  #...................
  # get t and c
  #...................
  for(l in 1:L){ # for each loci find time that pair coalesces
    pairs <- as.data.frame( expand.grid(parasites, parasites) )
    pairs <- pairs[pairs[,1] != pairs[,2], ]
    pairs$coaltime <- rep(NA, nrow(pairs))
    for(i in 1:nrow(pairs)){ # find when each pair coalesces

      p1 <- pairs[i,1] # first "pointer"
      p2 <- pairs[i,2] # second "pointer"

      for (g in length(anc):1) {
        p1 <- anc[[g]][cbind(p1, l)]
        p2 <- anc[[g]][cbind(p2, l)]
        if(is.na(pairs$coaltime[i])){
          if(p1 == p2) {
            pairs$coaltime[i] <- length(anc)-g + 1 # 1-based
          }
        }
      } # end for loop for g
    } # end for loop for pairs

    # if pairs do not have a common ancestor, then set to tlim
    if(any(is.na(pairs$coaltime))){
      pairs$coaltime[is.na(pairs$coaltime)] <- -1
      warning("Your simulation did not fully coalesce. Pairs that did not reach a common ancestor have been set to the t-limit of -1.")
    }


    # extract lightweight infromation for bvtree
    # make connections for parasites always go right
    for(i in 1:(length(parasites)-1)){
      mintime <- min( pairs[ pairs$Var1 == parasites[i] & pairs$Var2 > i, ]$coaltime )

      # only do look ahead when we aren't one right of root
      if(i == length(parasites)-1){
        mintime.forward <- Inf
      } else {
        # loop through all pairs ahead
        mintime.ahead <- c()
        for(j in (i+1):(length(parasites)-1)){
          mintime.ahead <- c(mintime.ahead, min( pairs[ pairs$Var1 == parasites[j] & pairs$Var2 > j, ]$coaltime ))
        }
        mintime.ahead <- min(mintime.ahead)
      }

      lineage <- pairs$Var2[ pairs$Var1 == parasites[i] & pairs$coaltime == mintime & pairs$Var2 > i ] # always make trees look "right" via pairs$Var2 > i

      if(length(lineage) > 1){ # multiple coal
        if(mintime >= mintime.ahead){ # look ahead to make sure we don't block future branches
          lineage <- max(lineage) # pick farthest right
        } else {
          lineage <- lineage[1] # pick first (most left)
        }
      }
      ARG[[l]]@c[i] <- lineage
      ARG[[l]]@t[i] <- mintime
    }

    # note last node must always be root (if we are always going right)
    ARG[[l]]@c[length(parasites)] <- -1
    ARG[[l]]@t[length(parasites)] <- -1


    #...................
    # get z
    #...................
    ord <- sort(ARG[[l]]@t)
    ord.un <- unique(ord)
    ord.un <- ord.un[ord.un != -1] # remove root
    ord.num <- 1:length(ord.un)

    ARG[[l]]@z <-  ARG[[l]]@t # temporary overwrite

    for(i in 1:length(ord.un)){
      ARG[[l]]@z[ ARG[[l]]@z == ord.un[i] ] <- ord.num[i]
    }
    # end section for Z fill in

    #...................
    # get coal time square matrix
    #...................
    coaltime.squaremat[[l]] <-  pairs %>%
      tidyr::spread(., key = "Var2", value = "coaltime") %>%
      dplyr::select(-c("Var1")) %>%
      stats::as.dist()

  } #end for loop for loci


  ARGlist <- list(ARG = ARG,
                  coal_times = coaltime.squaremat,
                  coi = swf$coi,
                  parasites = parasites)

  class(ARGlist) <- "ARGsim"

  return(ARGlist)
}





#' plot the bvtrees
#' @param ARGsim S4 object; The ARGsim object that was created from a structured Wright
#'        Fisher Simulation via `sim_structured_WF`` and then processed with `get_ARG`.
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
