#' internal class
#' @noRd
setClass("bvtree",
         slots=list(c="numeric", t="numeric", z="numeric"))


#' Find Coalescence
#' @param swf S4 object;
#' @importFrom magrittr %>%
#' @export
get_ARG <- function(swf){

  assert_custom_class(x = swf, c = "sWFsim")
  L <- length(swf$pos)
  anc <- swf$anc
  coi <- swf$coi

  ARG <- lapply(1:L, function(x) return(new("bvtree")))
  coaltime.squaremat <- lapply(1:L, function(x) return(matrix(NA, nrow = sum(coi), ncol = sum(coi))))

  #...................
  # get t and c
  #...................
  for(l in 1:L){ # for each loci find time that pair coalesces
    pairs <- as.data.frame( expand.grid(1:sum(coi), 1:sum(coi)) )
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

    # extract lightweight infromation for bvtree
    nodes <- sort( unique(pairs[,1]) )

    # make connections for nodes always go right
    for(i in 1:(length(nodes)-1)){
      mintime <- min( pairs[ pairs$Var1 == nodes[i] & pairs$Var2 > i, ]$coaltime )

      # only do look ahead when we aren't one right of root
      if(i == length(nodes)-1){
        mintime.forward <- Inf
      } else {
        # loop through all pairs ahead
        mintime.ahead <- c()
        for(j in (i+1):(length(nodes)-1)){
          mintime.ahead <- c(mintime.ahead, min( pairs[ pairs$Var1 == nodes[j] & pairs$Var2 > j, ]$coaltime ))
        }
        mintime.ahead <- min(mintime.ahead)
      }

      lineage <- pairs$Var2[ pairs$Var1 == nodes[i] & pairs$coaltime == mintime & pairs$Var2 > i ] # always make trees look "right" via pairs$Var2 > i

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
    ARG[[l]]@c[nodes[length(nodes)]] <- -1
    ARG[[l]]@t[nodes[length(nodes)]] <- -1


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
      as.dist(.)

  } #end for loop for loci


  ARGlist <- list(ARG = ARG,
                  coal_times = coaltime.squaremat)

  class(ARGlist) <- "ARGsim"

  return(ARGlist)
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
