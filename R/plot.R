#------------------------------------------------
# Expand a series of colours by interpolation to produce any number of colours
# from a given series. The pattern of interpolation is designed so that (n+1)th
# value contains the nth value plus one more colour, rather than being a
# completely different series. For example, running more_colours(5) and
# more_colours(4), the first 4 colours will be shared between the two series.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @noRd

more_colours <- function (n = 5, raw_cols = RColorBrewer::brewer.pal(10, "Paired")) {
  
  assert_single_pos_int(n, zero_allowed = FALSE)
  assert_string(raw_cols)
  assert_vector(raw_cols)
  
  my_palette <- colorRampPalette(raw_cols)
  if (n <= 2) {
    return(my_palette(3)[1:n])
  }
  n_steps <- ceiling(log(n - 1)/log(2))
  n_breaks <- 2^(1:n_steps) + 1
  s <- unlist(mapply(function(x) seq(0, 1, l = x), n_breaks, 
                     SIMPLIFY = FALSE))
  s <- s[!duplicated(s)]
  w <- match(s, seq(0, 1, l = n_breaks[n_steps]))
  w <- w[1:n]
  all_cols <- my_palette(n_breaks[n_steps])
  ret <- all_cols[w]
  return(ret)
}

#------------------------------------------------
#' plot the bvtrees
#' 
#' @description Produces a ggplot of bvtrees at every locus.
#' 
#' @param ARG S4 object; The ARG object that was created from a structured Wright
#'        Fisher Simulation via `sim_structured_WF` and then processed with `get_ARG`.
#' @param loci numeric vector; Loci to plot
#' 
#' @import ggplot2
#' @importFrom tibble tibble
#' @return ggplot object of geom_segments
#' @export

plot_coalescence_trees <- function(ARG, loci = NULL){

  # avoid "no visible binding" note
  root <- haplotype <- con <- time <- NULL
  
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

#------------------------------------------------
#' plot barcodes
#' 
#' @description Produces a ggplot of barcode counts, in which width
#'   indicates the number of times each barcode was seen, and all completely
#'   unique barcodes (only seen once) are indicated in grey.
#' 
#' @param hapmat a matrix of haplotypes.
#' @param coi a vector of COIs per sample.
#' 
#' @importFrom grDevices grey
#' @export

plot_barcodes <- function(hapmat, coi) {
  
  # avoid "no visible binding" note
  width <- hap_ID <- NULL
  
  # assertions
  assert_matrix(hapmat)
  assert_vector(coi)
  assert_eq(sum(coi), ncol(hapmat))
  
  # give each observed haplotype a unique index and compute the vector of these
  # indices
  haplist <- split(t(hapmat), f = 1:ncol(hapmat))
  hapvec <- match(haplist, unique(haplist))
  
  # split by COI and remove duplicate haplotypes within each host
  hap_raw <- split(hapvec, f = rep(1:length(coi), times = coi))
  hap_final <- mapply(unique, hap_raw, SIMPLIFY = FALSE)
  
  # subset to monoclonal
  w <- which(mapply(length, hap_final) == 1)
  if (length(w) == 0) {
    stop("no monoclonal samples to plot")
  }
  hap_mono <- unlist(hap_final[w])
  
  # tabulate
  tab1 <- table(hap_mono)
  tab1 <- sort(tab1, decreasing = TRUE)
  
  # make plotting data.frame. Unique haplotypes are counted together, but are coloured grey
  df_hap <- NULL
  if (any(tab1 > 1)) {
    df_hap <- rbind(df_hap,
                    data.frame(hap_ID = names(tab1)[tab1 > 1],
                               width = as.vector(tab1)[tab1 > 1])
    )
    df_hap$col <- more_colours(nrow(df_hap))
  }
  if (any(tab1 == 1)) {
    df_hap <- rbind.data.frame(df_hap,
                               list(hap_ID = "unique",
                                    width = sum(tab1 == 1),
                                    col = grey(0.8)))
  }
  
  # plot
  df_hap %>%
    ggplot(aes(width = width)) + theme_void() +
    geom_bar(aes(x = cumsum(width) - width / 2, y = 1, fill = hap_ID),
             col = "black", stat = "identity") +
    scale_fill_manual(values = df_hap$col, guide = "none") +
    ylim(c(-3, 3))
}
