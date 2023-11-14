#------------------------------------------------
# Expand a series of colors by interpolation to produce any number of colors
# from a given series. The pattern of interpolation is designed so that (n+1)th
# value contains the nth value plus one more color, rather than being a
# completely different series. For example, running more_colors(5) and
# more_colors(4), the first 4 colors will be shared between the two series.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @noRd

more_colors <- function (n = 5, raw_cols = RColorBrewer::brewer.pal(10, "Paired")) {
  
  goodegg::assert_single_pos_int(n, zero_allowed = FALSE)
  goodegg::assert_string(raw_cols)
  goodegg::assert_vector(raw_cols)
  
  my_palette <- grDevices::colorRampPalette(raw_cols)
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
#' plot barcodes
#' 
#' @description Produces a ggplot of barcode counts, in which width
#'   indicates the number of times each barcode was seen, and all completely
#'   unique barcodes (only seen once) are indicated in grey
#' 
#' @param hapmat matrix; a matrix of haplotypes
#' @param coi numeric vector; a vector of COIs per sample
#' 
#' @importFrom grDevices grey
#' @export

plot_barcodes <- function(hapmat, coi) {
  
  # avoid "no visible binding" note
  width <- hap_ID <- NULL
  
  # assertions
  goodegg::assert_matrix(hapmat)
  goodegg::assert_vector(coi)
  goodegg::assert_eq(sum(coi), ncol(hapmat))
  
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
    ggplot2::ggplot(aes(width = width)) + ggplot2::theme_void() +
    ggplot2::geom_bar(aes(x = cumsum(width) - width / 2, y = 1, fill = hap_ID),
             col = "black", stat = "identity") +
    ggplot2::scale_fill_manual(values = df_hap$col, guide = "none") +
    ggplot2::ylim(c(-3, 3))
}
