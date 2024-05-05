# Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#' The Structured Wright Fisher Model Output
#' This S3 class represents realization of the Discrete-Time Discrete-Loci Spatial Wright Fisher Malaria Model.
#' @name swfsim
#' @field pos vector; the genomic coordinates for chromosome and position of the sites
#' @field coi vector; the COI of each host 
#' @field recomb list; recombination blocks for each ancestral host for each generation 
#' @field parent_host1 list; parent for host 1 for each generation 
#' @field parent_host2 list; parent for host 2 for each generation 
#' @field parent_haplo1 list; haplotypes for parent 1 for each generation 
#' @field parent_haplo2 list; haplotypes for parent 2 for each generation 
#' @description The realization of the Discrete-Time Discrete-Loci Spatial Wright Fisher Malaria Model contains all of the information 
#' to create the ARG: each generation's parents, the resulting recombination events between parents, and the offspring haplotypes 
NULL


#' The Ancestral Recombination Graph 
#' This S3 class represents the ARG from the realized simulation
#' @name argraph
#' @description A list of \code{\link{bv_tree}} for each discrete-loci, which constitutes the ARG 
NULL



#' A Simplified Tree: bvtree 
#' This S3 class represents the bvtree which is a simple tree marginal representation  
#' @name bvtree
#' @field c vector; the node connection for each haplotype (each haplotype is an element in a vector)
#' @field t vector; the timing of the node connection (time to MRCA)
#' @field z vector; the order of coalescence for each set of haplotypes
#' @description  The \code{\link{bv_tree}} class is a lightweight representation of a marginal tree
NULL

