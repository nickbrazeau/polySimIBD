% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bvtree}
\alias{bvtree}
\title{A Simplified Tree: bvtree
This S3 class represents the bvtree which is a simple marginal tree representation}
\description{
The \code{bv_tree} class is a lightweight representation of a marginal tree
}
\section{Fields}{

\describe{
\item{\code{c}}{vector; the node connection for each haplotype (each haplotype is an element in a vector)}

\item{\code{t}}{vector; the timing of the node connection (time to MRCA)}

\item{\code{z}}{vector; the order of coalescence for each set of haplotypes}
}}

