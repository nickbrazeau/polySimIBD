% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_realized_pairwise_ibd}
\alias{get_realized_pairwise_ibd}
\title{Effective IBD by Loci from SWF Simulation for a Pair of Hosts}
\usage{
get_realized_pairwise_ibd(swf, host_index = NULL)
}
\arguments{
\item{swf}{result of forwards simulation using the function \code{sim_swf()}}

\item{host_index}{a vector of target hosts. Defaults to all hosts}
}
\value{
double of pairwise IBD
}
\description{
Assumes that the minimum realized COI between the pairs of host determines
the denominator for the between realized IBD
}
\details{
Only accepts a pair of hosts. Ignores mutations as interrupting IBD segments.
}
