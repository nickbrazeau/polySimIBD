% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_realized_coi}
\alias{get_realized_coi}
\title{Extract Effective COI by Loci from SWF Simulation for a Single Host}
\usage{
get_realized_coi(swf, host_index = NULL)
}
\arguments{
\item{swf}{result of forwards simulation using the function \code{sim_swf()}}

\item{host_index}{a vector of target hosts. Defaults to all hosts}
}
\value{
vector of effective COI by loci
}
\description{
Extract Effective COI by Loci from SWF Simulation for a Single Host
}
\details{
Only accepts a single host
}
