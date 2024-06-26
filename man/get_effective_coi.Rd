% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_results.R
\name{get_effective_coi}
\alias{get_effective_coi}
\title{Extract Effective COI by Loci from SWF Simulation for a Single Host}
\usage{
get_effective_coi(swf, host_index = NULL)
}
\arguments{
\item{swf}{result of forwards simulation using the function \code{sim_swf()}}

\item{host_index}{a vector of target hosts. Defaults to all hosts}
}
\value{
vector of effective COI by loci
}
\description{
From a single host in a SWF Simulation, extract the effective COI
for loci within the ARG. Effective COI is defined as the number of non-coalesced genomes
at the end of \code{tlim}. Note, this framework is an independent process for each recombination event
and thus will vary along the simulated genome (but not necessarily by locus).
}
\details{
Function limited to a single host per "realization"
}
