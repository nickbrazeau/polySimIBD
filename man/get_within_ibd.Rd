% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_results.R
\name{get_within_ibd}
\alias{get_within_ibd}
\title{Calculate Within-Host IBD}
\usage{
get_within_ibd(swf, host_index = NULL)
}
\arguments{
\item{swf}{result of forwards simulation using the function \code{sim_swf()}}

\item{host_index}{a vector of target hosts. Defaults to all hosts}
}
\value{
double of within-host IBD
}
\description{
The within-host IBD is calculated as the number of strains that have
coalesced within the \code{tlim} at each loci divided by the original (i.e. not effective)
COI. As an example, consider that there are three strains (i.e. parasites) within a host and
that the parasite genome has ten equidistant loci with a single recombination breakpoint at loci 5 (i.e.).
Within this framework, we consider at loci 1:5 if 2/3 strains have coalesced, the
within-host IBD for this section is 2/3. Next, for loci 6:10 if no strains have coalesced
the within-host IBD is 0. Combining these results with-weighting for respective length/portion of the
genome (weights here are equal and therefore negligible) the overall within-host IBD is:
\deqn{\frac{2 + 0}{Host_{COI} - 1}}, where one is subtracted from the Host-COI for self-comparison,
which gives (3-1) + (3-1) (for each loci). Note, because we consider self comparisons, the
denominator is always less than the true COI.
}
\details{
Function limited to a single host per "realization"
}
