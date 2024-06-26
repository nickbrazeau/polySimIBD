% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forward_sWF.R
\name{sim_swf}
\alias{sim_swf}
\title{The Structured Wright Fisher Model for IBD}
\usage{
sim_swf(pos, N, m, rho, mean_coi, tlim, migr_mat = 1, verbose = FALSE)
}
\arguments{
\item{pos}{vector; the genomic coordinates for chromosome and position of the sites}

\item{N}{integer vector; The number of individuals to consider in each deme}

\item{m}{numeric numeric; Probability of internal migration where m represents
the probability of moving from host_{origin} to host_{new} by m*(1-1/N) of each deme}

\item{rho}{numeric; expected recombination rate}

\item{mean_coi}{numeric vector; The lambda of a right-shifted Poisson process,
1 + Pos(lambda) representing the average COI of each deme}

\item{tlim}{numeric; the maximum number of generations to consider before
exiting gracefully if all samples have not coalesced}

\item{migr_mat}{numeric matrix; Migrations rates or probabilities between
destination and origin. Note, because this is a Wright-Fisher model, we are
drawing parents and therefore migration matrix is parameterized towards
"where one came from" versus "where one is headed": origin specified as
columns and destination in rows. Default value of 1 indicates non-spatial
model. Note, if using a probability matrix, rows must sum to 1
(valid marginal probability); otherwise, values will be assumed to be rates
and converted to probabilities.}

\item{verbose}{boolean}
}
\value{
Returns a list of length six that contains \enumerate{
\item pos: The simulated genetic coordinates
\item coi: The COI of each individual
\item recomb: A recombination list of length of tlim where each element contains
the recombination block -- as a boolean -- of the two parental haplotypes.
\item parent_host1: the parental host assignments for the "paternal" haplotype
\item parent_host1: the parental host assignments for the "maternal" haplotype
\item parent_haplo1 "paternal" haplotype assigment (as above)
\item parent_haplo2 "maternal" haplotype assigment (as above)
}
}
\description{
Simulate a population forwards with recombination that approximates
the Structured Wright Fisher Process and tracks haplotype identity by descent
where individuals represent demes, such that within a deme individual-level COI is
considered. The model is also extended to consider spatial demes that individual hosts
can move between.
}
\details{
Demes are assumed to be ordered throughout (i.e. the order needs to be
consistent between N, m, mean_coi, and the rows and columns of the migration matrix).

The migration matrix is assumed to be a distance matrix that is either
a rate or a probability. The program will coerce the matrix into a probability
distribution between origin and destination based on the row-sums.

This function is intended to be fed into the \link{get_arg}
function to summarize the simulation results,
}
