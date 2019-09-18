

makecrossover <- function(p1, p2, chrompos, rho){

  # TODO do this allow too much "backcrossing" ?

  breakpoint <- rexp(1, rho)

  if(breakpoint > max(chrompos$POS)){ # recombination doesn't happen
    m1 <- p1
    m2  <- p2

  } else {
    m1 <- new("simhaplo")
    m2 <- new("simhaplo")
    recombo.block <- chrompos$POS <= breakpoint
    # sample starting parent
    start <- sample(x=c(p1,p2), 1)[[1]]
    # determine end parent
    if(identical(p1, start)){
      end <- p2
    } else{
      end <- p1
    }
    # make m1 child
    m1@haploint[recombo.block] <- start@haploint[recombo.block]
    m1@haploint[!recombo.block] <- end@haploint[!recombo.block]

    # make m2 child
    m2@haploint[recombo.block] <- end@haploint[recombo.block]
    m2@haploint[!recombo.block] <- start@haploint[!recombo.block]



    # redraw breakpoint to see if we do back crossover
    newbreakpoint <- breakpoint + rexp(1, rho)
    if(newbreakpoint < max(chrompos$POS)){ # we cross back over
      offset <- runif(1, breakpoint, newbreakpoint) # need to offset to allow cross back over
      # start now "re-starts"
      backcross.block <- newbreakpoint <= chrompos$POS
      # back cross m1
      m1@haploint[backcross.block] <- start@haploint[backcross.block]
      # back cross m2
      m2@haploint[backcross.block] <- end@haploint[backcross.block]
    }

  } # end recombination ifelse
  children <- list(m1 =m1, m2=m2)
  return(children)
}
