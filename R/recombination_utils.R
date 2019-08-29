

makecrossover <- function(p1, p2, chrompos, rho){

  breakpoint <- rexp(1, rho)

  if(breakpoint > max(chrompos$POS)){ # recombination doesn't happen
    child <- sample(x=c(p1,p2), 1)[[1]] # remove list

  } else {
    child <- new("simhaplo")
    recombo.block <- chrompos$POS <= breakpoint
    # sample starting parent
    start <- sample(x=c(p1,p2), 1)[[1]]
    # determine end parent
    if(identical(p1, start)){
      end <- p2
    } else{
      end <- p1
    }
    # make child
    child@haplogt[recombo.block] <- start@haplogt[recombo.block]
    child@haplobit[recombo.block] <- start@haplobit[recombo.block]
    child@haplogt[!recombo.block] <- end@haplogt[!recombo.block]
    child@haplobit[!recombo.block] <- end@haplobit[!recombo.block]

    # redraw breakpoint to see if we do back crossover
    newbreakpoint <- breakpoint + rexp(1, rho)
    if(newbreakpoint < max(chrompos$POS)){ # we cross back over
      offset <- runif(1, breakpoint, newbreakpoint) # need to offset to allow cross back over
      # start now "re-starts"
      backcross.block <- newbreakpoint <= chrompos$POS
      # back cross child
      child@haplogt[backcross.block] <- start@haplogt[backcross.block]
      child@haplobit[backcross.block] <- start@haplobit[backcross.block]
    }

  } # end recombination ifelse

    return(child)
  }

