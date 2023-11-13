#' @title Convert `bvtree` to Newick Tree Format
#' @param bvtree S3 class; internal class for the `polySimIBD`: package \link{bvtree}
#' @inheritParams sim_swf
#' @description Recursively converts `bvtree` to Newick tree format for compatibility with other downstream packages
#' @details Note, the overall TMRCA is set to the \code{tlim}, which in Newick format looks inappropriately rooted 
#' @returns Newick String 
#' @export

bvtreeToNewick <- function(bvtree, tlim = 10) {
  #............................................................
  # checks
  #............................................................
  goodegg::assert_eq(x = "bvtree", y = class(bvtree)[[1]],
                     message = "bvtree must be of class bvtree from polySimIBD")    
  goodegg::assert_eq(x = "polySimIBD", y = attr(class(bvtree),"package"),
                     message = "bvtree must be made with polySimIBD")    
  goodegg::assert_single_int(tlim)
  
  #............................................................
  # core
  #............................................................
  buildNewick <- function(node, parentTime = tlim) {  # tlim for root time 
    children <- which(bvtree@c == (node-1)) # sim code 0-based 
    
    # Calculate branch length
    if (bvtree@t[node] == -1) {
      branchLength <- parentTime
    } else {
      branchLength <- bvtree@t[node]
    }
    
    # If the node is a leaf
    if (length(children) == 0) {
      return(paste0("Node", (node-1), ":", branchLength)) # c/w sim code 
    }
    
    # Recursively build Newick string for children
    childrenStrings <- sapply(children, function(child){buildNewick(node = child, parentTime = tlim)})
    childrenNewick <- paste(childrenStrings, collapse = ",")
    
    #......................
    # Format the Newick string for this node and it's parent
    #......................
    # parent coalescence time w/ child under bvtree is same > need to find that time and convert to integer
    # going to use it again to make equal total branch lengths for node(s)  
    parentcoaltime <- as.integer( substring(childrenNewick, regexpr(":", childrenNewick) + 1) ) 
    return(paste0("(", "Node", (node-1), ":", # parent base
                  parentcoaltime, # same as child TMRCA
                  ",", childrenNewick, "):",  # child base + TMRCA
                  parentTime - parentcoaltime # equal total branch length 
    ))
  } # end recursion fxn

# Find the ancestral nodes
ancestralNodes <- which(bvtree@c == -1) # R is 1-based, Cpp sim code is 0-based


#............................................................
# out
#............................................................
# Build Newick strings for each ancestral node and concatenate
newickStrings <- sapply(ancestralNodes, buildNewick, parentTime = tlim)
fullNewick <- paste(newickStrings, collapse = ",")
fullNewick <- paste0("(", fullNewick, ");") 
return(fullNewick)
}

