# News

## polySimIBD v1.1.0
Finalized documentation and additional functions for the spatial discrete-time discrete-loci structured Wright-Fisher Model. New functions include a C++ version of the IBD calculation (`get_bvibd`) and a function to convert marginal _bvtrees_ into Newick trees (`bvtreeToNewick`). 

## polySimIBD v1.0.0
Extension of the discrete-time discrete-loci structured Wright-Fisher Model for spatial "meta-populations". This model now assumes that there are spatial "meta-demes" and within each deme there are multiple hosts, each of who is considered as a "micro-deme" to model COI.


## polySimIBD v0.5.1
Add "get realized"" functions for COI and IBD. 

## polySimIBD v0.5.0
The original discrete-time discrete-loci structured Wright-Fisher Model. A single "meta-deme" is considered with each host considered a "micro-deme".
The full mathematical model is described in: Verity, Aydemir, Brazeau _et al._ 20202 Nat. Comms. ([PMC7192906](https://pubmed.ncbi.nlm.nih.gov/32355199/)).
