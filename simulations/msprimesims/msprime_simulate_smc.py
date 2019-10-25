#!/usr/local/bin/python3

import sys
import os
import argparse
import msprime
import tskit


# simulate tree under smc model
def wrap_msprime(sample_size, Ne, length, recombination_rate, model, outpath, repnum):
    tree_sequence = msprime.simulate(
      sample_size=sample_size, Ne=Ne, # note, it thinks diploid pop
      length=length, recombination_rate=recombination_rate,
      model=model)
    # recursively make out directories
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    # write out simulations to tab delimited file
    f = open(os.path.join(outpath, "msprimesim_"+repnum+".tab.txt"), "w")
    for tree in tree_sequence.trees():
        for u in tree.nodes():
            print(tree.index, tree.interval, u, tree.parent(u), tree.time(u), sep="\t", file = f )
    f.close()



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_size", required = True, type = int, help = "The number of sampled individual in the population. For more information, see msprime documentations")
    parser.add_argument("--Ne", required = True, type = float, help = "Effective population size. For more information, see msprime documentations")
    parser.add_argument("--length", required = True, type = float, help = "Length of chromsome. For more information, see msprime documentations")
    parser.add_argument("--recombination_rate", required = True, type = float, help = "Rate per base per generation. For more information, see msprime documentations")
    parser.add_argument("--model", required = True, type = str, help = " For more information, see msprime documentations")
    parser.add_argument("--outpath", required = True, type = str, help = "Path for tab-delimited file for msprime output")
    parser.add_argument("--replicates", required = True, type = int, help = "Number of model replicates to perform")
    return parser.parse_args()


def main():
    args = parse_args()
    for i in range(0, args.replicates):
        wrap_msprime(args.sample_size, args.Ne, args.length, args.recombination_rate, args.model, args.outpath, str(i))
    return 0


if __name__== '__main__':
    main()
