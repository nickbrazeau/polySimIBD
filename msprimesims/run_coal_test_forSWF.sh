#!/usr/bin/env bash

EXC=/Users/nickbrazeau/Documents/GitHub/polyARG/msprime_sanity_checks/
OUTDIR=/Users/nickbrazeau/Desktop/msprimesims/

# discrete Wright-Fisher
python3 $EXC/msprime_simulate_smc.py --sample_size 2 --Ne 5e2 --length 1e3 --recombination_rate 1e-4 --model dtwf --outpath $OUTDIR/discreteWF/smplsize2/ --replicates 1000
python3 $EXC/msprime_simulate_smc.py --sample_size 3 --Ne 5e2 --length 1e3 --recombination_rate 1e-4 --model dtwf --outpath $OUTDIR/discreteWF/smplsize3/ --replicates 1000
python3 $EXC/msprime_simulate_smc.py --sample_size 5 --Ne 5e2 --length 1e3 --recombination_rate 1e-4 --model dtwf --outpath $OUTDIR/discreteWF/smplsize5/ --replicates 1000
