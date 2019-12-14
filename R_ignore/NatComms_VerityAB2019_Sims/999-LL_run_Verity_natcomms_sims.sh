#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=49512
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

~/.linuxbrew/bin/R -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/polySimIBD"); source("R_ignore/NatComms_VerityAB2019_Sims/run_Verity_natcomms_sims.R")'
