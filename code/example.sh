#!/bin/bash

#SBATCH --time=23:59:59
#SBATCH --mem=4096

python3 ${WDIR}/sphere-assembly/example.py --number ${N} -S 250
