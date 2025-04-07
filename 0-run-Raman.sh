#! /bin/bash

#===============================================================================#
# date: 2025.04.05 wdc
# purpose: calculate Raman spectrum by using finite displacement method.
# usage  : ./0-run-Raman.sh
#===============================================================================#
# description:
# step0: calculat PHONON by using VASP
# step1: generate the displacement POSCAR by reading OUTCAR.phon
# step2: VASP calculate for each displacement POSCAR, pick one of them
# step3: calculate the raman spectrum by each OUTCAR
# 1. VASP_RAMAN_PARAMS='[first-mode]_[last-mode]_[nderiv]_[step-size]
# detailed structrue-tree in README.md
#===============================================================================#

set -e

first="01"
last=$(grep "f  =" OUTCAR.phon | tail -1 | awk '{ print $1 }')
phonon_POTIM="0.01"

## description 1. ##
export VASP_RAMAN_PARAMS=${first}"_"${last}"_2_"${phonon_POTIM}

echo ${VASP_RAMAN_PARAMS}

## step1: generate the displacement POSCAR by reading OUTCAR.phon
# chmod +x 1_raman_gen_POS.py
#./1_raman_gen_POS.py &> gen.out

## step2: VASP calculate for each displacement POSCAR, pick one of them
## [ method-I ]
#sbatch -J raman 2_raman_sbatch
## [ method-II ]
# chmod +x 2_raman_vasp.sh
#./2_raman_vasp.sh

## step3: calculate the raman spectrum by each OUTCAR
# chmod +x 3_raman_spectrum.py 4_plot_raman.py
#./3_raman_spectrum.py &> read.out
#./4_plot_raman.py
