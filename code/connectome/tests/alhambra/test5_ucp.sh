#!/bin/bash
#
#$ -S /bin/bash
#$ -N test5_ucp
#$ -o test5_ucp.out
#$ -j y
#$ -cwd
#$ -l NOParalela
#$ -M dsalazar@ugr.es
#$ -m b
#$ -m e
#$ -m a
#$ -m s
export R_LIBS=/home/UGR002/dsalazar/RLib

module load alhambra/R-3.3.2 #Si vas a usar R
module load alhambra/gsl-1.16
module load alhambra/fftw-3.3.3
R CMD BATCH test5_ucp.R
test5_ucp.Rout
