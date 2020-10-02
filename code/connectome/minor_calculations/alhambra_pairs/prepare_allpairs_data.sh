#!/bin/bash
#
#$ -S /bin/bash
#$ -N prepare_allpairs_data
#$ -o prepare_allpairs_data.out
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
R CMD BATCH prepare_allpairs_data.R
prepare_allpairs_data.Rout
