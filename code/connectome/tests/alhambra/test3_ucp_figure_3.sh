#!/bin/bash
#
#$ -S /bin/bash
#$ -N test3_ucp_figure_3
#$ -o test3_ucp_figure_3.out
#$ -j y
#$ -cwd
#$ -l NOParalela
#$ -M dsalazar@ugr.es
#$ -m b
#$ -m e
#$ -m a
#$ -m s
export R_LIBS=/home/UGR002/dsalazar/RLib

module load alhambra/R-3.2.1 #Si vas a usar R
module load alhambra/gsl-1.16
module load alhambra/fftw-3.3.3
R CMD BATCH test3_ucp_figure_3.R
test3_ucp_figure_3.Rout

