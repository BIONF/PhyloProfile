#!/bin/sh
#PBS -N qsub_R_test

cd $PBS_O_WORKDIR
echo installing shinyPackage ...

R CMD BATCH install.R