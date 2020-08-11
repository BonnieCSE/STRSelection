#!/bin/bash

#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/bonnieh/abc/logs/0810/jobname.out
#PBS -e /projects/ps-gymreklab/bonnieh/abc/logs/0810/jobname.err
#PBS -A gymreklab-group

PERIOD=per
OPTIMAL_ALLELE=opt
K_PARAM=k_param
THETA_PARAM=theta_param
FILENUM=filenum
OUT_FOLDER=out_folder
# Running command
time python /projects/ps-gymreklab/bonnieh/abc/ABC_lookup.py per opt num_sims k_param theta_param filenum out_folder
