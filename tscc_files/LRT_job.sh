#!/bin/bash

#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -o /projects/ps-gymreklab/bonnieh/lrt/logs/0810/jobname.out
#PBS -e /projects/ps-gymreklab/bonnieh/lrt/logs/0810/jobname.err
#PBS -A gymreklab-group

# Running command
time python /projects/ps-gymreklab/bonnieh/lrt/LRT_lookup.py per opt filenum s_vals num_sims
