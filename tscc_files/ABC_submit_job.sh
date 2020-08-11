#!/bin/bash
file="/projects/ps-gymreklab/bonnieh/abc/ABC_combos.txt"
while IFS=$'\t' read -r -a line; do
    PERIOD=${line[0]}
    OPTIMAL_ALLELE=${line[1]}
    K_PARAM=${line[2]}
    THETA_PARAM=${line[3]}
    FILENUM=${line[4]}
    OUT_FOLDER=${line[5]}
    JOBNAME="${PERIOD}_${OPTIMAL_ALLELE}_${OUT_FOLDER}_${FILENUM}"
    NUM_SIMS=${line[6]}
    
    
    echo Running ABC for $JOBNAME 
    cat ABC_job.sh | sed "s/jobname/${JOBNAME}/g" | sed "s#per#${PERIOD}#g" | sed "s#opt#${OPTIMAL_ALLELE}#g" |  sed "s#num_sims#${NUM_SIMS}#g" | sed "s#k_param#${K_PARAM}#g" | sed "s#theta_param#${THETA_PARAM}#g" | sed "s#filenum#${FILENUM}#g" | sed "s/out_folder/${OUT_FOLDER}/g" |qsub
done < "$file"
