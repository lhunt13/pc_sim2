#!/bin/bash
#$ -l mem_free=60G
#$ -l h_vmem=60G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-400

mkdir -p results
if [ ! -f "/results/run-$SGE_TASK_ID" ]; then
  Rscript simulation_study.R $SGE_TASK_ID 
fi