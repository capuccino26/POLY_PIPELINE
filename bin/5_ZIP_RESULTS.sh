#!/bin/bash
#$ -V
#$ -cwd
#$ -pe smp 1
#$ -N ZIP_RESULTS_STEREOPY

echo "==========================================="
echo "STEREOPY ZIP RESULTS"
echo "==========================================="
echo "START:"
date
echo "Job ID: $JOB_ID"
echo "Host: $HOSTNAME"
echo "Working directory: $(pwd)"
echo "CPUs allocated: $NSLOTS"
echo ""

# System resources check
echo "System Resources Available:"
free -h
echo "CPU Info:"
nproc
lscpu | grep "Model name"
echo ""

current_date=$(date +%Y%m%d)
zip_filename="RESULTS/RESULTS_${current_date}.zip"
zip -r "${zip_filename}" RESULTS/results_annotation RESULTS/results_ultimate RESULTS/results_validation

echo "==========================================="
echo "FINISHED:"
date
echo "==========================================="
date
