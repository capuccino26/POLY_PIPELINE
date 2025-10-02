#!/bin/bash
#$ -V
#$ -cwd
#$ -pe smp 2
#$ -l h_vmem=4G
#$ -N SETUP_STEREOPY_ENV
#$ -o setup_stereopy_$JOB_ID.out
#$ -e setup_stereopy_$JOB_ID.err

echo "==========================================="
echo "SETUP STEREOPY ENVIRONMENT"
echo "==========================================="
echo "START:"
date
echo "Job ID: $JOB_ID"
echo "Host: $HOSTNAME"
echo "Working directory: $(pwd)"
echo ""

# Load miniconda
echo "Loading miniconda3 module..."
module load miniconda3

echo "Conda version:"
conda --version
echo ""

# Check if env is already set
if conda env list | grep -q "^st "; then
    echo "WARNING: Environment 'st' already exists!"
    echo "Existing environments:"
    conda env list
    echo ""
    echo "If you want to recreate it, first remove with:"
    echo "conda env remove --name st"
    exit 0
fi

# Create conda env
echo "Creating conda environment 'st' with Python 3.8"
conda create --name st python=3.8 -y
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create conda environment"
    exit 1
fi
echo "Environment created successfully!"
echo ""

# Activate env
echo "Activating environment 'st'"
source activate st
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate conda environment"
    exit 1
fi
echo "Environment activated!"
echo ""

# Check python version
echo "Python version in environment:"
python --version
echo ""

# Install stereopy and dependencies
## STEREOPY REQUIRES SPECIFIC INSTALLATION PROCESS, CHECK DOCUMENTATION
echo "Installing Stereopy and dependencies, this may take several minutes"
conda install stereopy -c stereopy -c grst -c numba -c conda-forge -c bioconda -c fastai -c defaults -y

if [ $? -eq 0 ]; then
    echo "Stereopy installation completed successfully!"
else
    echo "ERROR: Stereopy installation failed!"
    echo "Trying alternative installation approach"
    
    # Try installing dependencies first
    echo "Installing base dependencies"
    conda install numpy pandas matplotlib seaborn scipy -c conda-forge -y
    
    # Try installing Stereopy again
    echo "Retrying Stereopy installation"
    conda install stereopy -c stereopy -c grst -c numba -c conda-forge -c bioconda -c fastai -c defaults -y
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to install Stereopy, check documentation for updates"
        conda deactivate
        exit 1
    fi
fi

echo ""

# Try installation
echo "Trying Stereopy installation"
python -c "
try:
    import stereo as st
    print('Stereopy successfully imported')
    print('Version:', st.__version__)
except ImportError as e:
    print('Failed to import Stereopy:', e)
    exit(1)

try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    print('Core dependencies successfully imported')
    print('pandas version:', pd.__version__)
    print('numpy version:', np.__version__)
except ImportError as e:
    print('Failed to import dependencies:', e)
    exit(1)

print('Environment is all set!')
"

if [ $? -eq 0 ]; then
    echo ""
    echo "==========================================="
    echo "Setup was completed sucessfully! To follow downstream analysis:"
    echo "  module load miniconda3"
    echo "  conda activate st"
    echo "==========================================="
    echo ""
else
    echo ""
    echo "==========================================="
    echo "SETUP FAILED! Check logs for errors"
    echo "==========================================="
    echo ""
fi

# Deactivate environment
conda deactivate

echo "==========================================="
echo "FINISHED:"
date
echo "==========================================="
