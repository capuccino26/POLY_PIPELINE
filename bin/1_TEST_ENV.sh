#!/bin/bash
#$ -V
#$ -cwd
#$ -pe smp 2
#$ -N TEST_STEREOPY
#$ -o test_stereopy_$JOB_ID.out
#$ -e test_stereopy_$JOB_ID.err

echo "==========================================="
echo "TEST STEREOPY ENVIRONMENT"
echo "==========================================="
echo "START:"
date
echo "Job ID: $JOB_ID"
echo "Host: $HOSTNAME"
echo "Working directory: $(pwd)"
echo ""

if [ -z "$ST_PYTHON" ] || [ ! -f "$ST_PYTHON" ]; then
    echo "=========================================================" >&2
    echo "ERROR: The ST_PYTHON variable was not defined during submission." >&2
    echo "=========================================================" >&2

    echo "========================================================="
    echo "ERROR: The ST_PYTHON variable was not defined during submission, check error log for more information!"
    echo "========================================================="
    
    EXAMPLE_PATH="/home/user/.conda/envs/st/bin/python"
    
    echo "1. Example Format:" >&2
    echo "   The expected value is the full path to the 'python' executable inside your 'st' conda environment." >&2
    echo "   Example: $EXAMPLE_PATH" >&2
    echo "" >&2

    echo "2. How to Find the Correct Path:" >&2
    echo "   If you don't know the exact path, run the following steps in your terminal (outside the qsub submission):" >&2
    echo "   a) Load the Conda module: 'module load miniconda3' (or similar)" >&2
    echo "   b) List the Conda environments: 'conda env list' " >&2
    echo "      (Look for the line listing the 'st' environment; it will show the directory path.)" >&2
    echo "      (If you don't find the 'st' environment run the script 'bin/0_SET_ENV.sh' to correctly set the environment.)" >&2
    echo "   c) The ST_PYTHON path will be: [Conda st Environment Path]/bin/python" >&2
    echo "" >&2
    echo "3. Correct Submission (SGE/qsub):" >&2
    echo "   qsub -v ST_PYTHON=/home/user/.conda/envs/st/bin/python bin/1_TEST_ENV.sh" >&2
    
    exit 1
fi

echo "Path for testing environment: $ST_PYTHON"
echo ""

echo "Testing Stereopy environment with explicit paths"
module load miniconda3

echo "Available environments:"
conda env list

echo "Using explicit Python path: $ST_PYTHON"
echo "Python version:"
$ST_PYTHON --version

echo "Testing Stereopy import:"
$ST_PYTHON -c "import stereo as st; print('Stereopy version:', st.__version__)"

echo "Testing other dependencies:"
$ST_PYTHON -c "import pandas as pd; print('Pandas version:', pd.__version__)"

echo "==========================================="
echo "FINISHED:"
date
echo "==========================================="
