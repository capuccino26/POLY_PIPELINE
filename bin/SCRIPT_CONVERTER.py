#!/usr/bin/env python3
import stereo as st
import os
import warnings
import h5py
import numpy as np
import time
from datetime import datetime
import glob
import pandas as pd
import sys

# Ignore warnings
warnings.filterwarnings('ignore')

# Parameters
INPUT_PATH = "INPUT/RAW/S1-4_stereoseq.h5ad"
OUTPUT_DIR = "INPUT/RAW/CONVERTER_RESULTS"
BIN_SIZE = 100

# Function for logging stamps
def log_step(message):
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()

# Function to get bin files
def get_available_bins_from_file(file_path, file_ext):
    if file_ext == '.gem':
        return "User Defined"
    bins = []
    try:
        with h5py.File(file_path, 'r') as f:
            if 'geneExp' in f:
                bins = [b.replace('bin', '') for b in f['geneExp'].keys()]
            elif 'cellBin' in f:
                bins.append('cellBin')
            else:
                bins.append('1 (default)')
    except:
        return "Unknown"
    return ", ".join(bins) if bins else "Standard"

def generate_detailed_stats(data, output_dir, base_name, requested_bin):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    summary_file = os.path.join(output_dir, f"{base_name}_SUMMARY_STATS.txt")
    current_bin = getattr(data, 'bin_size', 1)
    
    with open(summary_file, "w") as f:
        f.write("="*100 + "\n")
        f.write(" " * 15 + "GEF FILE STATISTICS\n")
        f.write("="*100 + "\n\n")
        if str(requested_bin) != str(current_bin) and requested_bin != 50:
            f.write(f"Requested Bin: {requested_bin} | Actual Bin Used: {current_bin}\n")
            f.write(f"Structured files (H5AD/GEF) often have immutable bins.\n\n")
        f.write(f"1. Data Object Summary:\n{str(data)}\n\n")
        f.write(f"2. Genes Total: {len(data.gene_names)}\n")
        f.write(f"3. Cells Total: {len(data.cell_names)}\n")
        f.write(f"4. Bin Size Used: {current_bin}\n")
        f.write(f"5. Matrix Shape: {data.exp_matrix.shape}\n")
        
    pd.DataFrame(data.gene_names, columns=['gene_name']).to_csv(os.path.join(output_dir, f"{base_name}_GENES.csv"))
    pd.DataFrame(data.cell_names, columns=['cell_name']).to_csv(os.path.join(output_dir, f"{base_name}_CELLS.csv"))

def convert_file(path):
    start_time = time.time()
    file_ext = os.path.splitext(path)[1].lower()
    if path.endswith('.gem.gz'): file_ext = '.gem'

    base_name = os.path.basename(path).split('.')[0]
    output_path = os.path.join(OUTPUT_DIR, f"{base_name}.gef")
    avail = get_available_bins_from_file(path, file_ext)
    
    log_step(f"\n--- Processing: {os.path.basename(path)} ---")
    log_step(f"Format: {file_ext.upper()} | Available bins in input: {avail}")

    try:
        if file_ext == '.gem':
            data = st.io.read_gem(file_path=path, bin_type='bins', bin_size=BIN_SIZE, is_sparse=True)
        elif file_ext == '.h5ad':
            data = st.io.read_h5ad(path)
        else:
            log_step(f"Skipping {path}: Unsupported format.")
            return

        if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
        st.io.write_mid_gef(data=data, output=output_path)
        generate_detailed_stats(data, OUTPUT_DIR, base_name, BIN_SIZE)
        log_step(f"Done! Created {output_path} in {time.time() - start_time:.2f}s")
    except Exception as e:
        log_step(f"[ERROR] Failed to process {path}: {e}")

if __name__ == "__main__":
    # 1. Check if INPUT_PATH is empty
    if not INPUT_PATH:
        log_step("Module SCRIPT_CONVERTER skipped: No INPUT_PATH provided.")
        sys.exit(0)

    # 2. Check if Path exists
    if not os.path.exists(INPUT_PATH):
        log_step(f"Module SCRIPT_CONVERTER Error: Path '{INPUT_PATH}' not found.")
        sys.exit(1)

    # 3. Process File or Directory
    if os.path.isfile(INPUT_PATH):
        convert_file(INPUT_PATH)
    elif os.path.isdir(INPUT_PATH):
        files = glob.glob(os.path.join(INPUT_PATH, "*.gem*")) + glob.glob(os.path.join(INPUT_PATH, "*.h5ad"))
        if not files:
            log_step(f"No compatible files (GEM/H5AD) found in {INPUT_PATH}")
        for f in files:
            convert_file(f)
