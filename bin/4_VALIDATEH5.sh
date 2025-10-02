#!/bin/bash
#$ -V
#$ -cwd
#$ -pe smp 8
#$ -l h_vmem=8G
#$ -N STEREOPY_DATA_VALIDATION
#$ -o stereopy_validation_$JOB_ID.out
#$ -e stereopy_validation_$JOB_ID.err

echo "==========================================="
echo "STEREOPY DATA VALIDATION ANALYSIS"
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

# Check if the variable was correctly defined
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
    echo "   qsub -v ST_PYTHON=\"/your/full/path/to/st/bin/python\" bin/1_TEST_ENV.sh" >&2
    
    exit 1
fi

# Load environment with explicit paths
echo "Loading miniconda3 module..."
module load miniconda3

ST_BIN_DIR=$(dirname "$ST_PYTHON")
export PATH="$ST_BIN_DIR:$PATH"

echo "Using Python: $ST_PYTHON"
$ST_PYTHON --version

echo "Verifying dependencies..."
$ST_PYTHON -c "
import sys
print(f'Python executable: {sys.executable}')
import stereo as st
print(f'Stereopy version: {st.__version__}')
import pandas as pd
print(f'Pandas version: {pd.__version__}')
import numpy as np
print(f'NumPy version: {np.__version__}')
import anndata as ad
print(f'AnnData version: {ad.__version__}')
import h5py
print(f'H5PY version: {h5py.__version__}')
"
echo ""

# Create the validation script
echo "Creating data validation script..."
cat > bin/data_validation_analysis.py << EOF
#!/usr/bin/env python3

import stereo as st
import warnings
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import psutil
import gc
from datetime import datetime
import sys
import json
import h5py
import anndata as ad
from pathlib import Path

# Set matplotlib backend for cluster
import matplotlib
matplotlib.use('Agg')

warnings.filterwarnings('ignore')

def log_memory_usage(step_name=""):
    """Enhanced memory monitoring"""
    try:
        memory = psutil.virtual_memory()
        print(f"[MEMORY {step_name}] RAM: {memory.percent:.1f}% ({memory.used/1e9:.1f}GB/{memory.total/1e9:.1f}GB)")
    except:
        print(f"[MEMORY {step_name}] Unable to get memory info")

def log_step(message):
    """Enhanced logging with timestamp"""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()

def get_file_size_info(file_path):
    """Get comprehensive file size information"""
    if not os.path.exists(file_path):
        return None
    
    size_bytes = os.path.getsize(file_path)
    size_mb = size_bytes / 1e6
    size_gb = size_bytes / 1e9
    
    return {
        'bytes': size_bytes,
        'mb': size_mb,
        'gb': size_gb,
        'exists': True
    }

def analyze_gef_file(gef_path, bin_size=100):
    """Analyze original GEF file structure and content"""
    log_step(f"Analyzing original GEF file: {gef_path}")
    
    try:
        # Get GEF file info first
        gef_info = st.io.read_gef_info(gef_path)
        log_step("GEF info retrieved successfully")
        
        # Load the data
        data = st.io.read_gef(file_path=gef_path, bin_size=bin_size)
        log_step(f"GEF data loaded: {data}")
        
        # Extract comprehensive statistics
        gef_stats = {
            'file_path': gef_path,
            'bin_size': bin_size,
            'n_cells': data.n_cells,
            'n_genes': data.n_genes,
            'data_shape': (data.n_cells, data.n_genes),
            'file_info': get_file_size_info(gef_path)
        }
        
        # Get spatial information if available
        if hasattr(data.cells, 'spatial'):
            spatial_coords = data.cells.spatial
            gef_stats.update({
                'spatial_dims': spatial_coords.shape,
                'spatial_x_range': (float(spatial_coords[:, 0].min()), float(spatial_coords[:, 0].max())),
                'spatial_y_range': (float(spatial_coords[:, 1].min()), float(spatial_coords[:, 1].max())),
                'spatial_mean_x': float(spatial_coords[:, 0].mean()),
                'spatial_mean_y': float(spatial_coords[:, 1].mean())
            })
        
        # Get gene information
        if hasattr(data.genes, 'gene_name'):
            gef_stats['gene_names_sample'] = data.genes.gene_name[:10].tolist()
            gef_stats['total_unique_genes'] = len(data.genes.gene_name)
        
        # Get cell information  
        if hasattr(data.cells, 'cell_name'):
            gef_stats['cell_names_sample'] = data.cells.cell_name[:10].tolist()
            gef_stats['total_unique_cells'] = len(data.cells.cell_name)
        
        # Expression data statistics
        if hasattr(data, 'exp_matrix'):
            exp_matrix = data.exp_matrix
            if hasattr(exp_matrix, 'toarray'):
                # Sparse matrix - sample for stats
                sample_data = exp_matrix[:1000, :100].toarray() if exp_matrix.shape[0] > 1000 else exp_matrix.toarray()
            else:
                sample_data = exp_matrix[:1000, :100] if exp_matrix.shape[0] > 1000 else exp_matrix
            
            gef_stats.update({
                'expression_matrix_shape': exp_matrix.shape,
                'expression_matrix_type': str(type(exp_matrix)),
                'is_sparse': hasattr(exp_matrix, 'toarray'),
                'sample_mean_expression': float(np.mean(sample_data)),
                'sample_max_expression': float(np.max(sample_data)),
                'sample_min_expression': float(np.min(sample_data)),
                'sample_nonzero_fraction': float(np.count_nonzero(sample_data) / sample_data.size)
            })
        
        log_step("GEF analysis completed successfully")
        return gef_stats, data
        
    except Exception as e:
        log_step(f"ERROR analyzing GEF file: {e}")
        return None, None

def analyze_h5ad_file(h5ad_path):
    """Analyze H5AD file structure and content"""
    log_step(f"Analyzing H5AD file: {h5ad_path}")
    
    try:
        # Load AnnData object
        adata = ad.read_h5ad(h5ad_path)
        log_step(f"H5AD data loaded: {adata.n_obs} cells × {adata.n_vars} genes")
        
        # Extract comprehensive statistics
        h5ad_stats = {
            'file_path': h5ad_path,
            'n_obs': adata.n_obs,
            'n_vars': adata.n_vars, 
            'data_shape': (adata.n_obs, adata.n_vars),
            'file_info': get_file_size_info(h5ad_path)
        }
        
        # Observation (cell) information
        h5ad_stats.update({
            'obs_columns': list(adata.obs.columns),
            'obs_dtypes': {col: str(dtype) for col, dtype in adata.obs.dtypes.items()},
            'obs_sample': adata.obs.head(3).to_dict() if len(adata.obs) > 0 else {}
        })
        
        # Variable (gene) information
        h5ad_stats.update({
            'var_columns': list(adata.var.columns),
            'var_dtypes': {col: str(dtype) for col, dtype in adata.var.dtypes.items()},
            'var_sample': adata.var.head(3).to_dict() if len(adata.var) > 0 else {}
        })
        
        # Unstructured data
        h5ad_stats['uns_keys'] = list(adata.uns.keys())
        
        # Multidimensional observations
        h5ad_stats['obsm_keys'] = list(adata.obsm.keys())
        if 'spatial' in adata.obsm:
            spatial_coords = adata.obsm['spatial']
            h5ad_stats.update({
                'spatial_dims': spatial_coords.shape,
                'spatial_x_range': (float(spatial_coords[:, 0].min()), float(spatial_coords[:, 0].max())),
                'spatial_y_range': (float(spatial_coords[:, 1].min()), float(spatial_coords[:, 1].max())),
                'spatial_mean_x': float(spatial_coords[:, 0].mean()),
                'spatial_mean_y': float(spatial_coords[:, 1].mean())
            })
        
        # Pairwise observations
        h5ad_stats['obsp_keys'] = list(adata.obsp.keys())
        
        # Expression data statistics
        X_data = adata.X
        if hasattr(X_data, 'toarray'):
            # Sparse matrix - sample for stats
            sample_data = X_data[:1000, :100].toarray() if X_data.shape[0] > 1000 else X_data.toarray()
        else:
            sample_data = X_data[:1000, :100] if X_data.shape[0] > 1000 else X_data
        
        h5ad_stats.update({
            'expression_matrix_shape': X_data.shape,
            'expression_matrix_type': str(type(X_data)),
            'is_sparse': hasattr(X_data, 'toarray'),
            'sample_mean_expression': float(np.mean(sample_data)),
            'sample_max_expression': float(np.max(sample_data)),
            'sample_min_expression': float(np.min(sample_data)),
            'sample_nonzero_fraction': float(np.count_nonzero(sample_data) / sample_data.size)
        })
        
        # Check for raw data
        if adata.raw is not None:
            h5ad_stats['has_raw_data'] = True
            h5ad_stats['raw_shape'] = adata.raw.shape
        else:
            h5ad_stats['has_raw_data'] = False
        
        # Clustering information
        clustering_columns = [col for col in adata.obs.columns if any(x in col.lower() for x in ['cluster', 'louvain', 'leiden'])]
        if clustering_columns:
            h5ad_stats['clustering_methods'] = clustering_columns
            for cluster_col in clustering_columns:
                h5ad_stats[f'{cluster_col}_unique_values'] = len(adata.obs[cluster_col].unique())
        
        # Annotation information
        annotation_columns = [col for col in adata.obs.columns if 'anno' in col.lower()]
        if annotation_columns:
            h5ad_stats['annotation_methods'] = annotation_columns
            for anno_col in annotation_columns:
                h5ad_stats[f'{anno_col}_unique_values'] = len(adata.obs[anno_col].unique())
        
        log_step("H5AD analysis completed successfully")
        return h5ad_stats, adata
        
    except Exception as e:
        log_step(f"ERROR analyzing H5AD file: {e}")
        return None, None

def compare_datasets(gef_stats, h5ad_stats):
    """Compare GEF and H5AD datasets for validation"""
    log_step("Performing comprehensive dataset comparison...")
    
    comparison_results = {
        'validation_timestamp': datetime.now().isoformat(),
        'overall_status': 'UNKNOWN'
    }
    
    validation_checks = []
    
    # Check 1: Cell count comparison
    gef_cells = gef_stats.get('n_cells', 0)
    h5ad_cells = h5ad_stats.get('n_obs', 0)
    
    cell_check = {
        'check_name': 'cell_count_validation',
        'gef_value': gef_cells,
        'h5ad_value': h5ad_cells,
        'match': gef_cells == h5ad_cells,
        'difference': abs(gef_cells - h5ad_cells),
        'relative_diff_percent': abs(gef_cells - h5ad_cells) / max(gef_cells, 1) * 100
    }
    validation_checks.append(cell_check)
    
    # Check 2: Gene count comparison (may differ due to filtering)
    gef_genes = gef_stats.get('n_genes', 0)
    h5ad_genes = h5ad_stats.get('n_vars', 0)
    
    gene_check = {
        'check_name': 'gene_count_validation',
        'gef_value': gef_genes,
        'h5ad_value': h5ad_genes,
        'match': gef_genes == h5ad_genes,
        'difference': abs(gef_genes - h5ad_genes),
        'relative_diff_percent': abs(gef_genes - h5ad_genes) / max(gef_genes, 1) * 100,
        'note': 'Gene counts may differ due to filtering during preprocessing'
    }
    validation_checks.append(gene_check)
    
    # Check 3: Spatial coordinates comparison
    spatial_check = {
        'check_name': 'spatial_coordinates_validation',
        'gef_has_spatial': 'spatial_x_range' in gef_stats,
        'h5ad_has_spatial': 'spatial_x_range' in h5ad_stats
    }
    
    if spatial_check['gef_has_spatial'] and spatial_check['h5ad_has_spatial']:
        gef_x_range = gef_stats['spatial_x_range']
        h5ad_x_range = h5ad_stats['spatial_x_range']
        gef_y_range = gef_stats['spatial_y_range']
        h5ad_y_range = h5ad_stats['spatial_y_range']
        
        spatial_check.update({
            'x_range_match': np.allclose(gef_x_range, h5ad_x_range, rtol=1e-5),
            'y_range_match': np.allclose(gef_y_range, h5ad_y_range, rtol=1e-5),
            'gef_x_range': gef_x_range,
            'h5ad_x_range': h5ad_x_range,
            'gef_y_range': gef_y_range,
            'h5ad_y_range': h5ad_y_range
        })
    
    validation_checks.append(spatial_check)
    
    # Check 4: Expression data comparison (sample-based)
    expression_check = {
        'check_name': 'expression_data_validation',
        'gef_mean': gef_stats.get('sample_mean_expression', 0),
        'h5ad_mean': h5ad_stats.get('sample_mean_expression', 0),
        'gef_max': gef_stats.get('sample_max_expression', 0),
        'h5ad_max': h5ad_stats.get('sample_max_expression', 0),
        'gef_nonzero_frac': gef_stats.get('sample_nonzero_fraction', 0),
        'h5ad_nonzero_frac': h5ad_stats.get('sample_nonzero_fraction', 0)
    }
    
    # Note: Expression values will differ due to normalization/transformation
    expression_check['note'] = 'Expression values expected to differ due to normalization and log transformation'
    validation_checks.append(expression_check)
    
    # Check 5: File size comparison
    gef_size = gef_stats.get('file_info', {}).get('mb', 0)
    h5ad_size = h5ad_stats.get('file_info', {}).get('mb', 0)
    
    size_check = {
        'check_name': 'file_size_comparison',
        'gef_size_mb': gef_size,
        'h5ad_size_mb': h5ad_size,
        'size_ratio': h5ad_size / max(gef_size, 1),
        'note': 'H5AD typically smaller due to compression and processing'
    }
    validation_checks.append(size_check)
    
    # Overall validation status
    critical_checks_passed = (
        cell_check['match'] and
        spatial_check.get('x_range_match', True) and
        spatial_check.get('y_range_match', True)
    )
    
    if critical_checks_passed:
        comparison_results['overall_status'] = 'PASSED'
    else:
        comparison_results['overall_status'] = 'FAILED'
    
    comparison_results['validation_checks'] = validation_checks
    comparison_results['summary'] = {
        'total_checks': len(validation_checks),
        'critical_checks_passed': critical_checks_passed,
        'cell_count_preserved': cell_check['match'],
        'spatial_coordinates_preserved': spatial_check.get('x_range_match', False) and spatial_check.get('y_range_match', False)
    }
    
    return comparison_results

# Function for graphical visualization of comparisons
def create_validation_visualizations(gef_stats, h5ad_stats, comparison_results, output_dir):
    log_step("Creating validation visualization plots")
    
    try:
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # Plot 1: Data dimensions comparison
        dimensions = ['Cells/Obs', 'Genes/Vars']
        gef_values = [gef_stats.get('n_cells', 0), gef_stats.get('n_genes', 0)]
        h5ad_values = [h5ad_stats.get('n_obs', 0), h5ad_stats.get('n_vars', 0)]
        
        x = np.arange(len(dimensions))
        width = 0.35
        
        axes[0, 0].bar(x - width/2, gef_values, width, label='GEF Original', alpha=0.8)
        axes[0, 0].bar(x + width/2, h5ad_values, width, label='H5AD Export', alpha=0.8)
        axes[0, 0].set_xlabel('Data Dimensions')
        axes[0, 0].set_ylabel('Count')
        axes[0, 0].set_title('Data Dimensions Comparison')
        axes[0, 0].set_xticks(x)
        axes[0, 0].set_xticklabels(dimensions)
        axes[0, 0].legend()
        
        # Add value labels on bars
        for i, v in enumerate(gef_values):
            axes[0, 0].text(i - width/2, v + max(gef_values) * 0.01, str(v), ha='center')
        for i, v in enumerate(h5ad_values):
            axes[0, 0].text(i + width/2, v + max(h5ad_values) * 0.01, str(v), ha='center')
        
        # Plot 2: File size comparison
        file_sizes = [
            gef_stats.get('file_info', {}).get('mb', 0),
            h5ad_stats.get('file_info', {}).get('mb', 0)
        ]
        file_labels = ['GEF Original', 'H5AD Export']
        
        bars = axes[0, 1].bar(file_labels, file_sizes, color=['skyblue', 'lightcoral'], alpha=0.8)
        axes[0, 1].set_ylabel('File Size (MB)')
        axes[0, 1].set_title('File Size Comparison')
        
        # Add value labels
        for bar, size in zip(bars, file_sizes):
            axes[0, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(file_sizes) * 0.01, 
                           f'{size:.1f} MB', ha='center')
        
        # Plot 3: Spatial range comparison (if available)
        if 'spatial_x_range' in gef_stats and 'spatial_x_range' in h5ad_stats:
            ranges_data = {
                'X Min': [gef_stats['spatial_x_range'][0], h5ad_stats['spatial_x_range'][0]],
                'X Max': [gef_stats['spatial_x_range'][1], h5ad_stats['spatial_x_range'][1]],
                'Y Min': [gef_stats['spatial_y_range'][0], h5ad_stats['spatial_y_range'][0]],
                'Y Max': [gef_stats['spatial_y_range'][1], h5ad_stats['spatial_y_range'][1]]
            }
            
            x = np.arange(len(ranges_data))
            gef_ranges = [ranges_data[key][0] for key in ranges_data]
            h5ad_ranges = [ranges_data[key][1] for key in ranges_data]
            
            axes[0, 2].bar(x - width/2, gef_ranges, width, label='GEF Original', alpha=0.8)
            axes[0, 2].bar(x + width/2, h5ad_ranges, width, label='H5AD Export', alpha=0.8)
            axes[0, 2].set_xlabel('Spatial Coordinates')
            axes[0, 2].set_ylabel('Coordinate Value')
            axes[0, 2].set_title('Spatial Range Comparison')
            axes[0, 2].set_xticks(x)
            axes[0, 2].set_xticklabels(list(ranges_data.keys()))
            axes[0, 2].legend()
        else:
            axes[0, 2].text(0.5, 0.5, 'Spatial data\nnot available', 
                           ha='center', va='center', transform=axes[0, 2].transAxes)
            axes[0, 2].set_title('Spatial Range Comparison')
        
        # Plot 4: Expression data comparison
        expr_metrics = ['Mean', 'Max', 'NonZero Fraction']
        gef_expr = [
            gef_stats.get('sample_mean_expression', 0),
            gef_stats.get('sample_max_expression', 0),
            gef_stats.get('sample_nonzero_fraction', 0)
        ]
        h5ad_expr = [
            h5ad_stats.get('sample_mean_expression', 0),
            h5ad_stats.get('sample_max_expression', 0),
            h5ad_stats.get('sample_nonzero_fraction', 0)
        ]
        
        x = np.arange(len(expr_metrics))
        axes[1, 0].bar(x - width/2, gef_expr, width, label='GEF Original', alpha=0.8)
        axes[1, 0].bar(x + width/2, h5ad_expr, width, label='H5AD Export', alpha=0.8)
        axes[1, 0].set_xlabel('Expression Metrics')
        axes[1, 0].set_ylabel('Value')
        axes[1, 0].set_title('Expression Data Comparison')
        axes[1, 0].set_xticks(x)
        axes[1, 0].set_xticklabels(expr_metrics)
        axes[1, 0].legend()
        axes[1, 0].text(0.5, 0.95, 'Note: Values differ due to normalization', 
                       ha='center', transform=axes[1, 0].transAxes, fontsize=8)
        
        # Plot 5: Validation checks summary
        check_names = []
        check_status = []
        for check in comparison_results['validation_checks']:
            check_names.append(check['check_name'].replace('_', '\n'))
            if 'match' in check:
                check_status.append(1 if check['match'] else 0)
            elif 'x_range_match' in check and 'y_range_match' in check:
                check_status.append(1 if (check.get('x_range_match', False) and check.get('y_range_match', False)) else 0)
            else:
                check_status.append(0.5)  # Neutral for info-only checks
        
        colors = ['green' if s == 1 else 'red' if s == 0 else 'yellow' for s in check_status]
        axes[1, 1].bar(range(len(check_names)), check_status, color=colors, alpha=0.7)
        axes[1, 1].set_xlabel('Validation Checks')
        axes[1, 1].set_ylabel('Pass Status')
        axes[1, 1].set_title('Validation Checks Summary')
        axes[1, 1].set_xticks(range(len(check_names)))
        axes[1, 1].set_xticklabels(check_names, rotation=45, ha='right')
        axes[1, 1].set_ylim(0, 1.2)
        
        # Plot 6: Overall validation status
        status_text = f"Overall Status: {comparison_results['overall_status']}"
        status_color = 'green' if comparison_results['overall_status'] == 'PASSED' else 'red'
        
        axes[1, 2].text(0.5, 0.7, status_text, ha='center', va='center', 
                       transform=axes[1, 2].transAxes, fontsize=16, fontweight='bold',
                       bbox=dict(boxstyle="round,pad=0.3", facecolor=status_color, alpha=0.3))
        
        summary = comparison_results['summary']
        summary_text = f"""
Cells Preserved: {'Pass' if summary['cell_count_preserved'] else 'Fail'}
Spatial Coords: {'Pass' if summary['spatial_coordinates_preserved'] else 'Fail'}
Critical Checks: {'Pass' if summary['critical_checks_passed'] else 'Fail'}
        """
        axes[1, 2].text(0.5, 0.3, summary_text, ha='center', va='center',
                       transform=axes[1, 2].transAxes, fontsize=12,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.3))
        axes[1, 2].set_title('Validation Summary')
        axes[1, 2].axis('off')
        
        plt.tight_layout()
        validation_plot_file = os.path.join(output_dir, 'comprehensive_validation_analysis.png')
        plt.savefig(validation_plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        log_step(f"Validation visualization saved to {validation_plot_file}")
        
    except Exception as e:
        log_step(f"Warning: Error creating validation plots: {e}")

# MAIN VALIDATION ANALYSIS

start_time = datetime.now()
log_step("Starting data validation analysis")
log_memory_usage("START")

# Set up validation output directory
validation_output_dir = 'RESULTS/results_validation'
os.makedirs(validation_output_dir, exist_ok=True)

# Function to check the .gef file input
def get_single_gef_file(datasets_dir='INPUT/datasets'):
    search_pattern = os.path.join(datasets_dir, '*.gef')
    if not os.path.isdir(datasets_dir):
        print("=" * 55, file=sys.stderr)
        print(f"FATAL ERROR: The required input directory '{datasets_dir}' was not found.", file=sys.stderr)
        print("Please create this directory and place your .gef file inside.", file=sys.stderr)
        print("=" * 55, file=sys.stderr)
        sys.exit(1)
    gef_files = glob.glob(search_pattern)
    num_gef_files = len(gef_files)
    if num_gef_files == 0:
        print("=" * 55, file=sys.stderr)
        print(f"FATAL ERROR: No .gef file found in the '{datasets_dir}' directory.", file=sys.stderr)
        print("Please ensure there is exactly one .gef file to be analyzed.", file=sys.stderr)
        print("=" * 55, file=sys.stderr)
        sys.exit(1)
    if num_gef_files > 1:
        print("=" * 55, file=sys.stderr)
        print(f"FATAL ERROR: Multiple .gef files found in '{datasets_dir}'.", file=sys.stderr)
        print("Files found:", file=sys.stderr)
        for file_path in gef_files:
            print(f"  - {os.path.basename(file_path)}", file=sys.stderr)
        print("\nACTION REQUIRED: Please ensure ONLY the .gef file intended for analysis is in the directory.", file=sys.stderr)
        print("=" * 55, file=sys.stderr)
        sys.exit(1)
    data_path = gef_files[0]
    print(f"Successfully detected input file: {os.path.basename(data_path)}")
    print(f"Data path set to: {data_path}")
    
    return data_path

# File paths
data_path = get_single_gef_file()
h5ad_search_dir = 'results_annotation/exports/h5ad'

# Find the most recent H5AD file
h5ad_files = []
if os.path.exists(h5ad_search_dir):
    for f in os.listdir(h5ad_search_dir):
        if f.endswith('.h5ad'):
            full_path = os.path.join(h5ad_search_dir, f)
            h5ad_files.append((full_path, os.path.getmtime(full_path)))

if not h5ad_files:
    log_step("ERROR: No H5AD files found in results_annotation/exports/h5ad/")
    log_step("Please run the annotation analysis first (bin/3_ANNOTATION_ANALYSIS.sh)")
    sys.exit(1)

# Get the most recent H5AD file
h5ad_path = sorted(h5ad_files, key=lambda x: x[1], reverse=True)[0][0]
log_step(f"Found H5AD file: {h5ad_path}")

# Verify files exist
if not os.path.exists(data_path):
    log_step(f"ERROR: GEF file not found at {data_path}")
    sys.exit(1)

if not os.path.exists(h5ad_path):
    log_step(f"ERROR: H5AD file not found at {h5ad_path}")
    sys.exit(1)

# Analyze original GEF file
log_step("="*80)
log_step("ANALYZING ORIGINAL GEF FILE")
log_step("="*80)
log_memory_usage("PRE_GEF_ANALYSIS")

gef_stats, gef_data = analyze_gef_file(data_path, bin_size=100)
if gef_stats is None:
    log_step("FATAL ERROR: Could not analyze GEF file")
    sys.exit(1)

log_memory_usage("POST_GEF_ANALYSIS")
gc.collect()

# Analyze H5AD file
log_step("="*80)
log_step("ANALYZING H5AD EXPORT FILE")
log_step("="*80)
log_memory_usage("PRE_H5AD_ANALYSIS")

h5ad_stats, h5ad_data = analyze_h5ad_file(h5ad_path)
if h5ad_stats is None:
    log_step("FATAL ERROR: Could not analyze H5AD file")
    sys.exit(1)

log_memory_usage("POST_H5AD_ANALYSIS")

log_step("="*80)
log_step("PERFORMING VALIDATION")
log_step("="*80)

comparison_results = compare_datasets(gef_stats, h5ad_stats)

# Create validation visualizations
create_validation_visualizations(gef_stats, h5ad_stats, comparison_results, validation_output_dir)

# Save detailed validation results
validation_results_file = os.path.join(validation_output_dir, 'detailed_validation_results.json')
full_validation_data = {
    'gef_analysis': gef_stats,
    'h5ad_analysis': h5ad_stats,
    'comparison_results': comparison_results
}

with open(validation_results_file, 'w') as f:
    json.dump(full_validation_data, f, indent=2, default=str)

log_step(f"Detailed validation results saved to {validation_results_file}")

# Generate validation report
log_step("="*80)
log_step("GENERATING VALIDATION REPORT")
log_step("="*80)

end_time = datetime.now()
total_duration = (end_time - start_time).total_seconds()

validation_report_file = os.path.join(validation_output_dir, 'VALIDATION_REPORT.txt')
with open(validation_report_file, 'w') as f:
    f.write("="*100 + "\n")
    f.write("GEF vs H5AD Export Comparison\n")
    f.write("="*100 + "\n\n")
    
    f.write(f"Validation completed: {end_time}\n")
    f.write(f"Processing time: {total_duration:.2f} seconds\n")
    f.write(f"Overall validation status: {comparison_results['overall_status']}\n\n")
    
    f.write("FILE INFORMATION:\n")
    f.write(f"Original GEF file: {data_path}\n")
    f.write(f"GEF file size: {gef_stats['file_info']['mb']:.1f} MB\n")
    f.write(f"H5AD file: {os.path.basename(h5ad_path)}\n")
    f.write(f"H5AD file size: {h5ad_stats['file_info']['mb']:.1f} MB\n")
    f.write(f"Size compression ratio: {h5ad_stats['file_info']['mb'] / max(gef_stats['file_info']['mb'], 1):.2f}\n\n")
    
    f.write("DATA DIMENSIONS:\n")
    f.write(f"Original GEF: {gef_stats['n_cells']} cells × {gef_stats['n_genes']} genes\n")
    f.write(f"H5AD Export: {h5ad_stats['n_obs']} cells × {h5ad_stats['n_vars']} genes\n")
    f.write(f"Cell preservation: {'✓ PASS' if gef_stats['n_cells'] == h5ad_stats['n_obs'] else '✗ FAIL'}\n")
    f.write(f"Gene preservation: {'✓ PASS' if gef_stats['n_genes'] == h5ad_stats['n_vars'] else '⚠ DIFFER (expected due to filtering)'}\n\n")
    
    if 'spatial_x_range' in gef_stats and 'spatial_x_range' in h5ad_stats:
        f.write("SPATIAL COORDINATES:\n")
        f.write(f"GEF X range: {gef_stats['spatial_x_range'][0]:.1f} to {gef_stats['spatial_x_range'][1]:.1f}\n")
        f.write(f"H5AD X range: {h5ad_stats['spatial_x_range'][0]:.1f} to {h5ad_stats['spatial_x_range'][1]:.1f}\n")
        f.write(f"GEF Y range: {gef_stats['spatial_y_range'][0]:.1f} to {gef_stats['spatial_y_range'][1]:.1f}\n")
        f.write(f"H5AD Y range: {h5ad_stats['spatial_y_range'][0]:.1f} to {h5ad_stats['spatial_y_range'][1]:.1f}\n")
        
        x_match = np.allclose(gef_stats['spatial_x_range'], h5ad_stats['spatial_x_range'], rtol=1e-5)
        y_match = np.allclose(gef_stats['spatial_y_range'], h5ad_stats['spatial_y_range'], rtol=1e-5)
        f.write(f"Spatial coordinates preserved: {'✓ PASS' if x_match and y_match else '✗ FAIL'}\n\n")
    
    f.write("H5AD SPECIFIC INFORMATION:\n")
    f.write(f"Observation columns: {', '.join(h5ad_stats['obs_columns'])}\n")
    f.write(f"Variable columns: {', '.join(h5ad_stats['var_columns'])}\n")
    f.write(f"Unstructured data keys: {', '.join(h5ad_stats['uns_keys'])}\n")
    f.write(f"Multidimensional obs keys: {', '.join(h5ad_stats['obsm_keys'])}\n")
    f.write(f"Pairwise obs keys: {', '.join(h5ad_stats['obsp_keys'])}\n")
    f.write(f"Has raw data: {'Yes' if h5ad_stats['has_raw_data'] else 'No'}\n")
    
    if 'clustering_methods' in h5ad_stats:
        f.write(f"Clustering methods: {', '.join(h5ad_stats['clustering_methods'])}\n")
    if 'annotation_methods' in h5ad_stats:
        f.write(f"Annotation methods: {', '.join(h5ad_stats['annotation_methods'])}\n")
    f.write("\n")
    
    f.write("DETAILED VALIDATION CHECKS:\n")
    for i, check in enumerate(comparison_results['validation_checks'], 1):
        f.write(f"{i}. {check['check_name'].replace('_', ' ').title()}:\n")
        
        if 'match' in check:
            f.write(f"   Status: {'PASS' if check['match'] else 'FAIL'}\n")
            f.write(f"   GEF value: {check['gef_value']}\n")
            f.write(f"   H5AD value: {check['h5ad_value']}\n")
            if 'difference' in check:
                f.write(f"   Difference: {check['difference']}\n")
        
        if 'note' in check:
            f.write(f"   Note: {check['note']}\n")
        f.write("\n")
    
    f.write("EXPRESSION DATA COMPARISON:\n")
    f.write(f"GEF sample mean: {gef_stats.get('sample_mean_expression', 'N/A'):.4f}\n")
    f.write(f"H5AD sample mean: {h5ad_stats.get('sample_mean_expression', 'N/A'):.4f}\n")
    f.write(f"GEF sample max: {gef_stats.get('sample_max_expression', 'N/A'):.4f}\n")
    f.write(f"H5AD sample max: {h5ad_stats.get('sample_max_expression', 'N/A'):.4f}\n")
    f.write(f"GEF nonzero fraction: {gef_stats.get('sample_nonzero_fraction', 'N/A'):.4f}\n")
    f.write(f"H5AD nonzero fraction: {h5ad_stats.get('sample_nonzero_fraction', 'N/A'):.4f}\n")
    f.write("Note: Expression values differ due to normalization and log transformation\n\n")
    
log_step(f"Comprehensive validation report saved to {validation_report_file}")

# Print summary to console
log_step("="*80)
log_step("VALIDATION SUMMARY")
log_step("="*80)
log_step(f"Overall Status: {comparison_results['overall_status']}")
log_step(f"Original GEF: {gef_stats['n_cells']} cells × {gef_stats['n_genes']} genes ({gef_stats['file_info']['mb']:.1f} MB)")
log_step(f"H5AD Export: {h5ad_stats['n_obs']} cells × {h5ad_stats['n_vars']} genes ({h5ad_stats['file_info']['mb']:.1f} MB)")

critical_issues = []
for check in comparison_results['validation_checks']:
    if check.get('match') == False and 'cell_count' in check['check_name']:
        critical_issues.append(f"Cell count mismatch: {check['gef_value']} → {check['h5ad_value']}")
    elif check.get('x_range_match') == False or check.get('y_range_match') == False:
        critical_issues.append("Spatial coordinates do not match")

if critical_issues:
    log_step("CRITICAL ISSUES FOUND:")
    for issue in critical_issues:
        log_step(f"{issue}")
else:
    log_step("All critical validation checks passed")

log_step(f"Gene count difference: {abs(gef_stats['n_genes'] - h5ad_stats['n_vars'])} (expected due to filtering)")
log_step(f"File size reduction: {((gef_stats['file_info']['mb'] - h5ad_stats['file_info']['mb']) / gef_stats['file_info']['mb'] * 100):.1f}%")

# System info
log_memory_usage("FINAL")
log_step(f"Validation completed in {total_duration:.2f} seconds")
log_step(f"Results saved to: {validation_output_dir}/")

EOF

echo ""
echo "==========================================="
echo "Data validation analysis script created successfully!"
echo "==========================================="
echo ""

# Execute the validation analysis
echo "==========================================="
echo "Starting comprehensive data validation"
echo "==========================================="
echo ""

$ST_PYTHON bin/data_validation_analysis.py

if [ $? -eq 0 ]; then
    echo ""
    echo "==========================================="
    echo "DATA VALIDATION COMPLETED SUCCESSFULLY!"
    echo "==========================================="
    echo ""
    echo "Validation results saved to: results_validation/"
    echo "==========================================="
    echo ""
else
    echo ""
    echo "==========================================="
    echo "DATA VALIDATION FAILED!"
    echo "==========================================="
    echo "Check error logs for details"
    echo "==========================================="
    exit 1
fi

echo ""
echo "==========================================="
echo "Final system status:"
free -h
echo ""
echo "DATA VALIDATION FINISHED:"
date
echo "==========================================="
