#!/bin/bash
#$ -V
#$ -cwd
#$ -pe smp 16
#$ -l h_vmem=12G
#$ -N STEREOPY_ANALYSIS
#$ -o stereopy_ultimate_$JOB_ID.out
#$ -e stereopy_ultimate_$JOB_ID.err

echo "==========================================="
echo "STEREOPY ANALYSIS"
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
    echo "   qsub -v ST_PYTHON=/home/user/.conda/envs/st/bin/python,MIN_COUNTS=50,MIN_GENES=5,PCT_COUNTS_MT=100,N_PCS=30 bin/2_DOC_ANALYSIS.sh" >&2
    
    exit 1
fi

# Load environment with explicit paths
echo "Loading miniconda3 module..."
module load miniconda3

ST_BIN_DIR=$(dirname "$ST_PYTHON")
export PATH="$ST_BIN_DIR:$PATH"

echo "Using Python: $ST_PYTHON"
$ST_PYTHON --version

echo "Verifying all dependencies..."
$ST_PYTHON -c "
import sys
print(f'Python executable: {sys.executable}')
import stereo as st
print(f'Stereopy version: {st.__version__}')
import pandas as pd
print(f'Pandas version: {pd.__version__}')
import numpy as np
print(f'NumPy version: {np.__version__}')
import scipy
print(f'SciPy version: {scipy.__version__}')
import matplotlib
print(f'Matplotlib version: {matplotlib.__version__}')
import seaborn as sns
print(f'Seaborn version: {sns.__version__}')
"
echo ""

# Set thresholds
MIN_COUNTS=${MIN_COUNTS:-20}
MIN_GENES=${MIN_GENES:-3}
PCT_COUNTS_MT=${PCT_COUNTS_MT:-5}
N_PCS=${N_PCS:-10}
echo "Filtering Parameters Used:"
echo "  Minimum counts: $MIN_COUNTS"
echo "  Minimum genes per cell: $MIN_GENES"
echo "  Mitochondrial percentage: $PCT_COUNTS_MT"
echo "Number of Principal Components used:"
echo "  $N_PCS"
echo "  This step can be inproved after first run. Ceck the Elbow Plot (plots/qc/pca_elbow_enhanced.png) and insert the value of the elbow as N_PCS"
echo "You can alter the parameters inline:"
echo "  qsub -v ST_PYTHON="/home/user/.conda/envs/st/bin/python",MIN_COUNTS=50,MIN_GENES=5,PCT_COUNTS_MT=100,N_PCS=30 bin/2_DOC_ANALYSIS.sh"
echo ""

# Generate analysis script
echo "Creating analysis script for Stereopy"
cat > bin/stereopy_ultimate_analysis.py << EOF
#!/usr/bin/env python3
"""
Stereopy Spatial Transcriptomics Analysis
"""
# Import dependencies
import stereo as st
import warnings
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import mannwhitneyu, ttest_ind
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import psutil
import gc
from datetime import datetime
import sys
from collections import defaultdict
import json
import glob

# Set matplotlib backend for cluster
import matplotlib
matplotlib.use('Agg')

warnings.filterwarnings('ignore')

# Function for monitoring memory usage
def log_memory_usage(step_name=""):
    try:
        memory = psutil.virtual_memory()
        swap = psutil.swap_memory()
        print(f"[MEMORY {step_name}] RAM: {memory.percent:.1f}% ({memory.used/1e9:.1f}GB/{memory.total/1e9:.1f}GB) | Swap: {swap.percent:.1f}%")
    except:
        print(f"[MEMORY {step_name}] Unable to get memory info")

# Function for logging stamps
def log_step(message):
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()

# Function to create checkpoints
def save_progress_checkpoint(data, output_dir, checkpoint_name):
    checkpoint_dir = os.path.join(output_dir, 'checkpoints')
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    checkpoint_file = os.path.join(checkpoint_dir, f'{checkpoint_name}.json')
    checkpoint_info = {
        'timestamp': datetime.now().isoformat(),
        'n_cells': data.n_cells if hasattr(data, 'n_cells') else 0,
        'n_genes': data.n_genes if hasattr(data, 'n_genes') else 0,
        'checkpoint': checkpoint_name
    }
    
    with open(checkpoint_file, 'w') as f:
        json.dump(checkpoint_info, f, indent=2)

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

# Script Initialization
start_time = datetime.now()
log_step("Starting Stereopy analysis")
log_memory_usage("START")

# Filepaths
data_path = get_single_gef_file()
output_dir = 'RESULTS/results_ultimate'

# Directory structure
directories = [
    output_dir,
    os.path.join(output_dir, 'plots'),
    os.path.join(output_dir, 'plots', 'qc'),
    os.path.join(output_dir, 'plots', 'clustering'), 
    os.path.join(output_dir, 'plots', 'marker_genes'),
    os.path.join(output_dir, 'marker_genes'),
    os.path.join(output_dir, 'marker_genes', 'complete_results'),
    os.path.join(output_dir, 'marker_genes', 'filtered_results'),
    os.path.join(output_dir, 'statistical_analysis'),
    os.path.join(output_dir, 'statistical_analysis', 'advanced'),
    os.path.join(output_dir, 'cluster_data'),
    os.path.join(output_dir, 'cluster_data', 'detailed'),
    os.path.join(output_dir, 'logs'),
    os.path.join(output_dir, 'checkpoints'),
    os.path.join(output_dir, 'exports')
]

for directory in directories:
    os.makedirs(directory, exist_ok=True)

log_step(f"Created {len(directories)} output folders")

# Enhanced analysis log
analysis_log_path = os.path.join(output_dir, 'logs', 'ultimate_analysis_log.txt')
with open(analysis_log_path, 'w') as f:
    f.write("="*100 + "\n")
    f.write("STEREOPY SPATIAL TRANSCRIPTOMICS ANALYSIS LOG\n")
    f.write("="*100 + "\n")
    f.write(f"Analysis start: {start_time}\n")
    f.write(f"Data file: {data_path}\n")
    f.write(f"Python: {sys.version}\n")
    f.write(f"Stereopy: {st.__version__}\n")
    f.write(f"System RAM: {psutil.virtual_memory().total/1e9:.1f} GB\n")
    f.write(f"CPU cores: {psutil.cpu_count()}\n")
    f.write("="*100 + "\n\n")

# Verify data file and load
if not os.path.exists(data_path):
    log_step(f"ERROR: Data file not found at {data_path}")
    sys.exit(1)

log_step("Loading data")
try:
    st.io.read_gef_info(data_path)
    data = st.io.read_gef(file_path=data_path, bin_size=100)
    log_step(f"Data loaded successfully: {data}")
    save_progress_checkpoint(data, output_dir, 'data_loaded')
    log_memory_usage("DATA_LOADED")
except Exception as e:
    log_step(f"ERROR loading data: {e}")
    sys.exit(1)

# Filtering parameters from bash
MIN_COUNTS = $MIN_COUNTS
MIN_GENES = $MIN_GENES
PCT_COUNTS_MT = $PCT_COUNTS_MT

# General statistics
initial_stats_file = os.path.join(output_dir, 'logs', 'comprehensive_data_stats.txt')
with open(initial_stats_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("COMPREHENSIVE DATA STATISTICS\n")
    f.write("="*80 + "\n\n")
    f.write("RAW DATA CHARACTERISTICS:\n")
    f.write(str(data) + "\n\n")
    
    if hasattr(data, 'cells') and hasattr(data.cells, 'shape'):
        f.write(f"Cells matrix shape: {data.cells.shape}\n")
    if hasattr(data, 'genes') and hasattr(data.genes, 'shape'):
        f.write(f"Genes matrix shape: {data.genes.shape}\n")
    
    f.write(f"\nFILTERING CRITERIA:\n")
    f.write(f"- Minimum counts per cell: {MIN_COUNTS}\n")
    f.write(f"- Minimum genes per cell: {MIN_GENES}\n")
    f.write(f"- Maximum mitochondrial percentage: {PCT_COUNTS_MT}%\n\n")

log_step("Starting enhanced preprocessing pipeline...")

# Cell filtering with detailed statistics
log_step("Performing cell filtering with detailed tracking...")
cells_before = data.n_cells
genes_before = data.n_genes

try:
    data.tl.filter_cells(
        min_counts=MIN_COUNTS,
        min_genes=MIN_GENES,
        pct_counts_mt=PCT_COUNTS_MT,
        inplace=True
    )
    
    cells_after = data.n_cells
    genes_after = data.n_genes
    
    log_step(f"Filtering complete: {cells_before} -> {cells_after} cells ({cells_before-cells_after} removed)")
    log_step(f"Genes retained: {genes_after} (from {genes_before})")
    
    # Update stats
    with open(initial_stats_file, 'a') as f:
        f.write("AFTER FILTERING:\n")
        f.write(str(data) + "\n")
        f.write(f"Cells removed: {cells_before - cells_after} ({(cells_before-cells_after)/cells_before*100:.1f}%)\n")
        f.write(f"Genes retained: {genes_after}\n\n")
    
    save_progress_checkpoint(data, output_dir, 'cells_filtered')
    log_memory_usage("FILTERED")
    
except Exception as e:
    log_step(f"ERROR in cell filtering: {e}")
    sys.exit(1)

# QC plots before normalization
log_step("Generating QC plots before normalization")
qc_dir = os.path.join(output_dir, 'plots', 'qc')
try:
    data.plt.violin(out_path=os.path.join(qc_dir, 'violin_pre_normalization.png'))
    data.plt.spatial_scatter(out_path=os.path.join(qc_dir, 'spatial_scatter_pre_norm.png'))
    data.plt.genes_count(out_path=os.path.join(qc_dir, 'genes_count_pre_norm.png'))
    log_step("Pre-normalization QC plots saved")
except Exception as e:
    log_step(f"Warning: Error in pre-normalization plots: {e}")

# Normalization
log_step("Performing normalization and log transformation")
try:
    data.tl.normalize_total(target_sum=10000)
    data.tl.log1p()
    save_progress_checkpoint(data, output_dir, 'normalized')
    log_memory_usage("NORMALIZED")
except Exception as e:
    log_step(f"ERROR in normalization: {e}")
    sys.exit(1)

# Enhanced QC plots after normalization
log_step("Generating comprehensive QC plots after normalization")
try:
    data.plt.violin(out_path=os.path.join(qc_dir, 'violin_post_normalization.png'))
    data.plt.spatial_scatter(out_path=os.path.join(qc_dir, 'spatial_scatter_post_norm.png'))
    data.plt.genes_count(out_path=os.path.join(qc_dir, 'genes_count_post_norm.png'))
    log_step("Post-normalization QC plots saved")
except Exception as e:
    log_step(f"Warning: Error in post-normalization plots: {e}")

# Create raw checkpoint
log_step("Creating enhanced raw data checkpoint")
try:
    data.tl.raw_checkpoint()
    save_progress_checkpoint(data, output_dir, 'raw_checkpoint')
    log_step("Raw checkpoint created")
except Exception as e:
    log_step(f"ERROR creating raw checkpoint: {e}")
    sys.exit(1)

# Highly variable genes identification
log_step("Identifying highly variable genes")
try:
    data.tl.highly_variable_genes(
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_top_genes=2000,
        res_key='highly_variable_genes'
    )

    data.plt.highly_variable_genes(
        res_key='highly_variable_genes',
        out_path=os.path.join(qc_dir, 'highly_variable_genes_enhanced.png')
    )
    
    # Save HVG list
    if 'highly_variable_genes' in data.tl.result:
        hvg_file = os.path.join(output_dir, 'exports', 'highly_variable_genes.csv')
        hvg_data = data.tl.result['highly_variable_genes']
        if isinstance(hvg_data, pd.DataFrame):
            hvg_data.to_csv(hvg_file, index=False)
            log_step(f"Highly variable genes list saved ({len(hvg_data)} genes)")
    
    save_progress_checkpoint(data, output_dir, 'hvg_identified')
    log_memory_usage("HVG_IDENTIFIED")
    
except Exception as e:
    log_step(f"ERROR in highly variable genes: {e}")
    sys.exit(1)

# Data scaling
log_step("Performing data scaling")
try:
    data.tl.scale(max_value=10, zero_center=False)
    save_progress_checkpoint(data, output_dir, 'scaled')
except Exception as e:
    log_step(f"ERROR in scaling: {e}")
    sys.exit(1)

# PCA
N_PCS = $N_PCS
log_step("Performing PCA analysis...")
try:
    data.tl.pca(
        use_highly_genes=True,
        n_pcs=N_PCS,
        res_key='pca'
    )

    data.plt.elbow(
        pca_res_key='pca',
        out_path=os.path.join(qc_dir, 'pca_elbow_enhanced.png')
    )
    
    # Save PCA results
    if 'pca_variance_ratio' in data.tl.result:
        pca_file = os.path.join(output_dir, 'exports', 'pca_variance_ratios.csv')
        pca_data = pd.DataFrame({
            'PC': range(1, len(data.tl.result['pca_variance_ratio']) + 1),
            'variance_ratio': data.tl.result['pca_variance_ratio'],
            'cumulative_variance': np.cumsum(data.tl.result['pca_variance_ratio'])
        })
        pca_data.to_csv(pca_file, index=False)
        log_step(f"PCA variance ratios saved ({len(pca_data)} components)")
    
    save_progress_checkpoint(data, output_dir, 'pca_complete')
    log_memory_usage("PCA_COMPLETE")
    
except Exception as e:
    log_step(f"ERROR in PCA: {e}")
    sys.exit(1)

# Neighborhood graph computation
log_step("Computing neighborhood graph...")
try:
    data.tl.neighbors(
        pca_res_key='pca',
        n_pcs=N_PCS,
        res_key='neighbors'
    )
    save_progress_checkpoint(data, output_dir, 'neighbors_computed')
    log_memory_usage("NEIGHBORS_COMPUTED")
except Exception as e:
    log_step(f"ERROR in neighbors computation: {e}")
    sys.exit(1)

# UMAP computation
log_step("Computing enhanced UMAP embedding...")
try:
    data.tl.umap(
        pca_res_key='pca',
        neighbors_res_key='neighbors',
        res_key='umap'
    )
    
    # Save UMAP coordinates
    if hasattr(data, 'tl') and hasattr(data.tl, 'result') and 'umap' in data.tl.result:
        umap_coords = data.tl.result['umap']
        if isinstance(umap_coords, np.ndarray):
            umap_df = pd.DataFrame(umap_coords, columns=['UMAP1', 'UMAP2'])
        elif isinstance(umap_coords, pd.DataFrame):
            umap_df = umap_coords
        else:
            umap_df = None
        if umap_df is not None:
            if hasattr(data, 'cells') and hasattr(data.cells, 'cell_name'):
                umap_df['cell_name'] = data.cells.cell_name
            
            umap_file = os.path.join(output_dir, 'exports', 'umap_coordinates.csv')
            umap_df.to_csv(umap_file, index=False)
            log_step("UMAP coordinates saved")
    
    save_progress_checkpoint(data, output_dir, 'umap_complete')
    log_memory_usage("UMAP_COMPLETE")
    
except Exception as e:
    log_step(f"ERROR in UMAP: {e}")
    sys.exit(1)

# Clustering
log_step("Performing enhanced Louvain clustering...")
try:
    data.tl.louvain(
        neighbors_res_key='neighbors',
        res_key='louvain'
    )
    
    clustering_plots_dir = os.path.join(output_dir, 'plots', 'clustering')
    individual_clusters_dir = os.path.join(clustering_plots_dir, 'individual_clusters')
    os.makedirs(individual_clusters_dir, exist_ok=True)
    
    # Generate cluster plots
    log_step("Generating enhanced cluster visualization plots")
    data.plt.cluster_scatter(
        res_key='louvain',
        out_path=os.path.join(clustering_plots_dir, 'louvain_spatial_clusters_enhanced.png')
    )
    data.plt.umap(
        res_key='umap',
        cluster_key='louvain',
        out_path=os.path.join(clustering_plots_dir, 'louvain_umap_clusters_enhanced.png')
    )
    log_step("Generating individual cluster plots")
    unique_clusters = data.cells['louvain'].unique()
    n_clusters = len(unique_clusters)
    log_step(f"Found {n_clusters} clusters. Generating individual plots...")
    
    for i, cluster in enumerate(sorted(unique_clusters)):
        log_step(f"Processing cluster {cluster} ({i+1}/{n_clusters})...")
        
        try:
            # Create a copy of cluster annotations
            cluster_highlight = data.cells['louvain'].copy()
            
            # Handle Categorical data properly
            if hasattr(cluster_highlight, 'cat'):
                if 'Other' not in cluster_highlight.cat.categories:
                    cluster_highlight = cluster_highlight.cat.add_categories(['Other'])
            
            # Set all other clusters to 'Other'
            cluster_highlight[cluster_highlight != cluster] = 'Other'
            
            # Create the plot
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Get spatial coordinates
            if hasattr(data, 'position'):
                spatial_coords = data.position
                x_coords = spatial_coords[:, 0]
                y_coords = spatial_coords[:, 1]
            elif 'x' in data.cells.columns and 'y' in data.cells.columns:
                x_coords = data.cells['x'].values
                y_coords = data.cells['y'].values
            else:
                if hasattr(data, 'obsm') and 'spatial' in data.obsm:
                    spatial_coords = data.obsm['spatial']
                    x_coords = spatial_coords[:, 0]
                    y_coords = spatial_coords[:, 1]
                else:
                    log_step(f"ERROR: Could not find spatial coordinates for cluster {cluster}")
                    continue
            
            # Plot points
            other_mask = cluster_highlight == 'Other'
            if other_mask.any():
                ax.scatter(x_coords[other_mask], y_coords[other_mask], 
                          c='lightgray', s=1, alpha=0.5, rasterized=True)
            cluster_mask = cluster_highlight == cluster
            if cluster_mask.any():
                colors = plt.cm.tab20(int(cluster) % 20)
                ax.scatter(x_coords[cluster_mask], y_coords[cluster_mask], 
                          c=[colors], s=2, alpha=0.8, rasterized=True)
            
            # Formatting
            ax.set_xlabel('Spatial X (μm)', fontsize=12)
            ax.set_ylabel('Spatial Y (μm)', fontsize=12)
            ax.set_title(f'Cluster {cluster} Spatial Distribution', fontsize=14, fontweight='bold')
            
            # Add legend
            legend_elements = [
                mpatches.Patch(color='lightgray', label='Other clusters'),
                mpatches.Patch(color=colors, label=f'Cluster {cluster}')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
            
            # Add scale bar
            scale_bar_length = 2000  # 2mm in micrometers
            x_range = x_coords.max() - x_coords.min()
            y_range = y_coords.max() - y_coords.min()
            
            # Position scale bar at bottom left
            scale_x = x_coords.min() + 0.05 * x_range
            scale_y = y_coords.min() + 0.05 * y_range
            
            ax.plot([scale_x, scale_x + scale_bar_length], [scale_y, scale_y], 
                   'k-', linewidth=3)
            ax.text(scale_x + scale_bar_length/2, scale_y - 0.02 * y_range, 
                   '2.0mm', ha='center', va='top', fontsize=10, fontweight='bold')
            
            # Add cluster statistics as text
            n_cells_cluster = cluster_mask.sum()
            total_cells = len(cluster_highlight)
            percentage = (n_cells_cluster / total_cells) * 100
            
            stats_text = f'Cells in cluster: {n_cells_cluster:,}\nPercentage: {percentage:.1f}%'
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # Set equal aspect ratio and adjust layout
            ax.set_aspect('equal')
            plt.tight_layout()
            
            # Save the plot
            output_file = os.path.join(individual_clusters_dir, f'cluster_{cluster}_spatial.png')
            plt.savefig(output_file, dpi=300, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            plt.close()
            
            log_step(f"Saved individual plot for cluster {cluster} to {output_file}")
            
        except Exception as e:
            log_step(f"ERROR generating plot for cluster {cluster}: {e}")
            continue
    
    # Generate a summary HTML file to view all individual clusters
    log_step("Generating HTML summary for individual cluster plots...")
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Individual Cluster Analysis</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            .cluster-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }}
            .cluster-item {{ border: 1px solid #ddd; padding: 10px; border-radius: 5px; }}
            .cluster-item img {{ width: 100%; height: auto; }}
            .cluster-title {{ font-weight: bold; margin-bottom: 10px; }}
            h1 {{ color: #333; text-align: center; }}
            .summary {{ background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
        </style>
    </head>
    <body>
        <h1>Individual Cluster Spatial Distribution</h1>
        <div class="summary">
            <h3>Analysis Summary</h3>
            <p><strong>Total clusters found:</strong> {n_clusters}</p>
            <p><strong>Analysis type:</strong> Louvain clustering with spatial visualization</p>
            <p><strong>Visualization method:</strong> Individual cluster highlighting (target cluster in color, others in gray)</p>
        </div>
        <div class="cluster-grid">
    """
    
    for cluster in sorted(unique_clusters):
        html_content += f"""
            <div class="cluster-item">
                <div class="cluster-title">Cluster {cluster}</div>
                <img src="individual_clusters/cluster_{cluster}_spatial.png" alt="Cluster {cluster}">
            </div>
        """
    
    html_content += """
        </div>
    </body>
    </html>
    """
    
    html_file = os.path.join(clustering_plots_dir, 'individual_clusters_summary.html')
    with open(html_file, 'w') as f:
        f.write(html_content)
    
    log_step(f"Generated HTML summary: {html_file}")
    log_step(f"Individual cluster plots saved to: {individual_clusters_dir}")
    
    save_progress_checkpoint(data, output_dir, 'clustering_complete')
    log_memory_usage("CLUSTERING_COMPLETE")
    
except Exception as e:
    log_step(f"ERROR in clustering: {e}")
    import traceback
    log_step(f"Full traceback: {traceback.format_exc()}")
    sys.exit(1)

# MARKER GENE ANALYSIS
## This step process all genes and all clusters (requires HPC)
log_step("="*80)
log_step("MARKER GENE ANALYSIS")
log_step("="*80)

marker_start_time = datetime.now()
log_step("Finding marker genes using t-test method (COMPLETE DATASET)...")

try:
    data.tl.find_marker_genes(
        cluster_res_key='louvain',
        method='t_test',
        use_highly_genes=False,  # Use ALL genes
        use_raw=True,            # Use raw count data
        res_key='marker_genes'
    )

    marker_end_time = datetime.now()
    marker_duration = (marker_end_time - marker_start_time).total_seconds()
    log_step(f"Marker gene analysis completed in {marker_duration:.1f} seconds")
    log_memory_usage("MARKER_GENES_FOUND")
    
except Exception as e:
    log_step(f"ERROR in marker gene analysis: {e}")
    sys.exit(1)

# Marker genes processing
marker_genes_dir = os.path.join(output_dir, 'marker_genes')
complete_results_dir = os.path.join(marker_genes_dir, 'complete_results')

if 'marker_genes' in data.tl.result:
    try:
        marker_result = data.tl.result['marker_genes']
        log_step(f"Processing COMPLETE marker genes results with {len(marker_result)} components")
        
        # Get all cluster comparison keys
        cluster_keys = [k for k in marker_result.keys() if k.endswith('.vs.rest')]
        log_step(f"Processing ALL {len(cluster_keys)} cluster comparisons - NO LIMITATIONS")
        
        # Process ALL clusters with ALL genes
        all_marker_genes = []
        cluster_summaries = []
        processing_stats = []
        
        log_step("Processing all clusters with complete gene sets")
        
        for i, cluster_key in enumerate(cluster_keys):
            cluster_num = cluster_key.split('.')[0]
            cluster_start_time = datetime.now()
            
            try:
                cluster_data = marker_result[cluster_key]
                
                if isinstance(cluster_data, pd.DataFrame) and len(cluster_data) > 0:
                    # Add cluster identifier
                    cluster_data = cluster_data.copy()
                    cluster_data['cluster'] = cluster_num
                    
                    # Save COMPLETE individual results (ALL genes)
                    cluster_file = os.path.join(complete_results_dir, f'cluster_{cluster_num}_complete_markers.csv')
                    cluster_data.to_csv(cluster_file, index=False)
                    
                    # Add ALL genes to combined analysis (no filtering)
                    all_marker_genes.append(cluster_data)
                    
                    # Statistics calculation
                    if 'pvalues_adj' in cluster_data.columns and 'log2fc' in cluster_data.columns:
                        sig_001 = cluster_data[cluster_data['pvalues_adj'] < 0.001]
                        sig_01 = cluster_data[cluster_data['pvalues_adj'] < 0.01]
                        sig_05 = cluster_data[cluster_data['pvalues_adj'] < 0.05]
                        fc_05 = cluster_data[abs(cluster_data['log2fc']) > 0.5]
                        fc_1 = cluster_data[abs(cluster_data['log2fc']) > 1]
                        fc_2 = cluster_data[abs(cluster_data['log2fc']) > 2]
                        
                        # Combined thresholds
                        stringent = cluster_data[
                            (cluster_data['pvalues_adj'] < 0.001) & 
                            (abs(cluster_data['log2fc']) > 1)
                        ]
                        moderate = cluster_data[
                            (cluster_data['pvalues_adj'] < 0.01) & 
                            (abs(cluster_data['log2fc']) > 0.5)
                        ]
                        lenient = cluster_data[
                            (cluster_data['pvalues_adj'] < 0.05) & 
                            (abs(cluster_data['log2fc']) > 0.25)
                        ]
                        
                        summary = {
                            'cluster': cluster_num,
                            'total_genes_tested': len(cluster_data),
                            'sig_p001': len(sig_001),
                            'sig_p01': len(sig_01),
                            'sig_p05': len(sig_05),
                            'fc_05': len(fc_05),
                            'fc_1': len(fc_1),
                            'fc_2': len(fc_2),
                            'stringent_markers': len(stringent),
                            'moderate_markers': len(moderate),
                            'lenient_markers': len(lenient),
                            'mean_log2fc': cluster_data['log2fc'].mean(),
                            'median_log2fc': cluster_data['log2fc'].median(),
                            'max_log2fc': cluster_data['log2fc'].max(),
                            'min_log2fc': cluster_data['log2fc'].min(),
                            'min_pval_adj': cluster_data['pvalues_adj'].min(),
                            'mean_pval_adj': cluster_data['pvalues_adj'].mean()
                        }
                        cluster_summaries.append(summary)
                    
                    # Processing time tracking
                    cluster_end_time = datetime.now()
                    processing_time = (cluster_end_time - cluster_start_time).total_seconds()
                    processing_stats.append({
                        'cluster': cluster_num,
                        'processing_time_seconds': processing_time,
                        'genes_processed': len(cluster_data)
                    })
                    
                    if (i + 1) % 10 == 0:
                        log_step(f"Processed {i + 1}/{len(cluster_keys)} clusters")
                        log_memory_usage(f"PROCESSED_{i+1}")
                
            except Exception as e:
                log_step(f"Error processing cluster {cluster_key}: {e}")
                continue
        
        # Combine ALL marker genes (complete dataset)
        if all_marker_genes:
            log_step("Combining COMPLETE marker genes dataset...")
            log_step("WARNING: This will create a very large file with ALL genes for ALL clusters")
            combined_complete_markers = pd.concat(all_marker_genes, ignore_index=True)
            complete_markers_file = os.path.join(marker_genes_dir, 'COMPLETE_all_marker_genes_no_limits.csv')
            combined_complete_markers.to_csv(complete_markers_file, index=False)
            
            log_step(f"COMPLETE marker genes dataset saved: {len(combined_complete_markers)} total entries")
            log_step(f"File size: {os.path.getsize(complete_markers_file) / 1e9:.2f} GB")
            
            # Create filtered versions for practical use
            log_step("Creating filtered versions for practical analysis...")
            
            # Stringent markers
            if 'pvalues_adj' in combined_complete_markers.columns and 'log2fc' in combined_complete_markers.columns:
                stringent_markers = combined_complete_markers[
                    (combined_complete_markers['pvalues_adj'] < 0.001) & 
                    (abs(combined_complete_markers['log2fc']) > 1)
                ]
                stringent_file = os.path.join(marker_genes_dir, 'stringent_markers.csv')
                stringent_markers.to_csv(stringent_file, index=False)
                log_step(f"Stringent markers saved: {len(stringent_markers)} genes")
                
                # Moderate markers
                moderate_markers = combined_complete_markers[
                    (combined_complete_markers['pvalues_adj'] < 0.01) & 
                    (abs(combined_complete_markers['log2fc']) > 0.5)
                ]
                moderate_file = os.path.join(marker_genes_dir, 'moderate_markers.csv')
                moderate_markers.to_csv(moderate_file, index=False)
                log_step(f"Moderate markers saved: {len(moderate_markers)} genes")
                
                # Top markers per cluster (for visualization)
                top_markers_per_cluster = []
                if 'scores' in combined_complete_markers.columns:
                    for cluster in combined_complete_markers['cluster'].unique():
                        cluster_markers = combined_complete_markers[
                            combined_complete_markers['cluster'] == cluster
                        ].nlargest(50, 'scores')
                        top_markers_per_cluster.append(cluster_markers)
                    
                    top_combined = pd.concat(top_markers_per_cluster, ignore_index=True)
                    top_file = os.path.join(marker_genes_dir, 'top50_per_cluster_for_viz.csv')
                    top_combined.to_csv(top_file, index=False)
                    log_step(f"Top 50 per cluster saved for visualization: {len(top_combined)} genes")
        
        # Save comprehensive cluster summaries
        if cluster_summaries:
            summary_df = pd.DataFrame(cluster_summaries)
            summary_file = os.path.join(marker_genes_dir, 'comprehensive_cluster_marker_summary.csv')
            summary_df.to_csv(summary_file, index=False)
            
            # Enhanced summary statistics
            detailed_summary_file = os.path.join(marker_genes_dir, 'detailed_marker_analysis_report.txt')
            with open(detailed_summary_file, 'w') as f:
                f.write("ULTIMATE MARKER GENE ANALYSIS REPORT\n")
                f.write("="*80 + "\n\n")
                f.write(f"Analysis completed: {datetime.now()}\n")
                f.write(f"Total clusters analyzed: {len(cluster_summaries)}\n")
                f.write(f"Processing time: {marker_duration:.1f} seconds\n\n")
                
                f.write("SIGNIFICANCE THRESHOLDS SUMMARY:\n")
                f.write(f"Total stringent markers (p<0.001, |FC|>1): {summary_df['stringent_markers'].sum()}\n")
                f.write(f"Total moderate markers (p<0.01, |FC|>0.5): {summary_df['moderate_markers'].sum()}\n")
                f.write(f"Total lenient markers (p<0.05, |FC|>0.25): {summary_df['lenient_markers'].sum()}\n\n")
                
                f.write("TOP 10 CLUSTERS BY STRINGENT MARKERS:\n")
                top_stringent = summary_df.nlargest(10, 'stringent_markers')[['cluster', 'stringent_markers']]
                f.write(str(top_stringent) + "\n\n")
                
                f.write("FOLD CHANGE STATISTICS:\n")
                f.write(f"Average max fold change: {summary_df['max_log2fc'].mean():.2f}\n")
                f.write(f"Average min fold change: {summary_df['min_log2fc'].mean():.2f}\n")
                f.write(f"Overall range: {summary_df['max_log2fc'].max():.2f} to {summary_df['min_log2fc'].min():.2f}\n\n")
                
                f.write("DETAILED CLUSTER STATISTICS:\n")
                f.write(str(summary_df) + "\n")
        
        # Save processing statistics
        if processing_stats:
            processing_df = pd.DataFrame(processing_stats)
            processing_file = os.path.join(output_dir, 'logs', 'cluster_processing_statistics.csv')
            processing_df.to_csv(processing_file, index=False)
            
            total_processing_time = processing_df['processing_time_seconds'].sum()
            avg_processing_time = processing_df['processing_time_seconds'].mean()
            log_step(f"Processing statistics: Total {total_processing_time:.1f}s, Average {avg_processing_time:.2f}s per cluster")
        
        log_memory_usage("ALL_MARKERS_PROCESSED")
        
    except Exception as e:
        log_step(f"ERROR processing marker genes: {e}")
        import traceback
        log_step(f"Traceback: {traceback.format_exc()}")

# ADVANCED STATISTICAL ANALYSIS

log_step("="*70)
log_step("ADVANCED STATISTICAL ANALYSIS")
log_step("="*70)

advanced_stats_dir = os.path.join(output_dir, 'statistical_analysis', 'advanced')

# Load the complete marker genes for advanced analysis
complete_markers_file = os.path.join(marker_genes_dir, 'COMPLETE_all_marker_genes_no_limits.csv')
if os.path.exists(complete_markers_file):
    try:
        log_step("Loading complete marker genes dataset for advanced analysis...")
        df_complete_markers = pd.read_csv(complete_markers_file)
        log_step(f"Loaded {len(df_complete_markers)} marker gene records for analysis")
        
        # Advanced statistical analyses
        log_step("Performing advanced statistical analyses...")
        
        # 1. Distribution analysis
        if 'log2fc' in df_complete_markers.columns and 'pvalues_adj' in df_complete_markers.columns:
            
            # Statistical tests for normality
            log_step("Testing data distributions...")
            
            # Test log2fc distribution
            from scipy.stats import shapiro, kstest
            lfc_sample = df_complete_markers['log2fc'].dropna().sample(min(5000, len(df_complete_markers)))
            shapiro_stat, shapiro_p = shapiro(lfc_sample)
            
            # Save statistical test results
            stats_results_file = os.path.join(advanced_stats_dir, 'distribution_tests.txt')
            with open(stats_results_file, 'w') as f:
                f.write("ADVANCED STATISTICAL DISTRIBUTION ANALYSIS\n")
                f.write("="*60 + "\n\n")
                f.write(f"Dataset size: {len(df_complete_markers)} records\n")
                f.write(f"Unique clusters: {df_complete_markers['cluster'].nunique()}\n")
                f.write(f"Unique genes: {df_complete_markers['genes'].nunique() if 'genes' in df_complete_markers.columns else 'Unknown'}\n\n")
                
                f.write("LOG2 FOLD CHANGE DISTRIBUTION:\n")
                f.write(f"Mean: {df_complete_markers['log2fc'].mean():.4f}\n")
                f.write(f"Median: {df_complete_markers['log2fc'].median():.4f}\n")
                f.write(f"Std: {df_complete_markers['log2fc'].std():.4f}\n")
                f.write(f"Skewness: {df_complete_markers['log2fc'].skew():.4f}\n")
                f.write(f"Kurtosis: {df_complete_markers['log2fc'].kurtosis():.4f}\n")
                f.write(f"Shapiro-Wilk test (sample): W={shapiro_stat:.4f}, p={shapiro_p:.2e}\n\n")
                
                f.write("ADJUSTED P-VALUES DISTRIBUTION:\n")
                f.write(f"Mean: {df_complete_markers['pvalues_adj'].mean():.6f}\n")
                f.write(f"Median: {df_complete_markers['pvalues_adj'].median():.6f}\n")
                f.write(f"Min: {df_complete_markers['pvalues_adj'].min():.2e}\n")
                f.write(f"Max: {df_complete_markers['pvalues_adj'].max():.6f}\n")
        
        # 2. Cluster comparison analysis
        if 'cluster' in df_complete_markers.columns:
            log_step("Performing inter-cluster comparison analysis...")
            
            cluster_comparison_file = os.path.join(advanced_stats_dir, 'cluster_comparisons.csv')
            cluster_stats = []
            
            for cluster in df_complete_markers['cluster'].unique():
                cluster_data = df_complete_markers[df_complete_markers['cluster'] == cluster]
                
                if len(cluster_data) > 0 and 'log2fc' in cluster_data.columns:
                    stats_dict = {
                        'cluster': cluster,
                        'n_genes': len(cluster_data),
                        'mean_lfc': cluster_data['log2fc'].mean(),
                        'median_lfc': cluster_data['log2fc'].median(),
                        'std_lfc': cluster_data['log2fc'].std(),
                        'min_lfc': cluster_data['log2fc'].min(),
                        'max_lfc': cluster_data['log2fc'].max(),
                        'q25_lfc': cluster_data['log2fc'].quantile(0.25),
                        'q75_lfc': cluster_data['log2fc'].quantile(0.75),
                    }
                    
                    if 'pvalues_adj' in cluster_data.columns:
                        stats_dict.update({
                            'mean_pval': cluster_data['pvalues_adj'].mean(),
                            'median_pval': cluster_data['pvalues_adj'].median(),
                            'min_pval': cluster_data['pvalues_adj'].min(),
                            'sig_001': (cluster_data['pvalues_adj'] < 0.001).sum(),
                            'sig_01': (cluster_data['pvalues_adj'] < 0.01).sum(),
                            'sig_05': (cluster_data['pvalues_adj'] < 0.05).sum()
                        })
                    
                    cluster_stats.append(stats_dict)
            
            if cluster_stats:
                cluster_stats_df = pd.DataFrame(cluster_stats)
                cluster_stats_df.to_csv(cluster_comparison_file, index=False)
                log_step(f"Cluster comparison statistics saved ({len(cluster_stats)} clusters)")
        
        # 3. Generate advanced visualizations
        log_step("Generating advanced statistical visualizations")
        
        # Create comprehensive plots
        if len(df_complete_markers) > 0:
            fig, axes = plt.subplots(2, 3, figsize=(20, 12))
            
            # Plot 1: Log2FC distribution
            if 'log2fc' in df_complete_markers.columns:
                axes[0,0].hist(df_complete_markers['log2fc'], bins=100, alpha=0.7, edgecolor='black')
                axes[0,0].axvline(0, color='red', linestyle='--', alpha=0.7)
                axes[0,0].set_xlabel('Log2 Fold Change')
                axes[0,0].set_ylabel('Frequency')
                axes[0,0].set_title('Distribution of Log2 Fold Changes')
            
            # Plot 2: P-value distribution
            if 'pvalues_adj' in df_complete_markers.columns:
                axes[0,1].hist(-np.log10(df_complete_markers['pvalues_adj'].replace(0, 1e-300)), 
                              bins=100, alpha=0.7, edgecolor='black')
                axes[0,1].set_xlabel('-log10(Adjusted P-value)')
                axes[0,1].set_ylabel('Frequency')
                axes[0,1].set_title('Distribution of -log10(Adj. P-values)')
            
            # Plot 3: Volcano plot (sample)
            if 'log2fc' in df_complete_markers.columns and 'pvalues_adj' in df_complete_markers.columns:
                sample_data = df_complete_markers.sample(min(10000, len(df_complete_markers)))
                axes[0,2].scatter(sample_data['log2fc'], 
                                -np.log10(sample_data['pvalues_adj'].replace(0, 1e-300)),
                                alpha=0.5, s=1)
                axes[0,2].axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.7)
                axes[0,2].axvline(0, color='red', linestyle='--', alpha=0.7)
                axes[0,2].set_xlabel('Log2 Fold Change')
                axes[0,2].set_ylabel('-log10(Adjusted P-value)')
                axes[0,2].set_title('Volcano Plot (Sample)')
            
            # Plot 4: Cluster-wise statistics
            if 'cluster' in df_complete_markers.columns and len(cluster_stats) > 0:
                cluster_means = [stat['mean_lfc'] for stat in cluster_stats]
                cluster_ids = [stat['cluster'] for stat in cluster_stats]
                axes[1,0].bar(range(len(cluster_means)), cluster_means)
                axes[1,0].set_xlabel('Cluster Index')
                axes[1,0].set_ylabel('Mean Log2FC')
                axes[1,0].set_title('Mean Log2FC by Cluster')
            
            # Plot 5: Significance by cluster
            if 'cluster' in df_complete_markers.columns and len(cluster_stats) > 0:
                sig_counts = [stat.get('sig_01', 0) for stat in cluster_stats]
                axes[1,1].bar(range(len(sig_counts)), sig_counts)
                axes[1,1].set_xlabel('Cluster Index')
                axes[1,1].set_ylabel('Significant Genes (p<0.01)')
                axes[1,1].set_title('Significant Genes by Cluster')
            
            # Plot 6: Log2FC vs Significance
            if 'log2fc' in df_complete_markers.columns and 'pvalues_adj' in df_complete_markers.columns:
                sample_data = df_complete_markers.sample(min(5000, len(df_complete_markers)))
                scatter = axes[1,2].scatter(abs(sample_data['log2fc']), 
                                          -np.log10(sample_data['pvalues_adj'].replace(0, 1e-300)),
                                          alpha=0.6, s=2, c=sample_data.get('scores', 1))
                axes[1,2].set_xlabel('|Log2 Fold Change|')
                axes[1,2].set_ylabel('-log10(Adjusted P-value)')
                axes[1,2].set_title('Effect Size vs Significance')
            
            plt.tight_layout()
            advanced_plots_file = os.path.join(advanced_stats_dir, 'comprehensive_statistical_analysis.png')
            plt.savefig(advanced_plots_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            log_step("Advanced statistical plots saved")
        
        log_memory_usage("ADVANCED_STATS_COMPLETE")
        
    except Exception as e:
        log_step(f"Warning: Error in advanced statistical analysis: {e}")

# VISUALIZATION

log_step("="*70)
log_step("VISUALIZATION SUITE")
log_step("="*70)

viz_marker_dir = os.path.join(output_dir, 'plots', 'marker_genes')

# Use the filtered dataset for visualizations (more manageable)
top_viz_file = os.path.join(marker_genes_dir, 'top50_per_cluster_for_viz.csv')
if os.path.exists(top_viz_file):
    try:
        log_step("Generating ultimate Stereopy visualizations...")
        
        # Standard Stereopy visualizations
        data.plt.marker_genes_text(
            res_key='marker_genes',
            markers_num=10,
            sort_key='scores',
            out_path=os.path.join(viz_marker_dir, 'ultimate_top15_marker_genes_text.png')
        )
        
        data.plt.marker_genes_scatter(
            res_key='marker_genes', 
            markers_num=10,
            out_path=os.path.join(viz_marker_dir, 'ultimate_top10_marker_genes_scatter.png')
        )
        
        log_step("Enhanced Stereopy visualizations completed")
        
    except Exception as e:
        log_step(f"Warning: Error in Stereopy visualizations: {e}")

# Generate volcano plots for ALL major clusters
try:
    log_step("Generating volcano plots for ALL major clusters...")
    volcano_dir = os.path.join(viz_marker_dir, 'volcano_plots_complete')
    os.makedirs(volcano_dir, exist_ok=True)
    
    if 'louvain' in data.cells:
        all_clusters = sorted(data.cells['louvain'].unique())
        log_step(f"Generating volcano plots for ALL {len(all_clusters)} clusters...")
        
        successful_plots = 0
        
        for i, cluster in enumerate(all_clusters):
            try:
                cluster_str = str(cluster)
                group_name = f'{cluster_str}.vs.rest'
                data.plt.marker_gene_volcano(
                    group_name=group_name,
                    res_key='marker_genes',
                    vlines=False,
                    out_path=os.path.join(volcano_dir, f'volcano_cluster_{cluster_str.zfill(2)}.png')
                )
                
                successful_plots += 1
                
                if (i + 1) % 20 == 0:
                    log_step(f"Generated {successful_plots}/{i + 1} volcano plots successfully")
                    
            except Exception as e:
                log_step(f"Warning: Error generating volcano plot for cluster {cluster}: {e}")
                try:
                    alt_group_name = f'{int(cluster)}.vs.rest'
                    data.plt.marker_gene_volcano(
                        group_name=alt_group_name,
                        res_key='marker_genes',
                        vlines=False,
                        out_path=os.path.join(volcano_dir, f'volcano_cluster_{str(cluster).zfill(2)}_alt.png')
                    )
                    successful_plots += 1
                    log_step(f"Alternative format worked for cluster {cluster}")
                except Exception as e2:
                    log_step(f"Both formats failed for cluster {cluster}: {e2}")
                continue
        
        log_step(f"Volcano plots completed: {successful_plots}/{len(all_clusters)} successful")
        
        # If no plots were generated, try a different approach
        if successful_plots == 0:
            log_step("No volcano plots generated. Trying alternative approach...")
            
            # Check what group names are actually available in marker_genes result
            if 'marker_genes' in data.tl.result:
                available_groups = [k for k in data.tl.result['marker_genes'].keys() if k.endswith('.vs.rest')]
                log_step(f"Available group names in results: {available_groups[:10]}...")
                
                # Try with first few available group names
                for j, group_name in enumerate(available_groups[:10]):
                    try:
                        cluster_id = group_name.replace('.vs.rest', '')
                        data.plt.marker_gene_volcano(
                            group_name=group_name,
                            res_key='marker_genes',
                            vlines=False,
                            out_path=os.path.join(volcano_dir, f'volcano_cluster_{cluster_id}_direct.png')
                        )
                        successful_plots += 1
                        log_step(f"Direct approach successful for {group_name}")
                    except Exception as e:
                        log_step(f"Direct approach failed for {group_name}: {e}")
                        continue
            
            log_step(f"Final volcano plots count: {successful_plots}")
    
except Exception as e:
    log_step(f"Warning: Error in volcano plot generation: {e}")

log_memory_usage("VISUALIZATIONS_COMPLETE")

# MARKER GENE FILTERING

log_step("="*70)
log_step("MARKER GENE FILTERING")
log_step("="*70)

filtered_results_dir = os.path.join(marker_genes_dir, 'filtered_results')

try:
    log_step("Applying advanced statistical filters to marker genes...")
    
    # Apply multiple filtering criteria
    filter_configs = [
        {'min_fold_change': 0.5, 'min_in_group_fraction': 0.2, 'max_out_group_fraction': 0.6, 'name': 'lenient'},
        {'min_fold_change': 1.0, 'min_in_group_fraction': 0.25, 'max_out_group_fraction': 0.5, 'name': 'moderate'},
        {'min_fold_change': 1.5, 'min_in_group_fraction': 0.3, 'max_out_group_fraction': 0.4, 'name': 'stringent'},
        {'min_fold_change': 2.0, 'min_in_group_fraction': 0.4, 'max_out_group_fraction': 0.3, 'name': 'very_stringent'}
    ]
    
    for config in filter_configs:
        try:
            filter_key = f"marker_genes_filtered_{config['name']}"
            
            data.tl.filter_marker_genes(
                marker_genes_res_key='marker_genes',
                min_fold_change=config['min_fold_change'],
                min_in_group_fraction=config['min_in_group_fraction'],
                max_out_group_fraction=config['max_out_group_fraction'],
                res_key=filter_key
            )
            
            log_step(f"Applied {config['name']} filtering (FC>{config['min_fold_change']})")
            
            # Export filtered results
            if filter_key in data.tl.result:
                filtered_result = data.tl.result[filter_key]
                
                if isinstance(filtered_result, dict):
                    cluster_keys = [k for k in filtered_result.keys() if k.endswith('.vs.rest')]
                    
                    all_filtered = []
                    for cluster_key in cluster_keys:
                        cluster_num = cluster_key.split('.')[0]
                        try:
                            cluster_data = filtered_result[cluster_key]
                            if isinstance(cluster_data, pd.DataFrame) and len(cluster_data) > 0:
                                cluster_data = cluster_data.copy()
                                cluster_data['cluster'] = cluster_num
                                all_filtered.append(cluster_data)
                        except:
                            continue
                    
                    if all_filtered:
                        combined_filtered = pd.concat(all_filtered, ignore_index=True)
                        
                        filter_file = os.path.join(filtered_results_dir, f'filtered_markers_{config["name"]}.csv')
                        combined_filtered.to_csv(filter_file, index=False)
                        
                        log_step(f"Saved {len(combined_filtered)} {config['name']} filtered markers")
            
        except Exception as e:
            log_step(f"Warning: Error in {config['name']} filtering: {e}")
            continue
    
    log_memory_usage("FILTERING_COMPLETE")
    
except Exception as e:
    log_step(f"Warning: Error in advanced filtering: {e}")

# FINAL REPORT

log_step("="*70)
log_step("GENERATING FINAL REPORT")
log_step("="*70)

end_time = datetime.now()
total_duration = (end_time - start_time).total_seconds()

# Generate ultimate final report
final_report_path = os.path.join(output_dir, 'ULTIMATE_ANALYSIS_REPORT.txt')
with open(final_report_path, 'w') as f:
    f.write("="*100 + "\n")
    f.write("STEREOPY SPATIAL TRANSCRIPTOMICS ANALYSIS REPORT\n")
    f.write("="*100 + "\n\n")
    f.write(f"Analysis completed: {end_time}\n")
    f.write(f"Total processing time: {total_duration/60:.1f} minutes ({total_duration:.1f} seconds)\n")
    f.write(f"Data file: {data_path}\n")
    f.write(f"Output directory: {output_dir}\n\n")
    
    # System utilization
    memory = psutil.virtual_memory()
    f.write("SYSTEM RESOURCES UTILIZED:\n")
    f.write(f"- Total RAM available: {memory.total/1e9:.1f} GB\n")
    f.write(f"- Peak memory usage: {memory.percent:.1f}% ({memory.used/1e9:.1f} GB)\n")
    f.write(f"- CPU cores utilized: {psutil.cpu_count()}\n")
    f.write(f"- Python version: {sys.version.split()[0]}\n")
    f.write(f"- Stereopy version: {st.__version__}\n\n")
    
    # Analysis results
    if 'louvain' in data.cells:
        cluster_counts = data.cells['louvain'].value_counts()
        f.write("LOUVAIN CLUSTERING ANALYSIS RESULTS:\n")
        f.write(f"- Total clusters identified: {len(cluster_counts)}\n")
        f.write(f"- Total cells analyzed: {cluster_counts.sum()}\n")
        f.write(f"- Largest cluster: {cluster_counts.max()} cells (Cluster {cluster_counts.idxmax()})\n")
        f.write(f"- Smallest cluster: {cluster_counts.min()} cells (Cluster {cluster_counts.idxmin()})\n")
        f.write(f"- Average cluster size: {cluster_counts.mean():.1f} cells\n")
        f.write(f"- Total genes analyzed: {data.n_genes}\n\n")
    
    # Key files summary
    f.write("KEY OUTPUT FILES:\n")
    f.write("- COMPLETE_all_marker_genes_no_limits.csv: Complete dataset\n")
    f.write("- stringent_markers.csv: High-confidence markers\n")
    f.write("- moderate_markers.csv: Standard analysis markers\n")
    f.write("- Individual cluster files: complete_results/cluster_XX_complete_markers.csv\n")
    f.write("- Advanced statistics: statistical_analysis/advanced/\n")
    f.write("- Complete visualizations: plots/marker_genes/volcano_plots_complete/\n")
    f.write("- Comprehensive logs: logs/\n\n")
    
    # Performance metrics
    f.write("PERFORMANCE METRICS:\n")
    f.write(f"- Marker gene analysis: {marker_duration:.1f} seconds\n")
    f.write(f"- Memory efficiency: Used only {memory.percent:.1f}% of available RAM\n")
    f.write(f"- Processing rate: ~{(cluster_counts.sum() * data.n_genes) / total_duration:.0f} gene-cell comparisons/second\n")
    f.write(f"- Data throughput: ~{total_duration/len(cluster_counts):.1f} seconds per cluster\n\n")

# Save processing summary
processing_summary_file = os.path.join(output_dir, 'logs', 'processing_summary.json')
processing_summary = {
    'analysis_type': 'ULTIMATE_NO_LIMITS',
    'start_time': start_time.isoformat(),
    'end_time': end_time.isoformat(),
    'total_duration_seconds': total_duration,
    'total_duration_minutes': total_duration/60,
    'data_file': data_path,
    'n_clusters': len(data.cells['louvain'].value_counts()) if 'louvain' in data.cells else 0,
    'n_cells': data.n_cells,
    'n_genes': data.n_genes,
    'marker_analysis_duration': marker_duration,
    'peak_memory_percent': psutil.virtual_memory().percent,
    'peak_memory_gb': psutil.virtual_memory().used/1e9,
    'system_ram_gb': psutil.virtual_memory().total/1e9,
    'cpu_cores': psutil.cpu_count(),
    'python_version': sys.version.split()[0],
    'stereopy_version': st.__version__
}

with open(processing_summary_file, 'w') as f:
    json.dump(processing_summary, f, indent=2)

log_step(f"Ultimate analysis report saved to {final_report_path}")
log_step(f"Processing summary saved to {processing_summary_file}")
log_memory_usage("FINAL")

# Final cleanup
gc.collect()

log_step("="*80)
log_step("ULTIMATE ANALYSIS COMPLETED SUCCESSFULLY!")
log_step(f"Total processing time: {total_duration/60:.1f} minutes")
log_step(f"Peak memory usage: {psutil.virtual_memory().percent:.1f}%")
log_step(f"Results location: {output_dir}")
log_step("="*80)

EOF

echo "==========================================="
echo "Script created successfully!"
echo "==========================================="

# Execute the analysis
echo "==========================================="
echo "Starting Stereopy analysis"
echo "==========================================="
echo ""

$ST_PYTHON bin/stereopy_ultimate_analysis.py

if [ $? -eq 0 ]; then
    echo ""
    echo "==========================================="
    echo "ANALYSIS COMPLETED SUCCESSFULLY!"
    echo "==========================================="
    echo "Complete results saved to: results_ultimate/"
    echo "==========================================="
    echo ""
else
    echo ""
    echo "==========================================="
    echo "ULTIMATE ANALYSIS FAILED!"
    echo "==========================================="
    echo "Check error logs for details"
    echo "==========================================="
    exit 1
fi

echo "==========================================="
echo "Final system status:"
free -h
echo "FINISHED:"
date
echo "==========================================="
