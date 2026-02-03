#!/bin/bash
#PBS -N STEREOPY_ANALYSIS
#PBS -P ji21
#PBS -q normalbw
#PBS -l walltime=01:00:00
#PBS -l ncpus=16
#PBS -l mem=128GB
#PBS -l jobfs=10GB
#PBS -l storage=scratch/ji21+gdata/ji21
#PBS -l wd
#PBS -o SPATIAL_ANALYSIS.out
#PBS -e SPATIAL_ANALYSIS.err

IMAGE="/g/data/ji21/users/pc1837/POLY_PIPELINE/stereopy_1.5.1.sif"

echo "==========================================="
echo "STEREOPY ANALYSIS"
echo "==========================================="
echo "START:"
date
echo "Job ID: $PBS_JOBID"
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

# Load environment with explicit paths
echo "Loading singularity"
module load singularity

echo "Verifying all dependencies"
singularity exec -B /g/data/ji21,/scratch/ji21 $IMAGE python3 -c "
# Function for logging stamps
import sys
from datetime import datetime
def log_step(message):
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f'[{timestamp}] {message}')
    sys.stdout.flush()
log_step(f'Python executable: {sys.executable}')
import stereo as st
log_step(f'Stereopy version: {st.__version__}')
import pandas as pd
log_step(f'Pandas version: {pd.__version__}')
import numpy as np
log_step(f'NumPy version: {np.__version__}')
import scipy
log_step(f'SciPy version: {scipy.__version__}')
import matplotlib
log_step(f'Matplotlib version: {matplotlib.__version__}')
import seaborn as sns
log_step(f'Seaborn version: {sns.__version__}')
"
echo ""

# Define analysis (1 for primary, 2 for secondary)
ANALYSIS=${ANALYSIS:-1}
if [[ ! "$ANALYSIS" =~ ^[0123]$ ]]; then
    echo "ERROR: ANALYSIS must be primary [1], secondary [2], network [3] or converter [0]"
    echo "Provided: ANALYSIS=$ANALYSIS"
    exit 1
fi

# Set thresholds
MIN_COUNTS=${MIN_COUNTS:-20}
MIN_GENES=${MIN_GENES:-3}
PCT_COUNTS_MT=${PCT_COUNTS_MT:-5}
N_PCS=${N_PCS:-10}
BIN_SIZE=${BIN_SIZE:-50}
MIN_X=${MIN_X:-}
MAX_X=${MAX_X:-}
MIN_Y=${MIN_Y:-}
MAX_Y=${MAX_Y:-}
HVG_MIN_MEAN=${HVG_MIN_MEAN:-0.0125}
HVG_MAX_MEAN=${HVG_MAX_MEAN:-3}
HVG_DISP=${HVG_DISP:-0.5}
HVG_TOP=${HVG_TOP:-2000}
INTEREST_GENES_PATH=${INTEREST_GENES_PATH:-"INPUT/interest_genes.txt"}
EXPRESSION_THR=${EXPRESSION_THR:-1.0}
INPUT_PATH=${INPUT_PATH:-""}
if [ -n "$INPUT_PATH" ]; then
    CLEAN_PATH="${INPUT_PATH%/}"
    
    if [ -d "$CLEAN_PATH" ]; then
        BASE_DIR="$CLEAN_PATH"
    else
        BASE_DIR=$(dirname "$CLEAN_PATH")
    fi
    
    OUTPUT_DIR="$BASE_DIR/CONVERTER_RESULTS"
else
    OUTPUT_DIR="./CONVERTER_RESULTS"
fi

echo "Filtering Parameters Used:"
echo "  Minimum counts: $MIN_COUNTS"
echo "  Minimum genes per cell: $MIN_GENES"
echo "  Mitochondrial percentage: $PCT_COUNTS_MT"
echo "  Bin size: $BIN_SIZE"
echo " HVG parameters: Min: $HVG_MIN_MEAN; Max: $HVG_MAX_MEAN; Dispersion: $HVG_DISP; Number of Top Genes: $HVG_TOP"
echo "Number of Principal Components used:"
echo "  $N_PCS"
echo "  This step can be inproved after first run. Ceck the Elbow Plot (PLOTS/QC/PCA_ELBOW.png) and insert the value of the elbow as N_PCS"
echo "You can alter the parameters inline:"
echo "  qsub -v ST_PYTHON="/home/user/.conda/envs/st/bin/python",MIN_COUNTS=50,MIN_GENES=5,PCT_COUNTS_MT=100,N_PCS=30 bin/2_DOC_ANALYSIS.sh"
echo ""

# Generate analysis script
echo "==========================================="
echo "Creating Primary Script!"
echo "==========================================="
cat > bin/SCRIPT_PRIMARY_ANALYSIS.py << EOF
#!/usr/bin/env python3
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
from matplotlib.colors import Normalize
import seaborn as sns
import psutil
import gc
from datetime import datetime
import sys
from collections import defaultdict
import json
import glob
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from scipy.sparse import issparse
import csv

# Set matplotlib backend for cluster
import matplotlib
matplotlib.use('Agg')

warnings.filterwarnings('ignore')

# Function for logging stamps
def log_step(message):
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()

# Function for monitoring memory usage
def log_memory_usage(step_name=""):
    try:
        memory = psutil.virtual_memory()
        swap = psutil.swap_memory()
        log_step(f"[MEMORY {step_name}] RAM: {memory.percent:.1f}% ({memory.used/1e9:.1f}GB/{memory.total/1e9:.1f}GB) | Swap: {swap.percent:.1f}%")
    except:
        log_step(f"[MEMORY {step_name}] Unable to get memory info")

# Function to create checkpoints
def save_progress_checkpoint(data, output_dir, checkpoint_name):
    checkpoint_dir = os.path.join(output_dir, 'CHECKPOINTS')
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    checkpoint_file = os.path.join(checkpoint_dir, f'{checkpoint_name.upper()}.json')
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
        log_step("=" * 100, file=sys.stderr)
        log_step(f"FATAL ERROR: The required input directory '{datasets_dir}' was not found.", file=sys.stderr)
        log_step("Please create this directory and place your .gef file inside.", file=sys.stderr)
        log_step("=" * 100, file=sys.stderr)
        sys.exit(1)
    gef_files = glob.glob(search_pattern)
    num_gef_files = len(gef_files)
    if num_gef_files == 0:
        log_step("=" * 100, file=sys.stderr)
        log_step(f"FATAL ERROR: No .gef file found in the '{datasets_dir}' directory.", file=sys.stderr)
        log_step("Please ensure there is exactly one .gef file to be analyzed.", file=sys.stderr)
        log_step("=" * 100, file=sys.stderr)
        sys.exit(1)
    if num_gef_files > 1:
        log_step("=" * 100, file=sys.stderr)
        log_step(f"FATAL ERROR: Multiple .gef files found in '{datasets_dir}'.", file=sys.stderr)
        log_step("Files found:", file=sys.stderr)
        for file_path in gef_files:
            log_step(f"  - {os.path.basename(file_path)}", file=sys.stderr)
        log_step("\nACTION REQUIRED: Please ensure ONLY the .gef file intended for analysis is in the directory.", file=sys.stderr)
        log_step("=" * 100, file=sys.stderr)
        sys.exit(1)
    data_path = gef_files[0]
    log_step(f"Successfully detected input file: {os.path.basename(data_path)}")
    log_step(f"Data path set to: {data_path}")
    
    return data_path

# Function for creating annotations by clustering method
def create_annotations_for_cluster_method(data, cluster_method='louvain'):
    annotations = {}
    
    if cluster_method not in data.cells:
        log_step(f"WARNING: Clustering method '{cluster_method}' not found in data")
        return annotations
    
    unique_clusters = sorted(data.cells[cluster_method].unique(), key=lambda x: int(x))
    cluster_sizes = data.cells[cluster_method].value_counts().sort_index()
    
    # 1. Alphabetical
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    alphabetical = {}
    for i, cluster in enumerate(unique_clusters):
        if i < len(alphabet):
            alphabetical[str(cluster)] = alphabet[i]
        else:
            alphabetical[str(cluster)] = f"cluster_{cluster}"
    
    # 2. Numerical with size description
    numerical = {}
    for cluster in unique_clusters:
        size = cluster_sizes[cluster]
        if size > 1000:
            size_desc = "large"
        elif size > 100:
            size_desc = "medium"
        else:
            size_desc = "small"
        numerical[str(cluster)] = f"{cluster_method.capitalize()}_{cluster}_{size_desc}"
    
    annotations['alphabetical'] = alphabetical
    annotations['numerical'] = numerical
    
    return annotations


# Function for analysis of genes of interest
def create_custom_gene_markers():
    file_path = "$INTEREST_GENES_PATH"
    gene_markers = {}

    if not os.path.exists(file_path):
        log_step(f"WARNING: File of genes of interest not found. Proceeding without genes of interest.")
        return gene_markers

    try:
        with open(file_path, mode='r', newline='', encoding='utf-8') as file:
            reader = csv.reader(file)
            
            try:
                next(reader) 
            except StopIteration:
                log_step(f"WARNING: The file {file_path} is empty. Proceeding without genes of interest")
                return gene_markers

            for row in reader:
                if not any(row):
                    continue

                marker_name = row[0].strip()
                loc_ids = [loc_id.strip() for loc_id in row[1:] if loc_id.strip()]
                
                if marker_name and loc_ids:
                    gene_markers[marker_name] = loc_ids
                elif marker_name and not loc_ids:
                    log_step(f"WARNING: The gene '{marker_name}' was ignored because it has no gene IDs.")

    except Exception as e:
        log_step(f"ERROR: Processing {file_path} failed. Proceeding without genes of interest. Error: {e}")
        return {}
        
    return gene_markers

# Function for analyzing custom genes
def apply_gene_interest_annotation(data, gene_markers, cluster_res_key='leiden', threshold=1.2):

    available_genes = set(data.genes.gene_name.tolist())
    cluster_annotations = {}
    cluster_scores = {}
    annotation_details = {}
    
    log_step("Checking gene availability")
    for category, genes in gene_markers.items():
        available_in_category = [gene for gene in genes if gene in available_genes]
        missing_in_category = [gene for gene in genes if gene not in available_genes]
        
        log_step(f"{category}: {len(available_in_category)}/{len(genes)} genes found")
        if missing_in_category:
            log_step(f"  Missing: {', '.join(missing_in_category)}")
    
    # For each cluster, calculate expression scores
    for cluster in data.cells[cluster_res_key].unique():
        cluster_mask = data.cells[cluster_res_key] == cluster
        cluster_cells_count = cluster_mask.sum()
        
        best_score = 0
        best_category = None
        category_scores = {}
        
        for category, genes in gene_markers.items():
            if not genes:
                continue
                
            valid_genes = [gene for gene in genes if gene in available_genes]
            if not valid_genes:
                continue
            
            # Calculate expression scores for this category
            gene_scores = []
            
            for gene in valid_genes:
                try:
                    gene_idx = data.genes.gene_name.tolist().index(gene)
                    
                    # Expression in this cluster
                    cluster_expr = data.exp_matrix[cluster_mask, gene_idx].mean()
                    
                    # Expression in other clusters
                    other_mask = ~cluster_mask
                    if other_mask.sum() > 0:
                        other_expr = data.exp_matrix[other_mask, gene_idx].mean()
                        
                        # Calculate fold change
                        if other_expr > 0:
                            fold_change = cluster_expr / other_expr
                        else:
                            fold_change = cluster_expr + 1
                    else:
                        fold_change = cluster_expr + 1
                    
                    gene_scores.append(fold_change)
                    
                except (ValueError, IndexError):
                    continue
            
            # Average score for this category
            if gene_scores:
                avg_score = sum(gene_scores) / len(gene_scores)
                category_scores[category] = avg_score
                
                if avg_score > best_score and avg_score > threshold:
                    best_score = avg_score
                    best_category = category
        
        # Assign annotation
        if best_category:
            cluster_annotations[str(cluster)] = f"{best_category}_enriched_cluster_{cluster}"
        else:
            cluster_annotations[str(cluster)] = f"Cluster_{cluster}"
        
        cluster_scores[str(cluster)] = best_score
        annotation_details[str(cluster)] = {
            'best_category': best_category,
            'best_score': best_score,
            'all_scores': category_scores,
            'cell_count': int(cluster_cells_count)
        }
    
    return cluster_annotations, cluster_scores, annotation_details

# Function for biological generic annotation
def create_biological_annotation_dict(marker_genes_df, n_top_genes=5):
    annotations = {}
    
    if 'cluster' in marker_genes_df.columns and 'genes' in marker_genes_df.columns:
        for cluster in marker_genes_df['cluster'].unique():
            cluster_markers = marker_genes_df[
                marker_genes_df['cluster'] == cluster
            ]
            
            # Get top marker genes for this cluster
            if 'scores' in cluster_markers.columns:
                top_genes = cluster_markers.nlargest(n_top_genes, 'scores')['genes'].tolist()
            else:
                top_genes = cluster_markers.head(n_top_genes)['genes'].tolist()
            
            # Create annotation based on top 5 genes
            gene_string = "_".join(top_genes[:5])
            annotations[str(cluster)] = f"Cluster_{cluster}_{gene_string}"
    
    return annotations

# Function for direct gene visualization
def create_direct_gene_visualization():
    try:
        log_step("Extracting data directly for gene-level visualization")
        spatial_coords = None
        if hasattr(data.cells, 'obs'):
            obs_df = data.cells.obs
            coord_cols = [col for col in obs_df.columns if any(term in col.lower() for term in ['x', 'y', 'spatial', 'coord', 'position'])]
            if len(coord_cols) >= 2:
                spatial_coords = obs_df[coord_cols[:2]].values
                log_step(f"Using coordinate columns from data.cells.obs: {coord_cols[:2]}")
            else:
                numeric_cols = obs_df.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) >= 2:
                    spatial_coords = obs_df[numeric_cols[:2]].values
                    log_step(f"Using numeric columns from data.cells.obs: {numeric_cols[:2].tolist()}")
        
        if spatial_coords is None and hasattr(data, 'position') and data.position is not None:
            spatial_coords = data.position
            log_step("Using data.position coordinates.")
        
        if spatial_coords is None:
            raise Exception("Could not find spatial coordinates.")
        
        log_step(f"Final spatial coordinates shape: {spatial_coords.shape}")
        
        base_clusters = data.cells['louvain']
        
        exp_matrix = data.exp_matrix
        if issparse(exp_matrix):
            exp_matrix = exp_matrix.toarray()
        
        gene_names = data.genes.gene_name.tolist()
        gene_markers = create_custom_gene_markers()
        print(gene_markers)

        log_step("Calculating gene enrichment by category")
        enrichment_results = {}
        cell_categories = pd.Series(['Other'] * len(base_clusters), index=base_clusters.index)
        enriched_clusters = defaultdict(list)
        
        for cluster_id in base_clusters.unique():
            cluster_mask = base_clusters == cluster_id
            cluster_cells = cluster_mask.sum()
            
            if cluster_cells < 10:
                continue
            
            cluster_scores = {}
            for category, genes in gene_markers.items():
                valid_genes_indices = [gene_names.index(gene) for gene in genes if gene in gene_names]
                if not valid_genes_indices:
                    continue
                
                cluster_expr = exp_matrix[cluster_mask, :][:, valid_genes_indices].mean()
                other_expr = exp_matrix[~cluster_mask, :][:, valid_genes_indices].mean()
                
                if other_expr > 0:
                    score = cluster_expr / other_expr
                else:
                    score = cluster_expr + 1
                
                cluster_scores[category] = float(score)
            
            enrichment_results[cluster_id] = cluster_scores
            
            for category, score in cluster_scores.items():
                if score > 1.2:
                    if cluster_id not in enriched_clusters[category]:
                        enriched_clusters[category].append(cluster_id)
                        
        for category, clusters in enriched_clusters.items():
            for cluster_id in clusters:
                mask = base_clusters == cluster_id
                cell_categories[mask] = category

        log_step(f"Calculated enrichment for {len(enrichment_results)} clusters")
        
        log_step("Creating combined category visualization.")
        
        enriched_categories_list = [cat for cat in sorted(enriched_clusters.keys())]
        colors = plt.cm.tab10(np.linspace(0, 1, max(1, len(enriched_categories_list))))
        
        fig, ax = plt.subplots(1, 1, figsize=(12, 10))
        if 'Other' in cell_categories.unique():
            other_mask = cell_categories == 'Other'
            ax.scatter(spatial_coords[other_mask, 0], spatial_coords[other_mask, 1],
                       c='#CCCCCC', s=5.0, alpha=0.4, label='Other clusters')
        
        for i, category in enumerate(enriched_categories_list):
            mask = cell_categories == category
            count = mask.sum()
            ax.scatter(spatial_coords[mask, 0], spatial_coords[mask, 1],
                       c=[colors[i]], s=8.0, alpha=0.95,
                       label=f'{category.upper()}: {count:,} cells')
        
        ax.set_title('Gene Interest Spatial Distribution\n(Dynamic Gene Categories)', fontweight='bold', fontsize=14)
        ax.set_xlabel('Spatial X')
        ax.set_ylabel('Spatial Y')
        ax.legend(markerscale=3, fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        plt.tight_layout(rect=[0, 0, 0.8, 1])
        main_plot = os.path.join(output_dir, 'PLOTS', 'INTEREST_ANALYSIS', 'INTEREST_ANALYSIS_COMBINED.png')
        plt.savefig(main_plot, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        log_step(f"Main combined plot saved: {main_plot}")

        log_step("Creating individual plots for each gene of interest")
        individual_gene_plot_dir = os.path.join(output_dir, 'PLOTS', 'INTEREST_ANALYSIS')
        os.makedirs(individual_gene_plot_dir, exist_ok=True)
        
        all_genes_of_interest = []
        for category, genes in gene_markers.items():
            all_genes_of_interest.extend(genes)
        
        for gene_name in sorted(list(set(all_genes_of_interest))):
            if gene_name not in gene_names:
                log_step(f"WARNING: Gene '{gene_name}' not found in expression data. Skipping individual plot.")
                continue
            
            gene_idx = gene_names.index(gene_name)
            gene_expression = exp_matrix[:, gene_idx]
            
            normalized_expression = gene_expression / (gene_expression.max() if gene_expression.max() > 0 else 1)
            
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            sc = ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                            c=gene_expression,
                            cmap='viridis',
                            s=6, alpha=0.8, vmin=0, vmax=gene_expression.max())
            
            ax.set_title(f'Spatial Expression of Gene: {gene_name}',
                         fontweight='bold', fontsize=14)
            ax.set_xlabel('Spatial X (μm)')
            ax.set_ylabel('Spatial Y (μm)')
            #ax.set_aspect('equal', adjustable='box')
            
            cbar = fig.colorbar(sc, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
            cbar.set_label('Gene Expression Level', rotation=270, labelpad=15)
            
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            
            plot_file = os.path.join(individual_gene_plot_dir, f'SPATIAL_EXPRESSION_{gene_name}.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            
            log_step(f"Individual plot for gene '{gene_name}' saved: {plot_file}")
            log_step(f"File size: {os.path.getsize(plot_file) / 1e6:.1f} MB")
            
        report_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'ENRICHMENT_ANALYSIS_REPORT.txt')
        with open(report_file, 'w') as f:
            f.write("GENE ENRICHMENT ANALYSIS REPORT\n")
            f.write("="*100 + "\n\n")
            f.write(f"Date: {datetime.now()}\n")
            f.write(f"Total cells: {len(cell_categories):,}\n")
            f.write(f"Total clusters: {base_clusters.nunique()}\n")
            f.write(f"Enrichment threshold: 1.2\n\n")
            
            f.write("GENE CATEGORIES:\n")
            for category, genes in gene_markers.items():
                f.write(f"\n{category.upper()}:\n")
                if genes:
                    for gene in genes:
                        f.write(f"  - {gene}\n")
                else:
                    f.write("  - (No genes defined)\n")
            
            f.write(f"\nENRICHMENT RESULTS BY CATEGORY:\n")
            f.write("="*100 + "\n")
            
            for category in enriched_categories_list:
                count = cell_categories.value_counts().get(category, 0)
                percentage = (count / len(cell_categories) * 100) if len(cell_categories) > 0 else 0
                clusters = enriched_clusters.get(category, [])
                
                f.write(f"\n{category.upper()}:\n")
                f.write(f"  Enriched cells: {count:,} ({percentage:.1f}%)\n")
                f.write(f"  Enriched clusters: {len(clusters)}\n")
                f.write(f"  Cluster IDs: {clusters}\n")
        
        log_step(f"Report saved: {report_file}")
        
        return {
            'spatial_coords': spatial_coords,
            'cell_categories': cell_categories,
            'enriched_clusters': enriched_clusters,
            'enrichment_results': enrichment_results
        }
    
    except Exception as e:
        log_step(f"ERROR in direct visualization: {e}")
        import traceback
        traceback.print_exc()
        return None

# Spatial visualization
def create_clean_spatial_visualization_from_direct_results(direct_results, data, output_dir):
    try:
        if not direct_results:
            log_step("No direct results available for spatial visualization")
            return
        
        log_step("Creating clean spatial visualization from direct analysis results")
        
        spatial_coords = direct_results['spatial_coords']
        cell_categories = direct_results['cell_categories']
        enriched_clusters = direct_results['enriched_clusters']
        stereopy_spatial = None
        if hasattr(data.cells, 'obs'):
            obs_df = data.cells.obs
            numeric_cols = obs_df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) >= 2:
                stereopy_spatial = obs_df[numeric_cols[:2]].values
        log_step(f"Using direct analysis spatial coordinates: {spatial_coords.shape}")
        fig, ax = plt.subplots(1, 1, figsize=(14, 12))
        other_mask = cell_categories == 'Other'
        if other_mask.sum() > 0:
            ax.scatter(spatial_coords[other_mask, 0], spatial_coords[other_mask, 1],
                      c='whitesmoke', s=5, alpha=0.4, edgecolors='lightgray', 
                      linewidths=0.05, rasterized=True)
        enriched_categories = [cat for cat in cell_categories.unique() if cat != 'Other']
        category_colors = {}
        base_colors = ['#d62728', '#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#8c564b']
        for i, category in enumerate(enriched_categories):
            category_colors[category] = base_colors[i % len(base_colors)]
        for category in enriched_categories:
            mask = cell_categories == category
            count = mask.sum()
            percentage = count / len(cell_categories) * 100
            ax.scatter(spatial_coords[mask, 0], spatial_coords[mask, 1],
                      c=category_colors[category], s=8, alpha=0.9,
                      label=f'{category.upper()}: {count:,} cells ({percentage:.1f}%)',
                      rasterized=True)
        ax.set_title('Spatial Gene Interest Enrichment\n(Clean View)', 
                    fontweight='bold', fontsize=16, pad=20)
        ax.set_xlabel('Spatial X (μm)', fontsize=12)
        ax.set_ylabel('Spatial Y (μm)', fontsize=12)
        ax.legend(frameon=False, fontsize=11, markerscale=2, 
                 bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #ax.set_aspect('equal', adjustable='box')
        plt.tight_layout()
        
        clean_plot = os.path.join(output_dir, 'PLOTS', 'INTEREST_ANALYSIS', 'INTEREST_ANALYSIS_CLEAN_SPATIAL.png')
        plt.savefig(clean_plot, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        log_step(f"Clean spatial plot saved: {clean_plot}")
        log_step(f"File size: {os.path.getsize(clean_plot) / 1e6:.1f} MB")
        
        # Genes of interest visualization
        fig, ax = plt.subplots(1, 1, figsize=(12, 10))
        for category in enriched_categories:
            mask = cell_categories == category
            count = mask.sum()
            
            ax.scatter(spatial_coords[mask, 0], spatial_coords[mask, 1],
                      c=category_colors[category], s=6, alpha=0.95,
                      label=f'{category.upper()}: {count:,} cells',
                      rasterized=True)
        
        ax.set_title('Gene of Interest Regions', 
                    fontweight='bold', fontsize=16, pad=20)
        ax.set_xlabel('Spatial X (μm)', fontsize=12)
        ax.set_ylabel('Spatial Y (μm)', fontsize=12)
        
        ax.legend(frameon=False, fontsize=12, markerscale=2)
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #ax.set_aspect('equal', adjustable='box')
        
        plt.tight_layout()
        
        genes_only_plot = os.path.join(output_dir, 'PLOTS', 'INTEREST_ANALYSIS', 'INTEREST_ANALYSIS_ONLY_SPATIAL.png')
        plt.savefig(genes_only_plot, dpi=300, bbox_inches='tight', facecolor='white')
        
        genes_only_pdf = os.path.join(output_dir, 'PLOTS', 'INTEREST_ANALYSIS', 'INTEREST_ANALYSIS_ONLY_SPATIAL.pdf')
        plt.savefig(genes_only_pdf, bbox_inches='tight', facecolor='white')
        plt.close()
        
        log_step(f"Genes-only plot saved: {genes_only_plot} and {genes_only_pdf}")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        if other_mask.sum() > 0:
            ax1.scatter(spatial_coords[other_mask, 0], spatial_coords[other_mask, 1],
                       c='#CCCCCC', s=2.0, alpha=0.4, rasterized=True)
        
        for category in enriched_categories:
            mask = cell_categories == category
            ax1.scatter(spatial_coords[mask, 0], spatial_coords[mask, 1],
                       c=category_colors[category], s=4.0, alpha=0.95,
                       label=f'{category.upper()}', rasterized=True)
        
        ax1.set_title('With Context\n(All Clusters)', fontweight='bold', fontsize=14)
        ax1.set_xlabel('Spatial X (μm)')
        ax1.set_ylabel('Spatial Y (μm)')
        ax1.legend(frameon=False, fontsize=10, markerscale=2)
        ax1.grid(True, alpha=0.2)
        try:
            x_range = float(spatial_coords[:, 0].max()) - float(spatial_coords[:, 0].min())
            y_range = float(spatial_coords[:, 1].max()) - float(spatial_coords[:, 1].min())
            aspect_ratio = x_range / y_range if y_range > 0 else 1
        except (ValueError, TypeError):
            aspect_ratio = 1  # fallback se conversão falhar
        ax1.set_aspect(aspect_ratio, adjustable='box')
        #ax1.set_aspect('equal', adjustable='box')
        for category in enriched_categories:
            mask = cell_categories == category
            ax2.scatter(spatial_coords[mask, 0], spatial_coords[mask, 1],
                       c=category_colors[category], s=4, alpha=0.95,
                       label=f'{category.upper()}', rasterized=True)
        
        ax2.set_title('Gene Interest Only\n(Clean View)', fontweight='bold', fontsize=14)
        ax2.set_xlabel('Spatial X (μm)')
        ax2.set_ylabel('Spatial Y (μm)')
        ax2.legend(frameon=False, fontsize=10, markerscale=2)
        ax2.grid(True, alpha=0.2)
        ax2.set_aspect(aspect_ratio, adjustable='box')
        #ax2.set_aspect('equal', adjustable='box')
        
        for ax in [ax1, ax2]:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        
        comparison_plot = os.path.join(output_dir, 'PLOTS', 'INTEREST_ANALYSIS', 'INTEREST_ANALYSIS_COMPARISON_SPATIAL.png')
        plt.savefig(comparison_plot, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        log_step(f"Comparison plot saved: {comparison_plot}")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        categories_for_plot = enriched_categories + ['Other']
        counts_for_plot = [cell_categories.value_counts().get(cat, 0) for cat in categories_for_plot]
        colors_for_plot = [category_colors.get(cat, 'lightgray') for cat in categories_for_plot]
        
        bars = ax1.bar(range(len(categories_for_plot)), counts_for_plot, 
                      color=colors_for_plot, alpha=0.8)
        
        ax1.set_xticks(range(len(categories_for_plot)))
        ax1.set_xticklabels([cat.upper() for cat in categories_for_plot], rotation=45)
        ax1.set_ylabel('Number of Cells')
        ax1.set_title('Cell Distribution by Gene Category', fontweight='bold')
        
        total_cells = sum(counts_for_plot)
        for bar, count in zip(bars, counts_for_plot):
            percentage = count / total_cells * 100
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts_for_plot)*0.01,
                    f'{count:,}\n({percentage:.1f}%)', ha='center', va='bottom', fontweight='bold')
        
        sizes = [cell_categories.value_counts().get(cat, 0) for cat in enriched_categories]
        colors_pie = [category_colors[cat] for cat in enriched_categories]
        labels_pie = [f'{cat.upper()}\n({cell_categories.value_counts().get(cat, 0):,} cells)' 
                     for cat in enriched_categories]
        
        if sum(sizes) > 0:
            wedges, texts, autotexts = ax2.pie(sizes, labels=labels_pie, colors=colors_pie,
                                              autopct='%1.1f%%', startangle=90)
            ax2.set_title('Gene Interest Categories\n(Enriched Cells Only)', fontweight='bold')
        
        plt.tight_layout()
        
        stats_plot = os.path.join(output_dir, 'PLOTS', 'INTEREST_ANALYSIS', 'INTEREST_ANALYSIS_ONLY_SPATIAL_STATISTCS.png')
        plt.savefig(stats_plot, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        log_step(f"Statistics plot saved: {stats_plot}")
        
        for category in enriched_categories:
            count = cell_categories.value_counts().get(category, 0)
            percentage = count / len(cell_categories) * 100
            clusters = len(enriched_clusters.get(category, []))
            log_step(f"  {category.upper()}: {count:,} cells ({percentage:.1f}%) in {clusters} clusters")
        
        return True
        
    except Exception as e:
        log_step(f"ERROR in clean spatial visualization: {e}")
        import traceback
        traceback.print_exc()
        return False

# Script Initialization
start_time = datetime.now()
log_step("Starting analysis")
log_memory_usage("START")
timestamp = start_time.strftime("%Y%m%d_%H%M")

# Filepaths
data_path = get_single_gef_file()
sample_id = os.path.splitext(os.path.basename(data_path))[0].replace('.', '_')
output_dir = f'RESULTS_{sample_id}_{timestamp}'
log_step(f"RESULT_FOLDER_PATH:{output_dir}")

# Directory structure
directories = [
    output_dir,
    os.path.join(output_dir, 'PLOTS'),
    os.path.join(output_dir, 'PLOTS', 'QC'),
    os.path.join(output_dir, 'PLOTS', 'INTEREST_ANALYSIS'),
    os.path.join(output_dir, 'PLOTS', 'ANNOTATION'),
    os.path.join(output_dir, 'INTEREST_ANALYSIS'),
    os.path.join(output_dir, 'INTEREST_ANALYSIS', 'COMPLETE'),
    os.path.join(output_dir, 'INTEREST_ANALYSIS', 'FILTERED'),
    os.path.join(output_dir, 'STATISTICAL_ANALYSIS'),
    os.path.join(output_dir, 'LOGS'),
    os.path.join(output_dir, 'EXPORTS'),
    os.path.join(output_dir, 'CLUSTERING'),
    os.path.join(output_dir, 'CLUSTERING', 'LOUVAIN'),
    os.path.join(output_dir, 'CLUSTERING', 'LEIDEN'),
    os.path.join(output_dir, 'CLUSTERING', 'SPATIAL_LEIDEN'),
    os.path.join(output_dir, 'CHECKPOINTS')
]

for directory in directories:
    os.makedirs(directory, exist_ok=True)

log_step(f"Created {len(directories)} output folders")

# Analysis log
analysis_log_path = os.path.join(output_dir, 'LOGS', 'ANALYSIS_LOG.txt')
with open(analysis_log_path, 'w') as f:
    f.write("="*100 + "\n")
    f.write("ANALYSIS LOG\n")
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
    #st.io.read_gef_info(data_path)
    data = st.io.read_gef(file_path=data_path, bin_size=int(float("$BIN_SIZE")))
    m_x, M_x = os.getenv('MIN_X'), os.getenv('MAX_X')
    m_y, M_y = os.getenv('MIN_Y'), os.getenv('MAX_Y')

    if all(v and v.strip() for v in [m_x, M_x, m_y, M_y]):
        log_step(f"Applying manual crop: X[{m_x}:{M_x}], Y[{m_y}:{M_y}]")
        data.sub_by_position(
            x_min=int(m_x), x_max=int(M_x), 
            y_min=int(m_y), y_max=int(M_y)
        )

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
N_PCS = $N_PCS
BIN_SIZE = $BIN_SIZE

# General statistics
initial_stats_file = os.path.join(output_dir, 'LOGS', 'DATA_STATS.txt')
with open(initial_stats_file, 'w') as f:
    f.write("="*100 + "\n")
    f.write("DATA STATISTICS\n")
    f.write("="*100 + "\n\n")
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

# Cell filtering with detailed statistics
log_step("(LOUVAIN) Performing cell filtering")
cells_before = data.n_cells
genes_before = data.n_genes
try:
    data.tl.filter_cells(min_counts=MIN_COUNTS)
    if hasattr(data, 'position'):
        x_coords = data.position[:, 0]
        y_coords = data.position[:, 1]
        actual_min_x = x_coords.min()
        actual_max_x = x_coords.max()
        actual_min_y = y_coords.min()
        actual_max_y = y_coords.max()
        log_step(f"Detected coordinates: X[{actual_min_x}, {actual_max_x}], Y[{actual_min_y}, {actual_max_y}]")
    MIN_X = int(os.getenv('MIN_X', '0')) if os.getenv('MIN_X') else None
    MAX_X = int(os.getenv('MAX_X', '0')) if os.getenv('MAX_X') else None
    MIN_Y = int(os.getenv('MIN_Y', '0')) if os.getenv('MIN_Y') else None
    MAX_Y = int(os.getenv('MAX_Y', '0')) if os.getenv('MAX_Y') else None
    if all(v is not None for v in [MIN_X, MAX_X, MIN_Y, MAX_Y]):
        log_step(f"Filtering positions: X[{MIN_X}, {MAX_X}], Y[{MIN_Y}, {MAX_Y}]")
        data.tl.filter_coordinates(min_x=MIN_X, max_x=MAX_X, min_y=MIN_Y, max_y=MAX_Y)
    else:
        log_step("Skipping coordinate filtering")

    cells_after = data.n_cells
    genes_after = data.n_genes
    
    log_step(f"(LOUVAIN) Filtering complete: {cells_before} -> {cells_after} cells ({cells_before-cells_after} removed)")
    log_step(f"(LOUVAIN) Genes retained: {genes_after} (from {genes_before})")
    
    # Update stats
    with open(initial_stats_file, 'a') as f:
        f.write("AFTER FILTERING:\n")
        f.write(str(data) + "\n")
        f.write(f"Cells removed: {cells_before - cells_after} ({(cells_before-cells_after)/cells_before*100:.1f}%)\n")
        f.write(f"Genes retained: {genes_after}\n\n")
    
    save_progress_checkpoint(data, output_dir, 'cells_filtered')
    log_memory_usage("FILTERED")
    
except Exception as e:
    log_step(f"(LOUVAIN) ERROR in cell filtering: {e}")
    sys.exit(1)

# QC plots before normalization
log_step("(LOUVAIN) Generating QC plots before normalization")
qc_dir = os.path.join(output_dir, 'PLOTS', 'QC')
try:
    data.plt.violin(out_path=os.path.join(qc_dir, 'PRE_NORM_VIOLIN.png'))
    data.plt.spatial_scatter(out_path=os.path.join(qc_dir, 'PRE_NORM_SPATIAL_SCATTER.png'))
    data.plt.genes_count(out_path=os.path.join(qc_dir, 'PRE_NORM_COUNT.png'))
    log_step("(LOUVAIN) Pre-normalization QC plots saved")
except Exception as e:
    log_step(f"(LOUVAIN) WARNING: Error in pre-normalization plots: {e}")

# Normalization
log_step("(LOUVAIN) Performing normalization and log transformation")
try:
    data.tl.normalize_total(target_sum=None)
    data.tl.log1p()
    save_progress_checkpoint(data, output_dir, 'normalized')
    log_memory_usage("NORMALIZED")
except Exception as e:
    log_step(f"(LOUVAIN) ERROR in normalization: {e}")
    sys.exit(1)

# QC plots after normalization
log_step("(LOUVAIN) Generating QC plots after normalization")
try:
    data.plt.violin(out_path=os.path.join(qc_dir, 'POS_NORM_VIOLIN.png'))
    data.plt.spatial_scatter(out_path=os.path.join(qc_dir, 'POS_NORM_SPATIAL_SCATTER.png'))
    data.plt.genes_count(out_path=os.path.join(qc_dir, 'POS_NORM_COUNT.png'))
    log_step("(LOUVAIN) Post-normalization QC plots saved")
except Exception as e:
    log_step(f"(LOUVAIN) WARNING: Error in post-normalization plots: {e}")

# Create raw checkpoint
log_step("(LOUVAIN) Creating raw data checkpoint")
try:
    data.tl.raw_checkpoint()
    save_progress_checkpoint(data, output_dir, 'RAW_CHECKPOINT')
    log_step("(LOUVAIN) Raw checkpoint created")
except Exception as e:
    log_step(f"(LOUVAIN) ERROR creating raw checkpoint: {e}")
    sys.exit(1)

# Highly variable genes identification
log_step("(GENERAL) Identifying highly variable genes")
try:
    data.tl.highly_variable_genes(
        min_mean=$HVG_MIN_MEAN,
        max_mean=$HVG_MAX_MEAN,
        min_disp=$HVG_DISP,
        n_top_genes=$HVG_TOP,
        res_key='highly_variable_genes'
    )

    data.plt.highly_variable_genes(
        res_key='highly_variable_genes',
        out_path=os.path.join(qc_dir, 'HIGHLY_VARIABLE_GENES.png')
    )
    
    # Save HVG list
    if 'highly_variable_genes' in data.tl.result:
        hvg_file = os.path.join(output_dir, 'EXPORTS', 'HIGHLY_VARIABLE_GENES.csv')
        hvg_data = data.tl.result['highly_variable_genes']
        if isinstance(hvg_data, pd.DataFrame):
            hvg_data.to_csv(hvg_file, index=False)
            log_step(f"(GENERAL) Highly variable genes list saved ({len(hvg_data)} genes)")
    
    save_progress_checkpoint(data, output_dir, 'POS_HVG')
    log_memory_usage("HVG_IDENTIFIED")
    
except Exception as e:
    log_step(f"(LOUVAIN) ERROR in highly variable genes: {e}")
    sys.exit(1)

# Data scaling
log_step("(LOUVAIN) Performing data scaling")
try:
    data.tl.scale(max_value=10, zero_center=False)
    save_progress_checkpoint(data, output_dir, 'SCALED')
except Exception as e:
    log_step(f"(LOUVAIN) ERROR in scaling: {e}")
    sys.exit(1)

# PCA
N_PCS = $N_PCS
log_step("(LOUVAIN) Performing PCA analysis")
try:
    data.tl.pca(
        use_highly_genes=False,
        n_pcs=N_PCS,
        res_key='pca'
    )
    data.tl.key_record
    data.plt.elbow(
        pca_res_key='pca',
        out_path=os.path.join(qc_dir, 'PCA_ELBOW.png')
    )
    
    # Save PCA results
    if 'pca_variance_ratio' in data.tl.result:
        pca_file = os.path.join(output_dir, 'EXPORTS', 'PCA_VARIANCE_RATIOS.csv')
        pca_data = pd.DataFrame({
            'PC': range(1, len(data.tl.result['pca_variance_ratio']) + 1),
            'variance_ratio': data.tl.result['pca_variance_ratio'],
            'cumulative_variance': np.cumsum(data.tl.result['pca_variance_ratio'])
        })
        pca_data.to_csv(pca_file, index=False)
        log_step(f"(LOUVAIN) PCA variance ratios saved ({len(pca_data)} components)")
    
    save_progress_checkpoint(data, output_dir, 'POS_PCA')
    log_memory_usage("PCA_COMPLETE")
    
except Exception as e:
    log_step(f"(LOUVAIN) ERROR in PCA: {e}")
    sys.exit(1)

# Neighborhood graph computation
log_step("(LOUVAIN) Computing neighborhood graph")
try:
    data.tl.neighbors(
        pca_res_key='pca',
        n_pcs=N_PCS,
        res_key='neighbors'
    )
    save_progress_checkpoint(data, output_dir, 'neighbors_computed')
    log_memory_usage("NEIGHBORS_COMPUTED")
except Exception as e:
    log_step(f"(LOUVAIN) ERROR in neighbors computation: {e}")
    sys.exit(1)

# UMAP computation
log_step("(LOUVAIN) Computing UMAP embedding")
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
            
            umap_file = os.path.join(output_dir, 'EXPORTS', 'UMAP_COORDINATES.csv')
            umap_df.to_csv(umap_file, index=False)
            log_step("(LOUVAIN) UMAP coordinates saved")
    
    save_progress_checkpoint(data, output_dir, 'umap_complete')
    log_memory_usage("UMAP_COMPLETE")
    
except Exception as e:
    log_step(f"(LOUVAIN) ERROR in UMAP: {e}")
    sys.exit(1)

# Spatial neighbors for spatial clustering
log_step("(LOUVAIN) Computing spatial neighborhood graph")
try:
    data.tl.spatial_neighbors(
        neighbors_res_key='neighbors',
        res_key='spatial_neighbors'
    )
    log_step("(LOUVAIN) Spatial neighbors computed successfully")
    save_progress_checkpoint(data, output_dir, 'spatial_neighbors_computed')
    log_memory_usage("SPATIAL_NEIGHBORS_COMPUTED")
except Exception as e:
    log_step(f"(LOUVAIN) WARNING: Error in spatial neighbors computation: {e}")

# Clustering
log_step("(LOUVAIN) Performing clustering")
try:
    data.tl.louvain(
        neighbors_res_key='neighbors',
        res_key='louvain'
    )
    
    # Generate cluster plots
    log_step("(LOUVAIN) Generating cluster visualization plots")
    data.plt.cluster_scatter(
        res_key='louvain',
        out_path=os.path.join(output_dir,'CLUSTERING', 'LOUVAIN_SPATIAL_CLUSTERS.png')
    )
# Spatial neighbors for spatial clustering
    log_step("(LOUVAIN) Computing spatial neighborhood graph")
    try:
        data.tl.spatial_neighbors(
            neighbors_res_key='neighbors',
            res_key='spatial_neighbors'
        )
        log_step("(LOUVAIN) Spatial neighbors computed successfully")
        save_progress_checkpoint(data, output_dir, 'spatial_neighbors_computed')
        log_memory_usage("SPATIAL_NEIGHBORS_COMPUTED")
    except Exception as e:
        log_step(f"(LOUVAIN) WARNING: Error in spatial neighbors computation: {e}")

    data.plt.umap(
        res_key='umap',
        cluster_key='louvain',
        out_path=os.path.join(output_dir,'CLUSTERING', 'LOUVAIN_UMAP_CLUSTERS.png')
    )
    log_step("(LOUVAIN) Generating individual cluster plots")
    unique_clusters = data.cells['louvain'].unique()
    n_clusters = len(unique_clusters)
    log_step(f"(LOUVAIN) Found {n_clusters} clusters. Generating individual plots")
    
except Exception as e:
    log_step(f"(LOUVAIN) ERROR in clustering: {e}")
    import traceback
    log_step(f"(LOUVAIN) Full traceback: {traceback.format_exc()}")
    sys.exit(1)

## MARKER GENE ANALYSIS
## This step process all genes and all clusters (requires HPC)
log_step("="*100)
log_step("(LOUVAIN) MARKER GENE ANALYSIS")
log_step("="*100)

marker_start_time = datetime.now()
log_step(" (LOUVAIN) Finding marker genes using t-test method")

try:
    data.tl.find_marker_genes(
        cluster_res_key='louvain',
        method='t_test',
        use_highly_genes=False,
        use_raw=True,
        res_key='marker_genes'
    )

    marker_end_time = datetime.now()
    marker_duration = (marker_end_time - marker_start_time).total_seconds()
    log_step(f"(LOUVAIN) Marker gene analysis completed in {marker_duration:.1f} seconds")
    log_memory_usage("MARKER_GENES_FOUND")
    
except Exception as e:
    log_step(f"(LOUVAIN) ERROR in marker gene analysis: {e}")
    sys.exit(1)

# Marker genes processing
marker_genes_dir = os.path.join(output_dir, 'INTEREST_ANALYSIS')
complete_results_dir = os.path.join(marker_genes_dir, 'COMPLETE')

if 'marker_genes' in data.tl.result:
    try:
        marker_result = data.tl.result['marker_genes']
        log_step(f"(LOUVAIN) Processing marker genes results with {len(marker_result)} components")
        
        # Get all cluster comparison keys
        cluster_keys = [k for k in marker_result.keys() if k.endswith('.vs.rest')]
        log_step(f"(LOUVAIN) Processing ALL {len(cluster_keys)} cluster comparisons")
        
        # Process ALL clusters with ALL genes
        all_marker_genes = []
        cluster_summaries = []
        processing_stats = []
        
        log_step("(LOUVAIN) Processing all clusters with complete gene sets")
        
        for i, cluster_key in enumerate(cluster_keys):
            cluster_num = cluster_key.split('.')[0]
            cluster_start_time = datetime.now()
            
            try:
                cluster_data = marker_result[cluster_key]
                
                if isinstance(cluster_data, pd.DataFrame) and len(cluster_data) > 0:
                    # Add cluster identifier
                    cluster_data = cluster_data.copy()
                    cluster_data['cluster'] = cluster_num
                    
                    # Save individual results (ALL genes)
                    cluster_file = os.path.join(complete_results_dir, f'LOUVAIN_CLUSTER_{cluster_num}_COMPLETE_MARKERS.csv')
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
                        log_step(f"(LOUVAIN) Processed {i + 1}/{len(cluster_keys)} clusters")
                        log_memory_usage(f"PROCESSED_{i+1}")
                
            except Exception as e:
                log_step(f"(LOUVAIN) ERROR processing cluster {cluster_key}: {e}")
                continue
        
        # Combine ALL marker genes (complete dataset)
        if all_marker_genes:
            log_step("(LOUVAIN) Combining COMPLETE marker genes dataset")
            log_step("WARNING: HIGH CPU AND MEMORY USAGE")
            combined_complete_markers = pd.concat(all_marker_genes, ignore_index=True)
            complete_markers_file = os.path.join(output_dir,'EXPORTS', 'LOUVAIN_COMPLETE_ALL_MARKERS.csv')
            combined_complete_markers.to_csv(complete_markers_file, index=False)
            
            log_step(f"(LOUVAIN) COMPLETE marker genes dataset saved: {len(combined_complete_markers)} total entries")
            log_step(f"File size: {os.path.getsize(complete_markers_file) / 1e9:.2f} GB")
            
            # Stringent markers
            if 'pvalues_adj' in combined_complete_markers.columns and 'log2fc' in combined_complete_markers.columns:
                stringent_markers = combined_complete_markers[
                    (combined_complete_markers['pvalues_adj'] < 0.001) & 
                    (abs(combined_complete_markers['log2fc']) > 1)
                ]
                stringent_file = os.path.join(marker_genes_dir, 'LOUVAIN_STRINGENT_MARKERS.csv')
                stringent_markers.to_csv(stringent_file, index=False)
                log_step(f"(LOUVAIN) Stringent markers saved: {len(stringent_markers)} genes (p<0.001 and l2fc >1)")
                
                # Moderate markers
                moderate_markers = combined_complete_markers[
                    (combined_complete_markers['pvalues_adj'] < 0.01) & 
                    (abs(combined_complete_markers['log2fc']) > 0.5)
                ]
                moderate_file = os.path.join(marker_genes_dir, 'LOUVAIN_MODERATE_MARKERS.csv')
                moderate_markers.to_csv(moderate_file, index=False)
                log_step(f"(LOUVAIN) Moderate markers saved: {len(moderate_markers)} genes (p<0.01 and l2fc > 0.5)")
                
                # Top markers per cluster
                top_markers_per_cluster = []
                if 'scores' in combined_complete_markers.columns:
                    for cluster in combined_complete_markers['cluster'].unique():
                        cluster_markers = combined_complete_markers[
                            combined_complete_markers['cluster'] == cluster
                        ].nlargest(50, 'scores')
                        top_markers_per_cluster.append(cluster_markers)
                    
                    top_combined = pd.concat(top_markers_per_cluster, ignore_index=True)
                    top_file = os.path.join(marker_genes_dir, 'LOUVAIN_TOP50_PER_CLUSTER.csv')
                    top_combined.to_csv(top_file, index=False)
                    log_step(f"(LOUVAIN) Top 50 per cluster saved: {len(top_combined)} genes")
        
        # Save cluster summaries
        if cluster_summaries:
            summary_df = pd.DataFrame(cluster_summaries)
            summary_file = os.path.join(marker_genes_dir, 'LOUVAIN_CLUSTER_SUMMARY.csv')
            summary_df.to_csv(summary_file, index=False)
            
            # Summary statistics
            detailed_summary_file = os.path.join(marker_genes_dir, 'LOUVAIN_MARKERS_REPORT.txt')
            with open(detailed_summary_file, 'w') as f:
                f.write("(LOUVAIN) MARKER GENE ANALYSIS REPORT\n")
                f.write("="*100 + "\n\n")
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
            processing_file = os.path.join(output_dir, 'LOGS', 'LOUVAIN_CLUSTER_STATISTICS.csv')
            processing_df.to_csv(processing_file, index=False)
            
            total_processing_time = processing_df['processing_time_seconds'].sum()
            avg_processing_time = processing_df['processing_time_seconds'].mean()
            log_step(f"(LOUVAIN) Processing statistics: Total {total_processing_time:.1f}s, Average {avg_processing_time:.2f}s per cluster")
        
        log_memory_usage("ALL_MARKERS_PROCESSED")
        
    except Exception as e:
        log_step(f"(LOUVAIN) ERROR processing marker genes: {e}")
        import traceback
        log_step(f"Traceback: {traceback.format_exc()}")

# STATISTICAL ANALYSIS

log_step("="*100)
log_step("(LOUVAIN) STATISTICAL ANALYSIS")
log_step("="*100)

advanced_stats_dir = os.path.join(output_dir, 'STATISTICAL_ANALYSIS')

# Load the complete marker genes for advanced analysis
complete_markers_file = os.path.join(output_dir,'EXPORTS', 'LOUVAIN_COMPLETE_ALL_MARKERS.csv')
if os.path.exists(complete_markers_file):
    try:
        df_complete_markers = pd.read_csv(complete_markers_file)
        log_step(f"(LOUVAIN) Loaded {len(df_complete_markers)} marker gene records for analysis")
        
        # Distribution analysis
        if 'log2fc' in df_complete_markers.columns and 'pvalues_adj' in df_complete_markers.columns:
            
            # Statistical tests for normality
            log_step("(LOUVAIN) Testing data distributions")
            
            # Test log2fc distribution
            from scipy.stats import shapiro, kstest
            lfc_sample = df_complete_markers['log2fc'].dropna().sample(min(5000, len(df_complete_markers)))
            shapiro_stat, shapiro_p = shapiro(lfc_sample)
            
            # Save statistical test results
            stats_results_file = os.path.join(advanced_stats_dir, 'LOUVAIN_DISTRIBUTION_TESTS.txt')
            with open(stats_results_file, 'w') as f:
                f.write("(LOUVAIN) STATISTICAL DISTRIBUTION ANALYSIS\n")
                f.write("="*100 + "\n\n")
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
        
        # Cluster comparison analysis
        if 'cluster' in df_complete_markers.columns:
            log_step("(LOUVAIN) Performing inter-cluster comparison analysis")
            
            cluster_comparison_file = os.path.join(advanced_stats_dir, 'LOUVAIN_CLUSTER_COMPARISONS.csv')
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
                log_step(f"(LOUVAIN) Cluster comparison statistics saved ({len(cluster_stats)} clusters)")
        
        # Create plots
        if len(df_complete_markers) > 0:
            fig, axes = plt.subplots(2, 3, figsize=(20, 12))
            
            # Log2FC distribution
            if 'log2fc' in df_complete_markers.columns:
                axes[0,0].hist(df_complete_markers['log2fc'], bins=100, alpha=0.7, edgecolor='black')
                axes[0,0].axvline(0, color='red', linestyle='--', alpha=0.7)
                axes[0,0].set_xlabel('Log2 Fold Change')
                axes[0,0].set_ylabel('Frequency')
                axes[0,0].set_title('Distribution of Log2 Fold Changes')
            
            # P-value distribution
            if 'pvalues_adj' in df_complete_markers.columns:
                axes[0,1].hist(-np.log10(df_complete_markers['pvalues_adj'].replace(0, 1e-300)), 
                              bins=100, alpha=0.7, edgecolor='black')
                axes[0,1].set_xlabel('-log10(Adjusted P-value)')
                axes[0,1].set_ylabel('Frequency')
                axes[0,1].set_title('Distribution of -log10(Adj. P-values)')
            
            # Volcano plot (sample)
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
            
            # Cluster-wise statistics
            if 'cluster' in df_complete_markers.columns and len(cluster_stats) > 0:
                cluster_means = [stat['mean_lfc'] for stat in cluster_stats]
                cluster_ids = [stat['cluster'] for stat in cluster_stats]
                axes[1,0].bar(range(len(cluster_means)), cluster_means)
                axes[1,0].set_xlabel('Cluster Index')
                axes[1,0].set_ylabel('Mean Log2FC')
                axes[1,0].set_title('Mean Log2FC by Cluster')
            
            # Significance by cluster
            if 'cluster' in df_complete_markers.columns and len(cluster_stats) > 0:
                sig_counts = [stat.get('sig_01', 0) for stat in cluster_stats]
                axes[1,1].bar(range(len(sig_counts)), sig_counts)
                axes[1,1].set_xlabel('Cluster Index')
                axes[1,1].set_ylabel('Significant Genes (p<0.01)')
                axes[1,1].set_title('Significant Genes by Cluster')
            
            # Log2FC vs Significance
            if 'log2fc' in df_complete_markers.columns and 'pvalues_adj' in df_complete_markers.columns:
                sample_data = df_complete_markers.sample(min(5000, len(df_complete_markers)))
                scatter = axes[1,2].scatter(abs(sample_data['log2fc']), 
                                          -np.log10(sample_data['pvalues_adj'].replace(0, 1e-300)),
                                          alpha=0.6, s=2, c=sample_data.get('scores', 1))
                axes[1,2].set_xlabel('|Log2 Fold Change|')
                axes[1,2].set_ylabel('-log10(Adjusted P-value)')
                axes[1,2].set_title('Effect Size vs Significance')
            
            plt.tight_layout()
            advanced_plots_file = os.path.join(advanced_stats_dir, 'LOUVAIN_STATISTICAL_ANALYSIS.png')
            plt.savefig(advanced_plots_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            log_step("(LOUVAIN) Statistical plots saved")
        
        log_memory_usage("ADVANCED_STATS_COMPLETE")
        
    except Exception as e:
        log_step(f"(LOUVAIN) WARNING: Error in advanced statistical analysis: {e}")

# VISUALIZATION

log_step("="*100)
log_step("(LOUVAIN) VISUALIZATIONS")
log_step("="*100)

viz_marker_dir = os.path.join(output_dir, 'PLOTS')

# Use the filtered dataset for visualizations
top_viz_file = os.path.join(marker_genes_dir, 'LOUVAIN_TOP50_PER_CLUSTER.csv')
if os.path.exists(top_viz_file):
    try:
        log_step("(LOUVAIN) Generating visualizations")
        
        # Standard Stereopy visualizations
        data.plt.marker_genes_text(
            res_key='marker_genes',
            markers_num=10,
            sort_key='scores',
            out_path=os.path.join(viz_marker_dir, 'LOUVAIN_MARKERS.png')
        )
        
        data.plt.marker_genes_scatter(
            res_key='marker_genes', 
            markers_num=10,
            out_path=os.path.join(viz_marker_dir, 'LOUVAIN_MARKERS_SCATTER.png')
        )
        
        log_step("(LOUVAIN) Visualizations completed")
        
    except Exception as e:
        log_step(f"(LOUVAIN) WARNING: Error in Stereopy visualizations: {e}")

# Generate volcano plots for ALL major clusters
try:
    log_step("(LOUVAIN) Generating volcano plots")
    volcano_dir = os.path.join(viz_marker_dir, 'LOUVAIN_VOLCANO_PLOTS_COMPLETE')
    os.makedirs(volcano_dir, exist_ok=True)
    
    if 'louvain' in data.cells:
        all_clusters = sorted(data.cells['louvain'].unique())
        log_step(f"(LOUVAIN) Generating volcano plots for ALL {len(all_clusters)} clusters")
        
        successful_plots = 0
        
        for i, cluster in enumerate(all_clusters):
            try:
                cluster_str = str(cluster)
                group_name = f'{cluster_str}.vs.rest'
                data.plt.marker_gene_volcano(
                    group_name=group_name,
                    res_key='marker_genes',
                    vlines=False,
                    out_path=os.path.join(volcano_dir, f'LOUVAIN_VOLCANO_{cluster_str.zfill(2)}.png')
                )
                
                successful_plots += 1
                
                if (i + 1) % 20 == 0:
                    log_step(f"(LOUVAIN) Generated {successful_plots}/{i + 1} volcano plots successfully")
                    
            except Exception as e:
                log_step(f"(LOUVAIN) WARNING: Error generating volcano plot for cluster {cluster}: {e}")
                try:
                    alt_group_name = f'{int(cluster)}.vs.rest'
                    data.plt.marker_gene_volcano(
                        group_name=alt_group_name,
                        res_key='marker_genes',
                        vlines=False,
                        out_path=os.path.join(volcano_dir, f'VOLCANO_{str(cluster).zfill(2)}_ALT.png')
                    )
                    successful_plots += 1
                    log_step(f"(LOUVAIN) Alternative format worked for cluster {cluster}")
                except Exception as e2:
                    log_step(f"(LOUVAIN) Both formats failed for cluster {cluster}: {e2}")
                continue
        
        log_step(f"(LOUVAIN) Volcano plots completed: {successful_plots}/{len(all_clusters)} successful")
        
        # If no plots were generated, try a different approach
        if successful_plots == 0:
            log_step("(LOUVAIN) No volcano plots generated. Trying alternative approach")
            
            # Check what group names are available in marker_genes result
            if 'marker_genes' in data.tl.result:
                available_groups = [k for k in data.tl.result['marker_genes'].keys() if k.endswith('.vs.rest')]
                log_step(f"(LOUVAIN) Available group names in results: {available_groups[:10]}")
                
                # Try with first few available group names
                for j, group_name in enumerate(available_groups[:10]):
                    try:
                        cluster_id = group_name.replace('.vs.rest', '')
                        data.plt.marker_gene_volcano(
                            group_name=group_name,
                            res_key='marker_genes',
                            vlines=False,
                            out_path=os.path.join(volcano_dir, f'VOLCANO_{cluster_id}_DIRECT.png')
                        )
                        successful_plots += 1
                        log_step(f"(LOUVAIN) Direct approach successful for {group_name}")
                    except Exception as e:
                        log_step(f"(LOUVAIN) Direct approach failed for {group_name}: {e}")
                        continue
            
            log_step(f"(LOUVAIN) Final volcano plots count: {successful_plots}")
    
except Exception as e:
    log_step(f"(LOUVAIN) WARNING: Error in volcano plot generation: {e}")

log_memory_usage("VISUALIZATIONS_COMPLETE")

# CUSTOM GENE ANNOTATION - LOUVAIN
log_step("(LOUVAIN) Applying custom gene interest annotation")
custom_gene_markers = create_custom_gene_markers()

log_step(custom_gene_markers)
custom_annotations_louvain, custom_scores_louvain, annotation_details_louvain = apply_gene_interest_annotation(
    data, custom_gene_markers, cluster_res_key='louvain', threshold=1.2
)

# MARKER GENE FILTERING

log_step("="*100)
log_step("(LOUVAIN) MARKER GENE FILTERING")
log_step("="*100)

filtered_results_dir = os.path.join(marker_genes_dir, 'FILTERED')

try:
    log_step("(LOUVAIN) Statistical filtering of marker genes")
    
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
            
            log_step(f"(LOUVAIN) Applied {config['name']} filtering (FC>{config['min_fold_change']})")
            
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
                        
                        filter_file = os.path.join(filtered_results_dir, f'LOUVAIN_FILTERED_MARKERS_{config["name"].upper()}.csv')
                        combined_filtered.to_csv(filter_file, index=False)
                        
                        log_step(f"(LOUVAIN) Saved {len(combined_filtered)} {config['name']} filtered markers")
            
        except Exception as e:
            log_step(f"(LOUVAIN) WARNING: Error in {config['name']} filtering: {e}")
            continue
    
    log_memory_usage("FILTERING_COMPLETE")
    
except Exception as e:
    log_step(f"(LOUVAIN) WARNING: Error in filtering: {e}")

# LEIDEN CLUSTERING
log_step("="*100)
log_step("LEIDEN CLUSTERING")
log_step("="*100)

try:
    data.tl.leiden(
        neighbors_res_key='neighbors',
        res_key='leiden'
    )
    log_step("(LEIDEN) Clustering completed")
    
    # Plot Leiden clusters
    data.plt.cluster_scatter(
        res_key='leiden',
        out_path=os.path.join(output_dir,'CLUSTERING', 'LEIDEN_CLUSTERS.png')
    )
    log_step("(SPATIAL_LEIDEN) Spatial scatter plot saved")
    
    # Plot UMAP with Leiden clusters
    data.plt.umap(
        res_key='umap',
        cluster_key='leiden',
        out_path=os.path.join(output_dir,'CLUSTERING', 'LEIDEN_UMAP_CLUSTERS.png')
    )
    log_step("(LEIDEN) UMAP plot saved")

    # MARKER GENE ANALYSIS
    ## This step process all genes and all clusters (requires HPC)
    log_step("="*100)
    log_step("(LEIDEN) MARKER GENE ANALYSIS")
    log_step("="*100)

    marker_start_time = datetime.now()
    log_step("(LEIDEN) Finding marker genes using t-test method")

    # Generate marker genes for the spatial leiden clusters
    data.tl.find_marker_genes(
            cluster_res_key='leiden',
            method='t_test',
            use_highly_genes=False,
            use_raw=True,
            res_key='ldn_marker_genes'
            )

    marker_end_time = datetime.now()
    marker_duration = (marker_end_time - marker_start_time).total_seconds()
    log_step(f"(LEIDEN) Marker gene analysis completed in {marker_duration:.1f} seconds")
    log_memory_usage("MARKER_GENES_FOUND")
    
    # Plot marker genes
    data.plt.marker_genes_scatter(res_key='ldn_marker_genes', markers_num=5, out_path=os.path.join(output_dir, 'PLOTS', 'LEIDEN_MARKERS_SCATTER.png'))

    # Generate individual Leiden cluster plots
    log_step("(LEIDEN) Generating individual Leiden cluster plots")
    leiden_unique_clusters = sorted(data.cells['leiden'].unique())
    leiden_n_clusters = len(leiden_unique_clusters)
    log_step(f"(LEIDEN) Found {leiden_n_clusters} clusters.")
    save_progress_checkpoint(data, output_dir, 'LEIDEN')
    log_memory_usage("LEIDEN_COMPLETE")
    
except Exception as e:
    log_step(f"(LEIDEN) ERROR in Leiden clustering: {e}")
    import traceback
    log_step(f"(LEIDEN) Full traceback: {traceback.format_exc()}")
    sys.exit(1)

if 'ldn_marker_genes' in data.tl.result:
    try:
        marker_result = data.tl.result['ldn_marker_genes']
        log_step(f"(LEIDEN) Processing marker genes results with {len(marker_result)} components")
        
        # Get all cluster comparison keys
        cluster_keys = [k for k in marker_result.keys() if k.endswith('.vs.rest')]
        log_step(f"(LEIDEN) Processing ALL {len(cluster_keys)} cluster comparisons")
        
        # Process ALL clusters with ALL genes
        all_marker_genes = []
        cluster_summaries = []
        processing_stats = []
        
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
                    cluster_file = os.path.join(complete_results_dir, f'LEIDEN_CLUSTER_{cluster_num}_COMPLETE_MARKERS.csv')
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
                        log_step(f"(LEIDEN) Processed {i + 1}/{len(cluster_keys)} clusters")
                        log_memory_usage(f"PROCESSED_{i+1}")
                
            except Exception as e:
                log_step(f"(LEIDEN) ERROR processing cluster {cluster_key}: {e}")
                continue
        
        # Combine ALL marker genes (complete dataset)
        if all_marker_genes:
            log_step("(LEIDEN) Combining COMPLETE marker genes dataset")
            log_step("(LEIDEN) WARNING: HIGH CPU AND MEMORY USAGE")
            combined_complete_markers = pd.concat(all_marker_genes, ignore_index=True)
            complete_markers_file = os.path.join(output_dir,'EXPORTS', 'LEIDEN_COMPLETE_ALL_MARKERS.csv')
            combined_complete_markers.to_csv(complete_markers_file, index=False)
            
            log_step(f"(LEIDEN) Marker genes dataset saved: {len(combined_complete_markers)} total entries")
            log_step(f"(LEIDEN) File size: {os.path.getsize(complete_markers_file) / 1e9:.2f} GB")
            
            # Create filtered versions
            log_step("(LEIDEN) Creating filtered versions for practical analysis")
            
            # Stringent markers
            if 'pvalues_adj' in combined_complete_markers.columns and 'log2fc' in combined_complete_markers.columns:
                stringent_markers = combined_complete_markers[
                    (combined_complete_markers['pvalues_adj'] < 0.001) & 
                    (abs(combined_complete_markers['log2fc']) > 1)
                ]
                stringent_file = os.path.join(marker_genes_dir, 'LEIDEN_STRINGENT_MARKERS.csv')
                stringent_markers.to_csv(stringent_file, index=False)
                log_step(f"(LEIDEN) Stringent markers saved: {len(stringent_markers)} genes (p<0.001 and l2fc >1)")
                
                # Moderate markers
                moderate_markers = combined_complete_markers[
                    (combined_complete_markers['pvalues_adj'] < 0.01) & 
                    (abs(combined_complete_markers['log2fc']) > 0.5)
                ]
                moderate_file = os.path.join(marker_genes_dir, 'LEIDEN_MODERATE_MARKERS.csv')
                moderate_markers.to_csv(moderate_file, index=False)
                log_step(f"(LEIDEN) Moderate markers saved: {len(moderate_markers)} genes (p<0.01 and l2fc >0.5)")
                
                # Top markers per cluster
                top_markers_per_cluster = []
                if 'scores' in combined_complete_markers.columns:
                    for cluster in combined_complete_markers['cluster'].unique():
                        cluster_markers = combined_complete_markers[
                            combined_complete_markers['cluster'] == cluster
                        ].nlargest(50, 'scores')
                        top_markers_per_cluster.append(cluster_markers)
                    
                    top_combined = pd.concat(top_markers_per_cluster, ignore_index=True)
                    top_file = os.path.join(marker_genes_dir, 'LEIDEN_TOP50_PER_CLUSTER.csv')
                    top_combined.to_csv(top_file, index=False)
                    log_step(f"(LEIDEN) Top 50 per cluster saved: {len(top_combined)} genes")
        
        # Save cluster summaries
        if cluster_summaries:
            summary_df = pd.DataFrame(cluster_summaries)
            summary_file = os.path.join(marker_genes_dir, 'LEIDEN_CLUSTER_SUMMARY.csv')
            summary_df.to_csv(summary_file, index=False)
            
            # Summary statistics
            detailed_summary_file = os.path.join(marker_genes_dir, 'LEIDEN_MARKERS_REPORT.txt')
            with open(detailed_summary_file, 'w') as f:
                f.write("(LEIDEN) MARKER GENE ANALYSIS REPORT\n")
                f.write("="*100 + "\n\n")
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
                
                f.write("CLUSTER STATISTICS:\n")
                f.write(str(summary_df) + "\n")
        
        # Save processing statistics
        if processing_stats:
            processing_df = pd.DataFrame(processing_stats)
            processing_file = os.path.join(output_dir, 'LOGS', 'LEIDEN_CLUSTER_STATISTICS.csv')
            processing_df.to_csv(processing_file, index=False)
            
            total_processing_time = processing_df['processing_time_seconds'].sum()
            avg_processing_time = processing_df['processing_time_seconds'].mean()
            log_step(f"(LEIDEN) Processing statistics: Total {total_processing_time:.1f}s, Average {avg_processing_time:.2f}s per cluster")
        
        log_memory_usage("ALL_MARKERS_PROCESSED")
        
    except Exception as e:
        log_step(f"(LEIDEN) ERROR processing marker genes: {e}")
        import traceback
        log_step(f"(LEIDEN) Traceback: {traceback.format_exc()}")

# STATISTICAL ANALYSIS

log_step("="*100)
log_step("(LEIDEN) STATISTICAL ANALYSIS")
log_step("="*100)

# Load the complete marker genes for advanced analysis
complete_markers_file = os.path.join(output_dir,'EXPORTS', 'LEIDEN_COMPLETE_ALL_MARKERS.csv')
if os.path.exists(complete_markers_file):
    try:
        df_complete_markers = pd.read_csv(complete_markers_file)
        log_step(f"(LEIDEN) Loaded {len(df_complete_markers)} marker gene records for analysis")
        
        # Distribution analysis
        if 'log2fc' in df_complete_markers.columns and 'pvalues_adj' in df_complete_markers.columns:
            
            # Test log2fc distribution
            from scipy.stats import shapiro, kstest
            lfc_sample = df_complete_markers['log2fc'].dropna().sample(min(5000, len(df_complete_markers)))
            shapiro_stat, shapiro_p = shapiro(lfc_sample)
            
            # Save statistical test results
            stats_results_file = os.path.join(advanced_stats_dir, 'LEIDEN_DISTRIBUTION_TESTS.txt')
            with open(stats_results_file, 'w') as f:
                f.write("(LEIDEN) STATISTICAL DISTRIBUTION ANALYSIS\n")
                f.write("="*100 + "\n\n")
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
        
        # Cluster comparison analysis
        if 'cluster' in df_complete_markers.columns:
            log_step("(LEIDEN) Performing inter-cluster comparison analysis")
            
            cluster_comparison_file = os.path.join(advanced_stats_dir, 'LEIDEN_CLUSTER_COMPARISONS.csv')
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
                log_step(f"(LEIDEN) Cluster comparison statistics saved ({len(cluster_stats)} clusters)")
        
        # Create plots
        if len(df_complete_markers) > 0:
            fig, axes = plt.subplots(2, 3, figsize=(20, 12))
            
            # Log2FC distribution
            if 'log2fc' in df_complete_markers.columns:
                axes[0,0].hist(df_complete_markers['log2fc'], bins=100, alpha=0.7, edgecolor='black')
                axes[0,0].axvline(0, color='red', linestyle='--', alpha=0.7)
                axes[0,0].set_xlabel('Log2 Fold Change')
                axes[0,0].set_ylabel('Frequency')
                axes[0,0].set_title('Distribution of Log2 Fold Changes')
            
            # P-value distribution
            if 'pvalues_adj' in df_complete_markers.columns:
                axes[0,1].hist(-np.log10(df_complete_markers['pvalues_adj'].replace(0, 1e-300)), 
                              bins=100, alpha=0.7, edgecolor='black')
                axes[0,1].set_xlabel('-log10(Adjusted P-value)')
                axes[0,1].set_ylabel('Frequency')
                axes[0,1].set_title('Distribution of -log10(Adj. P-values)')
            
            # Volcano plot (sample)
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
            
            # Cluster-wise statistics
            if 'cluster' in df_complete_markers.columns and len(cluster_stats) > 0:
                cluster_means = [stat['mean_lfc'] for stat in cluster_stats]
                cluster_ids = [stat['cluster'] for stat in cluster_stats]
                axes[1,0].bar(range(len(cluster_means)), cluster_means)
                axes[1,0].set_xlabel('Cluster Index')
                axes[1,0].set_ylabel('Mean Log2FC')
                axes[1,0].set_title('Mean Log2FC by Cluster')
            
            # Significance by cluster
            if 'cluster' in df_complete_markers.columns and len(cluster_stats) > 0:
                sig_counts = [stat.get('sig_01', 0) for stat in cluster_stats]
                axes[1,1].bar(range(len(sig_counts)), sig_counts)
                axes[1,1].set_xlabel('Cluster Index')
                axes[1,1].set_ylabel('Significant Genes (p<0.01)')
                axes[1,1].set_title('Significant Genes by Cluster')
            
            # Log2FC vs Significance
            if 'log2fc' in df_complete_markers.columns and 'pvalues_adj' in df_complete_markers.columns:
                sample_data = df_complete_markers.sample(min(5000, len(df_complete_markers)))
                scatter = axes[1,2].scatter(abs(sample_data['log2fc']), 
                                          -np.log10(sample_data['pvalues_adj'].replace(0, 1e-300)),
                                          alpha=0.6, s=2, c=sample_data.get('scores', 1))
                axes[1,2].set_xlabel('|Log2 Fold Change|')
                axes[1,2].set_ylabel('-log10(Adjusted P-value)')
                axes[1,2].set_title('Effect Size vs Significance')
            
            plt.tight_layout()
            advanced_plots_file = os.path.join(advanced_stats_dir, 'LEIDEN_STATISTICAL_ANALYSIS.png')
            plt.savefig(advanced_plots_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            log_step("(LEIDEN) Statistical plots saved")
        
        log_memory_usage("ADVANCED_STATS_COMPLETE")
        
    except Exception as e:
        log_step(f"(LEIDEN) WARNING: Error in statistical analysis: {e}")

# VISUALIZATION

log_step("="*100)
log_step("(LEIDEN) VISUALIZATIONS")
log_step("="*100)

# Use the filtered dataset for visualizations
top_viz_file = os.path.join(marker_genes_dir, 'LEIDEN_TOP50_PER_CLUSTER.csv')
if os.path.exists(top_viz_file):
    try:
        # Standard Stereopy visualizations
        data.plt.marker_genes_text(
            res_key='ldn_marker_genes',
            markers_num=10,
            sort_key='scores',
            out_path=os.path.join(viz_marker_dir, 'LEIDEN_MARKERS.png')
        )
        
        data.plt.marker_genes_scatter(
            res_key='ldn_marker_genes', 
            markers_num=10,
            out_path=os.path.join(viz_marker_dir, 'LEIDEN_MARKERS_SCATTER.png')
        )
        
        log_step("(LEIDEN) Visualizations completed")
        
    except Exception as e:
        log_step(f"(LEIDEN) WARNING: Error in visualizations: {e}")

# Generate volcano plots for ALL major clusters
try:
    log_step("(LEIDEN) Generating volcano plots")
    volcano_dir = os.path.join(viz_marker_dir, 'LEIDEN_VOLCANO_PLOTS_COMPLETE')
    os.makedirs(volcano_dir, exist_ok=True)
    
    if 'leiden' in data.cells:
        all_clusters = sorted(data.cells['leiden'].unique())
        log_step(f"(LEIDEN) Generating volcano plots for ALL {len(all_clusters)} clusters")
        
        successful_plots = 0
        
        for i, cluster in enumerate(all_clusters):
            try:
                cluster_str = str(cluster)
                group_name = f'{cluster_str}.vs.rest'
                data.plt.marker_gene_volcano(
                    group_name=group_name,
                    res_key='ldn_marker_genes',
                    vlines=False,
                    out_path=os.path.join(volcano_dir, f'LEIDEN_VOLCANO_CLUSTER_{cluster_str.zfill(2)}.png')
                )
                
                successful_plots += 1
                
                if (i + 1) % 20 == 0:
                    log_step(f"(LEIDEN) Generated {successful_plots}/{i + 1} volcano plots successfully")
                    
            except Exception as e:
                log_step(f"(LEIDEN) WARNING: Error generating volcano plot for cluster {cluster}: {e}")
                try:
                    alt_group_name = f'{int(cluster)}.vs.rest'
                    data.plt.marker_gene_volcano(
                        group_name=alt_group_name,
                        res_key='ldn_marker_genes',
                        vlines=False,
                        out_path=os.path.join(volcano_dir, f'LEIDEN_VOLCANO_CLUSTER_{str(cluster).zfill(2)}_ALT.png')
                    )
                    successful_plots += 1
                    log_step(f"(LEIDEN) Alternative format worked for cluster {cluster}")
                except Exception as e2:
                    log_step(f"(LEIDEN) Both formats failed for cluster {cluster}: {e2}")
                continue
        
        log_step(f"(LEIDEN) Volcano plots completed: {successful_plots}/{len(all_clusters)} successful")
        
        # If no plots were generated, try a different approach
        if successful_plots == 0:
            log_step("(LEIDEN) No volcano plots generated. Trying alternative approach")
            
            # Check what group names are available in marker_genes result
            if 'ldn_marker_genes' in data.tl.result:
                available_groups = [k for k in data.tl.result['ldn_marker_genes'].keys() if k.endswith('.vs.rest')]
                log_step(f"(LEIDEN) Available group names in results: {available_groups[:10]}")
                
                # Try with first few available group names
                for j, group_name in enumerate(available_groups[:10]):
                    try:
                        cluster_id = group_name.replace('.vs.rest', '')
                        data.plt.marker_gene_volcano(
                            group_name=group_name,
                            res_key='ldn_marker_genes',
                            vlines=False,
                            out_path=os.path.join(volcano_dir, f'LEIDEN_VOLCANO_CLUSTER_{cluster_id}_DIRECT.png')
                        )
                        successful_plots += 1
                        log_step(f"(LEIDEN) Direct approach successful for {group_name}")
                    except Exception as e:
                        log_step(f"(LEIDEN) Direct approach failed for {group_name}: {e}")
                        continue
            
            log_step(f"(LEIDEN) {successful_plots} volcano plots generated")
    
except Exception as e:
    log_step(f"(LEIDEN) WARNING: Error in volcano plot generation: {e}")

log_memory_usage("VISUALIZATIONS_COMPLETE")

# CUSTOM GENE ANNOTATION
log_step("(LEIDEN) Applying custom gene interest annotation")

custom_gene_markers = create_custom_gene_markers()
log_step(custom_gene_markers)
custom_annotations, custom_scores, annotation_details = apply_gene_interest_annotation(
    data, custom_gene_markers, cluster_res_key='leiden', threshold=1.2
)

# Create detailed report for custom gene annotation
custom_report_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'LEIDEN_INTEREST_ANALYSIS_REPORT.txt')
with open(custom_report_file, 'w') as f:
    f.write("(LEIDEN) CUSTOM GENES OF INTEREST ANNOTATION REPORT\n")
    f.write("="*100 + "\n\n")
    
    f.write("GENES OF INTEREST DEFINED:\n")
    for category, genes in custom_gene_markers.items():
        f.write(f"\n{category.upper()}:\n")
        if genes:
            for gene in genes:
                f.write(f"  - {gene}\n")
        else:
            f.write("  - (No genes defined yet)\n")
    
    f.write(f"\nANNOTATION RESULTS (threshold: 1.2):\n")
    f.write("=" * 100 + "\n")
    
    for cluster in sorted(annotation_details.keys(), key=lambda x: int(x)):
        details = annotation_details[cluster]
        f.write(f"\nCLUSTER {cluster}:\n")
        f.write(f"  Annotation: {custom_annotations[cluster]}\n")
        f.write(f"  Cell count: {details['cell_count']}\n")
        f.write(f"  Best category: {details['best_category'] or 'None'}\n")
        f.write(f"  Best score: {details['best_score']:.3f}\n")

       	if details['all_scores']:
            f.write("  Category scores:\n")
            for category, score in details['all_scores'].items():
                f.write(f"    {category}: {score:.3f}\n")

log_step(f"(LEIDEN) Custom gene interest annotation completed for {len(custom_annotations)} clusters")

# TOP 50 MARKER GENES
marker_genes_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'LEIDEN_TOP50_PER_CLUSTER.csv')
biological_annotations = {}

if os.path.exists(marker_genes_file):
    try:
        marker_df = pd.read_csv(marker_genes_file)
        log_step(f"(LEIDEN) Loaded {len(marker_df)} marker gene records")
        
        # Create biological annotations based on marker genes
        biological_annotations = create_biological_annotation_dict(marker_df, n_top_genes=5)
        log_step(f"(LEIDEN) Created annotations for {len(biological_annotations)} clusters")
        
        # Save biological annotation mapping
        bio_anno_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'LEIDEN_INTEREST_ANALYSIS.json')
        with open(bio_anno_file, 'w') as f:
            json.dump(biological_annotations, f, indent=2)
        
        # Create marker-based annotation report
        marker_report_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'LEIDEN_MARKERS_REPORT.txt')
        with open(marker_report_file, 'w') as f:
            f.write("MARKER-BASED BIOLOGICAL ANNOTATION REPORT\n")
            f.write("="*100 + "\n\n")
            
            for cluster in sorted(biological_annotations.keys(), key=lambda x: int(x)):
                f.write(f"CLUSTER {cluster}:\n")
                f.write(f"Annotation: {biological_annotations[cluster]}\n")
                
                # Get top markers for this cluster
                cluster_markers = marker_df[marker_df['cluster'] == int(cluster)]
                if len(cluster_markers) > 0:
                    if 'scores' in cluster_markers.columns:
                        top_markers = cluster_markers.nlargest(10, 'scores')
                    else:
                        top_markers = cluster_markers.head(10)
                    
                    f.write("Top marker genes:\n")
                    for _, row in top_markers.iterrows():
                        gene = row['genes'] if 'genes' in row else 'Unknown'
                        score = row['scores'] if 'scores' in row else 'N/A'
                        f.write(f"  - {gene}: {score}\n")
                
                f.write("\n")
        
    except Exception as e:
        log_step(f"(LEIDEN) WARNING: Error loading marker genes for annotation: {e}")

# MARKER GENE FILTERING

log_step("="*100)
log_step("(LEIDEN) MARKER GENE FILTERING")
log_step("="*100)

try:
    # Apply multiple filtering criteria
    filter_configs = [
        {'min_fold_change': 0.5, 'min_in_group_fraction': 0.2, 'max_out_group_fraction': 0.6, 'name': 'lenient'},
        {'min_fold_change': 1.0, 'min_in_group_fraction': 0.25, 'max_out_group_fraction': 0.5, 'name': 'moderate'},
        {'min_fold_change': 1.5, 'min_in_group_fraction': 0.3, 'max_out_group_fraction': 0.4, 'name': 'stringent'},
        {'min_fold_change': 2.0, 'min_in_group_fraction': 0.4, 'max_out_group_fraction': 0.3, 'name': 'very_stringent'}
    ]
    
    for config in filter_configs:
        try:
            filter_key = f"LEIDEN_marker_genes_filtered_{config['name']}"
            
            data.tl.filter_marker_genes(
                marker_genes_res_key='ldn_marker_genes',
                min_fold_change=config['min_fold_change'],
                min_in_group_fraction=config['min_in_group_fraction'],
                max_out_group_fraction=config['max_out_group_fraction'],
                res_key=filter_key
            )
            
            log_step(f"(LEIDEN) Applied {config['name']} filtering (FC>{config['min_fold_change']})")
            
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
                        
                        filter_file = os.path.join(filtered_results_dir, f'LEIDEN_FILTERED_MARKERS_{config["name"].upper()}.csv')
                        combined_filtered.to_csv(filter_file, index=False)
                        
                        log_step(f"(LEIDEN) Saved {len(combined_filtered)} {config['name']} filtered markers")
            
        except Exception as e:
            log_step(f"(LEIDEN) WARNING: Error in {config['name']} filtering: {e}")
            continue
    
    log_memory_usage("FILTERING_COMPLETE")
    
except Exception as e:
    log_step(f"(LEIDEN) WARNING: Error in advanced filtering: {e}")

# SPATIAL LEIDEN CLUSTERING
log_step("="*100)
log_step("SPATIAL LEIDEN CLUSTERING")
log_step("="*100)

try:
    data.tl.leiden(
        neighbors_res_key='spatial_neighbors',
        res_key='spatial_leiden'
    )
    log_step("(SPATIAL LEIDEN) Clustering completed")
    
    # Plot Leiden clusters
    data.plt.cluster_scatter(
        res_key='spatial_leiden',
        out_path=os.path.join(output_dir,'CLUSTERING', 'SPATIAL_LEIDEN_CLUSTERS.png')
    )
    log_step("(SPATIAL_LEIDEN) Spatial scatter plot saved")
    
    # Plot UMAP with Leiden clusters
    data.plt.umap(
        res_key='umap',
        cluster_key='spatial_leiden',
        out_path=os.path.join(output_dir,'CLUSTERING', 'SPATIAL_LEIDEN_UMAP_CLUSTERS.png')
    )
    log_step("(SPATIAL LEIDEN) UMAP plot saved")

    # MARKER GENE ANALYSIS
    ## This step process all genes and all clusters (requires HPC)
    log_step("="*100)
    log_step("(SPATIAL LEIDEN) MARKER GENE ANALYSIS")
    log_step("="*100)

    marker_start_time = datetime.now()
    log_step("(SPATIAL LEIDEN) Finding marker genes using t-test method")

    # Generate marker genes for the spatial leiden clusters
    data.tl.find_marker_genes(
            cluster_res_key='spatial_leiden',
            method='t_test',
            use_highly_genes=False,
            use_raw=True,
            res_key='marker_genes_spldn'
            )

    marker_end_time = datetime.now()
    marker_duration = (marker_end_time - marker_start_time).total_seconds()
    log_step(f"(SPATIAL LEIDEN) Marker gene analysis completed in {marker_duration:.1f} seconds")
    log_memory_usage("MARKER_GENES_FOUND")
    
    # Plot marker genes
    data.plt.marker_genes_scatter(res_key='marker_genes_spldn', markers_num=5, out_path=os.path.join(output_dir, 'PLOTS', 'SPATIAL_LEIDEN_MARKERS_SCATTER.png'))

    # Generate individual Spatial Leiden cluster plots
    log_step("(SPATIAL LEIDEN) Generating individual Leiden cluster plots")
    spatial_leiden_unique_clusters = sorted(data.cells['spatial_leiden'].unique())
    spatial_leiden_n_clusters = len(spatial_leiden_unique_clusters)
    log_step(f"(SPATIAL LEIDEN) Found {spatial_leiden_n_clusters} clusters.")
    
except Exception as e:
    log_step(f"(SPATIAL LEIDEN) ERROR in Spatial Leiden clustering: {e}")
    import traceback
    log_step(f"(SPATIAL LEIDEN) Full traceback: {traceback.format_exc()}")

if 'marker_genes_spldn' in data.tl.result:
    try:
        marker_result = data.tl.result['marker_genes_spldn']
        log_step(f"(SPATIAL LEIDEN) Processing marker genes results with {len(marker_result)} components")
        
        # Get all cluster comparison keys
        cluster_keys = [k for k in marker_result.keys() if k.endswith('.vs.rest')]
        log_step(f"(SPATIAL LEIDEN) Processing ALL {len(cluster_keys)} cluster comparisons")
        
        # Process ALL clusters with ALL genes
        all_marker_genes = []
        cluster_summaries = []
        processing_stats = []
        
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
                    cluster_file = os.path.join(complete_results_dir, f'SPATIAL_LEIDEN_CLUSTER_{cluster_num}_COMPLETE_MARKERS.csv')
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
                        log_step(f"(SPATIAL LEIDEN) Processed {i + 1}/{len(cluster_keys)} clusters")
                        log_memory_usage(f"PROCESSED_{i+1}")
                
            except Exception as e:
                log_step(f"(SPATIAL LEIDEN) ERROR processing cluster {cluster_key}: {e}")
                continue
        
        # Combine ALL marker genes (complete dataset)
        if all_marker_genes:
            log_step("(SPATIAL LEIDEN) Combining COMPLETE marker genes dataset")
            log_step("(SPATIAL LEIDEN) WARNING: HIGH CPU AND MEMORY USAGE")
            combined_complete_markers = pd.concat(all_marker_genes, ignore_index=True)
            complete_markers_file = os.path.join(output_dir,'EXPORTS', 'SPATIAL_LEIDEN_COMPLETE_ALL_MARKERS.csv')
            combined_complete_markers.to_csv(complete_markers_file, index=False)
            
            log_step(f"(SPATIAL LEIDEN) Marker genes dataset saved: {len(combined_complete_markers)} total entries")
            log_step(f"(SPATIAL LEIDEN) File size: {os.path.getsize(complete_markers_file) / 1e9:.2f} GB")
            
            # Create filtered versions
            log_step("(SPATIAL LEIDEN) Creating filtered versions for practical analysis")
            
            # Stringent markers
            if 'pvalues_adj' in combined_complete_markers.columns and 'log2fc' in combined_complete_markers.columns:
                stringent_markers = combined_complete_markers[
                    (combined_complete_markers['pvalues_adj'] < 0.001) & 
                    (abs(combined_complete_markers['log2fc']) > 1)
                ]
                stringent_file = os.path.join(marker_genes_dir, 'SPATIAL_LEIDEN_STRINGENT_MARKERS.csv')
                stringent_markers.to_csv(stringent_file, index=False)
                log_step(f"(SPATIAL LEIDEN) Stringent markers saved: {len(stringent_markers)} genes (p<0.001 and l2fc >1)")
                
                # Moderate markers
                moderate_markers = combined_complete_markers[
                    (combined_complete_markers['pvalues_adj'] < 0.01) & 
                    (abs(combined_complete_markers['log2fc']) > 0.5)
                ]
                moderate_file = os.path.join(marker_genes_dir, 'SPATIAL_LEIDEN_MODERATE_MARKERS.csv')
                moderate_markers.to_csv(moderate_file, index=False)
                log_step(f"(SPATIAL LEIDEN) Moderate markers saved: {len(moderate_markers)} genes (p<0.01 and l2fc >0.5)")
                
                # Top markers per cluster
                top_markers_per_cluster = []
                if 'scores' in combined_complete_markers.columns:
                    for cluster in combined_complete_markers['cluster'].unique():
                        cluster_markers = combined_complete_markers[
                            combined_complete_markers['cluster'] == cluster
                        ].nlargest(50, 'scores')
                        top_markers_per_cluster.append(cluster_markers)
                    
                    top_combined = pd.concat(top_markers_per_cluster, ignore_index=True)
                    top_file = os.path.join(marker_genes_dir, 'SPATIAL_LEIDEN_TOP50_PER_CLUSTER.csv')
                    top_combined.to_csv(top_file, index=False)
                    log_step(f"(SPATIAL LEIDEN) Top 50 per cluster saved: {len(top_combined)} genes")
        
        # Save cluster summaries
        if cluster_summaries:
            summary_df = pd.DataFrame(cluster_summaries)
            summary_file = os.path.join(marker_genes_dir, 'SPATIAL_LEIDEN_CLUSTER_SUMMARY.csv')
            summary_df.to_csv(summary_file, index=False)
            
            # Summary statistics
            detailed_summary_file = os.path.join(marker_genes_dir, 'SPATIAL_LEIDEN_MARKERS_REPORT.txt')
            with open(detailed_summary_file, 'w') as f:
                f.write("(SPATIAL LEIDEN) MARKER GENE ANALYSIS REPORT\n")
                f.write("="*100 + "\n\n")
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
                
                f.write("CLUSTER STATISTICS:\n")
                f.write(str(summary_df) + "\n")
        
        # Save processing statistics
        if processing_stats:
            processing_df = pd.DataFrame(processing_stats)
            processing_file = os.path.join(output_dir, 'LOGS', 'SPATIAL_LEIDEN_CLUSTER_STATISTICS.csv')
            processing_df.to_csv(processing_file, index=False)
            
            total_processing_time = processing_df['processing_time_seconds'].sum()
            avg_processing_time = processing_df['processing_time_seconds'].mean()
            log_step(f"(SPATIAL LEIDEN) Processing statistics: Total {total_processing_time:.1f}s, Average {avg_processing_time:.2f}s per cluster")
        
        log_memory_usage("ALL_MARKERS_PROCESSED")
        
    except Exception as e:
        log_step(f"(SPATIAL LEIDEN) ERROR processing marker genes: {e}")
        import traceback
        log_step(f"(SPATIAL LEIDEN) Traceback: {traceback.format_exc()}")

# STATISTICAL ANALYSIS

log_step("="*100)
log_step("(SPATIAL LEIDEN) STATISTICAL ANALYSIS")
log_step("="*100)

# Load the complete marker genes for advanced analysis
complete_markers_file = os.path.join(output_dir,'EXPORTS', 'SPATIAL_LEIDEN_COMPLETE_ALL_MARKERS.csv')
if os.path.exists(complete_markers_file):
    try:
        df_complete_markers = pd.read_csv(complete_markers_file)
        log_step(f"(SPATIAL LEIDEN) Loaded {len(df_complete_markers)} marker gene records for analysis")
        
        # Distribution analysis
        if 'log2fc' in df_complete_markers.columns and 'pvalues_adj' in df_complete_markers.columns:
            
            # Test log2fc distribution
            from scipy.stats import shapiro, kstest
            lfc_sample = df_complete_markers['log2fc'].dropna().sample(min(5000, len(df_complete_markers)))
            shapiro_stat, shapiro_p = shapiro(lfc_sample)
            
            # Save statistical test results
            stats_results_file = os.path.join(advanced_stats_dir, 'SPATIAL_LEIDEN_DISTRIBUTION_TESTS.txt')
            with open(stats_results_file, 'w') as f:
                f.write("(SPATIAL LEIDEN) STATISTICAL DISTRIBUTION ANALYSIS\n")
                f.write("="*100 + "\n\n")
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
        
        # Cluster comparison analysis
        if 'cluster' in df_complete_markers.columns:
            log_step("(SPATIAL LEIDEN) Performing inter-cluster comparison analysis")
            
            cluster_comparison_file = os.path.join(advanced_stats_dir, 'SPATIAL_LEIDEN_CLUSTER_COMPARISONS.csv')
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
                log_step(f"(SPATIAL LEIDEN) Cluster comparison statistics saved ({len(cluster_stats)} clusters)")
        
        # Create plots
        if len(df_complete_markers) > 0:
            fig, axes = plt.subplots(2, 3, figsize=(20, 12))
            
            # Log2FC distribution
            if 'log2fc' in df_complete_markers.columns:
                axes[0,0].hist(df_complete_markers['log2fc'], bins=100, alpha=0.7, edgecolor='black')
                axes[0,0].axvline(0, color='red', linestyle='--', alpha=0.7)
                axes[0,0].set_xlabel('Log2 Fold Change')
                axes[0,0].set_ylabel('Frequency')
                axes[0,0].set_title('Distribution of Log2 Fold Changes')
            
            # P-value distribution
            if 'pvalues_adj' in df_complete_markers.columns:
                axes[0,1].hist(-np.log10(df_complete_markers['pvalues_adj'].replace(0, 1e-300)), 
                              bins=100, alpha=0.7, edgecolor='black')
                axes[0,1].set_xlabel('-log10(Adjusted P-value)')
                axes[0,1].set_ylabel('Frequency')
                axes[0,1].set_title('Distribution of -log10(Adj. P-values)')
            
            # Volcano plot (sample)
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
            
            # Cluster-wise statistics
            if 'cluster' in df_complete_markers.columns and len(cluster_stats) > 0:
                cluster_means = [stat['mean_lfc'] for stat in cluster_stats]
                cluster_ids = [stat['cluster'] for stat in cluster_stats]
                axes[1,0].bar(range(len(cluster_means)), cluster_means)
                axes[1,0].set_xlabel('Cluster Index')
                axes[1,0].set_ylabel('Mean Log2FC')
                axes[1,0].set_title('Mean Log2FC by Cluster')
            
            # Significance by cluster
            if 'cluster' in df_complete_markers.columns and len(cluster_stats) > 0:
                sig_counts = [stat.get('sig_01', 0) for stat in cluster_stats]
                axes[1,1].bar(range(len(sig_counts)), sig_counts)
                axes[1,1].set_xlabel('Cluster Index')
                axes[1,1].set_ylabel('Significant Genes (p<0.01)')
                axes[1,1].set_title('Significant Genes by Cluster')
            
            # Log2FC vs Significance
            if 'log2fc' in df_complete_markers.columns and 'pvalues_adj' in df_complete_markers.columns:
                sample_data = df_complete_markers.sample(min(5000, len(df_complete_markers)))
                scatter = axes[1,2].scatter(abs(sample_data['log2fc']), 
                                          -np.log10(sample_data['pvalues_adj'].replace(0, 1e-300)),
                                          alpha=0.6, s=2, c=sample_data.get('scores', 1))
                axes[1,2].set_xlabel('|Log2 Fold Change|')
                axes[1,2].set_ylabel('-log10(Adjusted P-value)')
                axes[1,2].set_title('Effect Size vs Significance')
            
            plt.tight_layout()
            advanced_plots_file = os.path.join(advanced_stats_dir, 'SPATIAL_LEIDEN_STATISTICAL_ANALYSIS.png')
            plt.savefig(advanced_plots_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            log_step("(SPATIAL LEIDEN) Statistical plots saved")
        
        log_memory_usage("ADVANCED_STATS_COMPLETE")
        
    except Exception as e:
        log_step(f"(SPATIAL LEIDEN) WARNING: Error in statistical analysis: {e}")

# VISUALIZATION

log_step("="*100)
log_step("(SPATIAL LEIDEN) VISUALIZATIONS")
log_step("="*100)

# Use the filtered dataset for visualizations
top_viz_file = os.path.join(marker_genes_dir, 'SPATIAL_LEIDEN_TOP50_PER_CLUSTER.csv')
if os.path.exists(top_viz_file):
    try:
        # Standard Stereopy visualizations
        data.plt.marker_genes_text(
            res_key='marker_genes_spldn',
            markers_num=10,
            sort_key='scores',
            out_path=os.path.join(viz_marker_dir, 'SPATIAL_LEIDEN_MARKERS.png')
        )
        
        data.plt.marker_genes_scatter(
            res_key='marker_genes_spldn', 
            markers_num=10,
            out_path=os.path.join(viz_marker_dir, 'SPATIAL_LEIDEN_MARKERS_SCATTER.png')
        )
        
        log_step("(SPATIAL LEIDEN) Visualizations completed")
        
    except Exception as e:
        log_step(f"(SPATIAL LEIDEN) WARNING: Error in visualizations: {e}")

# Generate volcano plots for ALL major clusters
try:
    log_step("(SPATIAL LEIDEN) Generating volcano plots")
    volcano_dir = os.path.join(viz_marker_dir, 'SPATIAL_LEIDEN_VOLCANO_PLOTS_COMPLETE')
    os.makedirs(volcano_dir, exist_ok=True)
    
    if 'spatial_leiden' in data.cells:
        all_clusters = sorted(data.cells['spatial_leiden'].unique())
        log_step(f"(SPATIAL LEIDEN) Generating volcano plots for ALL {len(all_clusters)} clusters")
        
        successful_plots = 0
        
        for i, cluster in enumerate(all_clusters):
            try:
                cluster_str = str(cluster)
                group_name = f'{cluster_str}.vs.rest'
                data.plt.marker_gene_volcano(
                    group_name=group_name,
                    res_key='marker_genes_spldn',
                    vlines=False,
                    out_path=os.path.join(volcano_dir, f'SPATIAL_LEIDEN_VOLCANO_CLUSTER_{cluster_str.zfill(2)}.png')
                )
                
                successful_plots += 1
                
                if (i + 1) % 20 == 0:
                    log_step(f"(SPATIAL LEIDEN) Generated {successful_plots}/{i + 1} volcano plots successfully")
                    
            except Exception as e:
                log_step(f"(SPATIAL LEIDEN) WARNING: Error generating volcano plot for cluster {cluster}: {e}")
                try:
                    alt_group_name = f'{int(cluster)}.vs.rest'
                    data.plt.marker_gene_volcano(
                        group_name=alt_group_name,
                        res_key='marker_genes_spldn',
                        vlines=False,
                        out_path=os.path.join(volcano_dir, f'SPATIAL_LEIDEN_VOLCANO_CLUSTER_{str(cluster).zfill(2)}_ALT.png')
                    )
                    successful_plots += 1
                    log_step(f"(SPATIAL LEIDEN) Alternative format worked for cluster {cluster}")
                except Exception as e2:
                    log_step(f"(SPATIAL LEIDEN) Both formats failed for cluster {cluster}: {e2}")
                continue
        
        log_step(f"(SPATIAL LEIDEN) Volcano plots completed: {successful_plots}/{len(all_clusters)} successful")
        
        # If no plots were generated, try a different approach
        if successful_plots == 0:
            log_step("(SPATIAL LEIDEN) No volcano plots generated. Trying alternative approach")
            
            # Check what group names are available in marker_genes result
            if 'marker_genes_spldn' in data.tl.result:
                available_groups = [k for k in data.tl.result['marker_genes_spldn'].keys() if k.endswith('.vs.rest')]
                log_step(f"(SPATIAL LEIDEN) Available group names in results: {available_groups[:10]}")
                
                # Try with first few available group names
                for j, group_name in enumerate(available_groups[:10]):
                    try:
                        cluster_id = group_name.replace('.vs.rest', '')
                        data.plt.marker_gene_volcano(
                            group_name=group_name,
                            res_key='marker_genes_spldn',
                            vlines=False,
                            out_path=os.path.join(volcano_dir, f'SPATIAL_LEIDEN_VOLCANO_CLUSTER_{cluster_id}_DIRECT.png')
                        )
                        successful_plots += 1
                        log_step(f"(SPATIAL LEIDEN) Direct approach successful for {group_name}")
                    except Exception as e:
                        log_step(f"(SPATIAL LEIDEN) Direct approach failed for {group_name}: {e}")
                        continue
            
            log_step(f"(SPATIAL LEIDEN) {successful_plots} volcano plots generated")
    
except Exception as e:
    log_step(f"(SPATIAL LEIDEN) WARNING: Error in volcano plot generation: {e}")

log_memory_usage("VISUALIZATIONS_COMPLETE")

# CUSTOM GENE ANNOTATION
log_step("(SPATIAL LEIDEN) Applying custom gene interest annotation")

custom_gene_markers = create_custom_gene_markers()
log_step(custom_gene_markers)
custom_annotations, custom_scores, annotation_details = apply_gene_interest_annotation(
    data, custom_gene_markers, cluster_res_key='spatial_leiden', threshold=1.2
)

# Create detailed report for custom gene annotation
custom_report_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'SPATIAL_LEIDEN_INTEREST_ANALYSIS_REPORT.txt')
with open(custom_report_file, 'w') as f:
    f.write("(SPATIAL LEIDEN) CUSTOM GENES OF INTEREST ANNOTATION REPORT\n")
    f.write("="*100 + "\n\n")
    
    f.write("GENES OF INTEREST DEFINED:\n")
    for category, genes in custom_gene_markers.items():
        f.write(f"\n{category.upper()}:\n")
        if genes:
            for gene in genes:
                f.write(f"  - {gene}\n")
        else:
            f.write("  - (No genes defined yet)\n")
    
    f.write(f"\nANNOTATION RESULTS (threshold: 1.2):\n")
    f.write("=" * 100 + "\n")
    
    for cluster in sorted(annotation_details.keys(), key=lambda x: int(x)):
        details = annotation_details[cluster]
        f.write(f"\nCLUSTER {cluster}:\n")
        f.write(f"  Annotation: {custom_annotations[cluster]}\n")
        f.write(f"  Cell count: {details['cell_count']}\n")
        f.write(f"  Best category: {details['best_category'] or 'None'}\n")
        f.write(f"  Best score: {details['best_score']:.3f}\n")

       	if details['all_scores']:
            f.write("  Category scores:\n")
            for category, score in details['all_scores'].items():
                f.write(f"    {category}: {score:.3f}\n")

log_step(f"(SPATIAL LEIDEN) Custom gene interest annotation completed for {len(custom_annotations)} clusters")

# TOP 50 MARKER GENES
marker_genes_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'SPATIAL_LEIDEN_TOP50_PER_CLUSTER.csv')
biological_annotations = {}

if os.path.exists(marker_genes_file):
    try:
        marker_df = pd.read_csv(marker_genes_file)
        log_step(f"(SPATIAL LEIDEN) Loaded {len(marker_df)} marker gene records")
        
        # Create biological annotations based on marker genes
        biological_annotations = create_biological_annotation_dict(marker_df, n_top_genes=5)
        log_step(f"(SPATIAL LEIDEN) Created annotations for {len(biological_annotations)} clusters")
        
        # Save biological annotation mapping
        bio_anno_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'SPATIAL_LEIDEN_INTEREST_ANALYSIS.json')
        with open(bio_anno_file, 'w') as f:
            json.dump(biological_annotations, f, indent=2)
        
        # Create marker-based annotation report
        marker_report_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'SPATIAL_LEIDEN_MARKERS_REPORT.txt')
        with open(marker_report_file, 'w') as f:
            f.write("MARKER-BASED BIOLOGICAL ANNOTATION REPORT\n")
            f.write("="*100 + "\n\n")
            
            for cluster in sorted(biological_annotations.keys(), key=lambda x: int(x)):
                f.write(f"CLUSTER {cluster}:\n")
                f.write(f"Annotation: {biological_annotations[cluster]}\n")
                
                # Get top markers for this cluster
                cluster_markers = marker_df[marker_df['cluster'] == int(cluster)]
                if len(cluster_markers) > 0:
                    if 'scores' in cluster_markers.columns:
                        top_markers = cluster_markers.nlargest(10, 'scores')
                    else:
                        top_markers = cluster_markers.head(10)
                    
                    f.write("Top marker genes:\n")
                    for _, row in top_markers.iterrows():
                        gene = row['genes'] if 'genes' in row else 'Unknown'
                        score = row['scores'] if 'scores' in row else 'N/A'
                        f.write(f"  - {gene}: {score}\n")
                
                f.write("\n")
        
    except Exception as e:
        log_step(f"(SPATIAL LEIDEN) WARNING: Error loading marker genes for annotation: {e}")

# MARKER GENE FILTERING

log_step("="*100)
log_step("(SPATIAL LEIDEN) MARKER GENE FILTERING")
log_step("="*100)

try:
    # Apply multiple filtering criteria
    filter_configs = [
        {'min_fold_change': 0.5, 'min_in_group_fraction': 0.2, 'max_out_group_fraction': 0.6, 'name': 'lenient'},
        {'min_fold_change': 1.0, 'min_in_group_fraction': 0.25, 'max_out_group_fraction': 0.5, 'name': 'moderate'},
        {'min_fold_change': 1.5, 'min_in_group_fraction': 0.3, 'max_out_group_fraction': 0.4, 'name': 'stringent'},
        {'min_fold_change': 2.0, 'min_in_group_fraction': 0.4, 'max_out_group_fraction': 0.3, 'name': 'very_stringent'}
    ]
    
    for config in filter_configs:
        try:
            filter_key = f"SPATIAL_LEIDEN_marker_genes_filtered_{config['name']}"
            
            data.tl.filter_marker_genes(
                marker_genes_res_key='marker_genes_spldn',
                min_fold_change=config['min_fold_change'],
                min_in_group_fraction=config['min_in_group_fraction'],
                max_out_group_fraction=config['max_out_group_fraction'],
                res_key=filter_key
            )
            
            log_step(f"(SPATIAL LEIDEN) Applied {config['name']} filtering (FC>{config['min_fold_change']})")
            
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
                        
                        filter_file = os.path.join(filtered_results_dir, f'SPATIAL_LEIDEN_FILTERED_MARKERS_{config["name"].upper()}.csv')
                        combined_filtered.to_csv(filter_file, index=False)
                        
                        log_step(f"(SPATIAL LEIDEN) Saved {len(combined_filtered)} {config['name']} filtered markers")
            
        except Exception as e:
            log_step(f"(SPATIAL LEIDEN) WARNING: Error in {config['name']} filtering: {e}")
            continue
    
    log_memory_usage("FILTERING_COMPLETE")
    
except Exception as e:
    log_step(f"(SPATIAL LEIDEN) WARNING: Error in advanced filtering: {e}")

# CLUSTER ANNOTATION FOR ALL METHODS
log_step("="*100)
log_step("CREATING ANNOTATIONS FOR ALL CLUSTERING METHODS")
log_step("="*100)

clustering_methods = ['louvain', 'leiden', 'spatial_leiden']
all_annotations = {}

for cluster_method in clustering_methods:
    if cluster_method not in data.cells:
        log_step(f"SKIPPING: {cluster_method} - clustering not found in data")
        continue
    
    log_step(f"\nGenerating annotations for {cluster_method.upper()} clustering")
    
    n_clusters = len(data.cells[cluster_method].unique())
    log_step(f"  Total clusters: {n_clusters}")
    
    method_annotations = create_annotations_for_cluster_method(data, cluster_method)
    all_annotations[cluster_method] = method_annotations
    
    # Save annotation files
    for anno_type, anno_dict in method_annotations.items():
        anno_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', f'{cluster_method}_{anno_type}_annotations'.upper() + '.json')
        with open(anno_file, 'w') as f:
            json.dump(anno_dict, f, indent=2)
        log_step(f"  Saved {cluster_method} {anno_type} annotations")

# APPLY ANNOTATIONS AND CREATE VISUALIZATIONS
log_step("="*100)
log_step("APPLYING ANNOTATIONS AND CREATING VISUALIZATIONS")
log_step("="*100)

for cluster_method in clustering_methods:
    if cluster_method not in data.cells or cluster_method not in all_annotations:
        continue
    
    log_step(f"\nApplying and visualizing {cluster_method.upper()} annotations")
    
    for anno_type, anno_dict in all_annotations[cluster_method].items():
        try:
            res_key = f'anno_{cluster_method}_{anno_type}'
            
            # Apply annotation
            data.tl.annotation(
                annotation_information=anno_dict,
                cluster_res_key=cluster_method,
                res_key=res_key
            )
            
            # Create spatial visualization
            plot_file = os.path.join(output_dir, 'PLOTS', 'ANNOTATION', 
                                    f'{cluster_method}_{anno_type}_spatial'.upper() + '.png')
            data.plt.cluster_scatter(res_key=res_key, out_path=plot_file)
            
            # Create UMAP visualization
            umap_plot_file = os.path.join(output_dir, 'PLOTS', 'ANNOTATION', 
                                         f'{cluster_method}_{anno_type}_umap'.upper() + '.png')
            data.plt.umap(res_key='umap', cluster_key=res_key, out_path=umap_plot_file)
            
            log_step(f"  {anno_type}: visualizations created")
            
        except Exception as e:
            log_step(f"  WARNING: Error with {anno_type}: {e}")
            continue

# CUSTOM GENE ANNOTATION
log_step("(SPATIAL LEIDEN) Applying custom gene interest annotation")

custom_gene_markers = create_custom_gene_markers()
log_step(custom_gene_markers)
custom_annotations, custom_scores, annotation_details = apply_gene_interest_annotation(
    data, custom_gene_markers, cluster_res_key='spatial_leiden', threshold=1.2
)

# Create detailed report for custom gene annotation
custom_report_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'LOUVAIN_INTEREST_ANALYSIS_REPORT.txt')
with open(custom_report_file, 'w') as f:
    f.write("(LOUVAIN) CUSTOM GENES OF INTEREST ANNOTATION REPORT\n")
    f.write("="*100 + "\n\n")
    
    f.write("GENES OF INTEREST DEFINED:\n")
    for category, genes in custom_gene_markers.items():
        f.write(f"\n{category.upper()}:\n")
        if genes:
            for gene in genes:
                f.write(f"  - {gene}\n")
        else:
            f.write("  - (No genes defined yet)\n")
    
    f.write(f"\nANNOTATION RESULTS (threshold: 1.2):\n")
    f.write("=" * 100 + "\n")
    
    for cluster in sorted(annotation_details.keys(), key=lambda x: int(x)):
        details = annotation_details[cluster]
        f.write(f"\nCLUSTER {cluster}:\n")
        f.write(f"  Annotation: {custom_annotations[cluster]}\n")
        f.write(f"  Cell count: {details['cell_count']}\n")
        f.write(f"  Best category: {details['best_category'] or 'None'}\n")
        f.write(f"  Best score: {details['best_score']:.3f}\n")
        
        if details['all_scores']:
            f.write("  Category scores:\n")
            for category, score in details['all_scores'].items():
                f.write(f"    {category}: {score:.3f}\n")

log_step(f"Custom gene interest annotation completed for {len(custom_annotations)} clusters")

# TOP50 CLUSTER
marker_genes_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'LOUVAIN_TOP50_PER_CLUSTER.csv')
biological_annotations = {}

if os.path.exists(marker_genes_file):
    try:
        log_step("Loading marker genes for biological annotation")
        marker_df = pd.read_csv(marker_genes_file)
        log_step(f"Loaded {len(marker_df)} marker gene records")
        
        # Create biological annotations based on marker genes
        biological_annotations = create_biological_annotation_dict(marker_df, n_top_genes=5)
        log_step(f"Created biological annotations for {len(biological_annotations)} clusters")
        
        # Save biological annotation mapping
        bio_anno_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'LOUVAIN_INTEREST_ANALYSIS.json')
        with open(bio_anno_file, 'w') as f:
            json.dump(biological_annotations, f, indent=2)
        
        # Create detailed marker-based annotation report
        marker_report_file = os.path.join(output_dir, 'INTEREST_ANALYSIS', 'LOUVAIN_MARKERS_REPORT.txt')
        with open(marker_report_file, 'w') as f:
            f.write("MARKER-BASED BIOLOGICAL ANNOTATION REPORT\n")
            f.write("="*100 + "\n\n")
            
            for cluster in sorted(biological_annotations.keys(), key=lambda x: int(x)):
                f.write(f"CLUSTER {cluster}:\n")
                f.write(f"Annotation: {biological_annotations[cluster]}\n")
                
                # Get top markers for this cluster
                cluster_markers = marker_df[marker_df['cluster'] == int(cluster)]
                if len(cluster_markers) > 0:
                    if 'scores' in cluster_markers.columns:
                        top_markers = cluster_markers.nlargest(10, 'scores')
                    else:
                        top_markers = cluster_markers.head(10)
                    
                    f.write("Top marker genes:\n")
                    for _, row in top_markers.iterrows():
                        gene = row['genes'] if 'genes' in row else 'Unknown'
                        score = row['scores'] if 'scores' in row else 'N/A'
                        f.write(f"  - {gene}: {score}\n")
                
                f.write("\n")
        
    except Exception as e:
        log_step(f"WARNING: Error loading marker genes for annotation: {e}")

# INDIVIDUAL CLUSTER PLOTS FOR ALL METHODS
log_step("="*100)
log_step("CREATING INDIVIDUAL CLUSTER PLOTS FOR ALL METHODS")
log_step("="*100)

for cluster_method in clustering_methods:
    if cluster_method not in data.cells:
        continue
    
    log_step(f"\nGenerating individual plots for {cluster_method.upper()}")
    
    try:
        # Get spatial coordinates
        if hasattr(data, 'position'):
            spatial_coords = data.position
            x_coords = spatial_coords[:, 0]
            y_coords = spatial_coords[:, 1]
        elif 'x' in data.cells.columns and 'y' in data.cells.columns:
            x_coords = data.cells['x'].values
            y_coords = data.cells['y'].values
        else:
            log_step(f"  ERROR: Could not find spatial coordinates for {cluster_method}")
            continue
        
        individual_dir = os.path.join(output_dir, 'CLUSTERING', cluster_method.upper())
        os.makedirs(individual_dir, exist_ok=True)
        
        unique_clusters = sorted(data.cells[cluster_method].unique())
        n_clusters = len(unique_clusters)
        
        log_step(f"  Generating {n_clusters} individual cluster plots...")
        
        for i, cluster in enumerate(unique_clusters):
            try:
                cluster_highlight = data.cells[cluster_method].copy()
                
                if hasattr(cluster_highlight, 'cat'):
                    if 'Other' not in cluster_highlight.cat.categories:
                        cluster_highlight = cluster_highlight.cat.add_categories(['Other'])
                
                cluster_highlight[cluster_highlight != cluster] = 'Other'
                
                fig, ax = plt.subplots(figsize=(10, 8))
                
                other_mask = cluster_highlight == 'Other'
                if other_mask.any():
                    ax.scatter(x_coords[other_mask], y_coords[other_mask], 
                              c='#CCCCCC', s=2.0, alpha=0.4, rasterized=True, edgecolors='none')
                
                cluster_mask = cluster_highlight == cluster
                if cluster_mask.any():
                    colors = plt.cm.Dark2(int(cluster) % 8)
                    ax.scatter(x_coords[cluster_mask], y_coords[cluster_mask], 
                              c=[colors], s=4.0, alpha=0.95, rasterized=True, edgecolors='white', linewidths=0.1)
                
                ax.set_xlabel('Spatial X (μm)', fontsize=12)
                ax.set_ylabel('Spatial Y (μm)', fontsize=12)
                ax.set_title(f'{cluster_method.capitalize()} Cluster {cluster} Spatial Distribution', 
                           fontsize=14, fontweight='bold')
                
                legend_elements = [
                    mpatches.Patch(color='lightgray', label='Other clusters'),
                    mpatches.Patch(color=colors, label=f'Cluster {cluster}')
                ]
                ax.legend(handles=legend_elements, loc='upper right', fontsize=10, bbox_to_anchor=(1.15, 1), frameon=True, fancybox=True)
                
                # Scale bar
                scale_bar_length = 2000
                x_range = x_coords.max() - x_coords.min()
                y_range = y_coords.max() - y_coords.min()
                scale_x_left = x_coords.min() + 0.05 * x_range
                scale_x_right = x_coords.max() - 0.05 * x_range - scale_bar_length
                scale_y = y_coords.min() + 0.05 * y_range
                scale_x = scale_x_right

                ax.plot([scale_x, scale_x + scale_bar_length], [scale_y, scale_y], 
                       'k-', linewidth=3)
                ax.text(scale_x + scale_bar_length/2, scale_y - 0.02 * y_range, 
                       '2.0mm', ha='center', va='top', fontsize=10, fontweight='bold',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, pad=0.3))
                
                # Statistics
                n_cells_cluster = cluster_mask.sum()
                total_cells = len(cluster_highlight)
                percentage = (n_cells_cluster / total_cells) * 100
                
                stats_text = f'Cells: {n_cells_cluster:,}\nPercentage: {percentage:.1f}%'
                ax.text(1.15, 0.5, stats_text, transform=ax.transAxes, 
                       verticalalignment='center', fontsize=10,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))
                
                ax.set_aspect('equal')
                plt.tight_layout()
                
                summary_file = os.path.join(individual_dir, f'{cluster_method}_CLUSTER_{cluster}_SPATIAL'.upper() + '.png')
                plt.savefig(summary_file, dpi=300, bbox_inches='tight', 
                           facecolor='white', edgecolor='none', pad_inches=0.2)
                plt.close()
                
                if (i + 1) % 10 == 0:
                    log_step(f"  Processed {i + 1}/{n_clusters} clusters")
                
            except Exception as e:
                log_step(f"  WARNING: Error with cluster {cluster}: {e}")
                continue
        
        # Generate HTML summary
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>{cluster_method.capitalize()} Individual Cluster Analysis</title>
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
            <h1>{cluster_method.capitalize()} Individual Cluster Spatial Distribution</h1>
            <div class="summary">
                <h3>Analysis Summary</h3>
                <p><strong>Total clusters found:</strong> {n_clusters}</p>
                <p><strong>Analysis type:</strong> {cluster_method.capitalize()} clustering with spatial visualization</p>
            </div>
            <div class="cluster-grid">
        """
        
        for cluster in unique_clusters:
            img_filename = f"{cluster_method.upper()}_CLUSTER_{cluster}_SPATIAL.png"
            html_content += f"""
                <div class="cluster-item">
                    <div class="cluster-title">{cluster_method.capitalize()} Cluster {cluster}</div>
                    <img src="{img_filename}" alt="{cluster_method.capitalize()} Cluster {cluster}">
                </div>
            """
        
        html_content += """
            </div>
        </body>
        </html>
        """
        
        html_file = os.path.join(individual_dir, f'{cluster_method.upper()}_SUMMARY.html')
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        log_step(f"  HTML summary created: {cluster_method.upper()}_SUUMARY.html")
        
    except Exception as e:
        log_step(f"  ERROR in individual plots for {cluster_method}: {e}")
        continue

# DATA EXPORT
log_step("="*100)
log_step("DATA EXPORT")
log_step("="*100)

export_start_time = datetime.now()

try:
    log_step("Converting to AnnData format")
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    h5ad_filename = f"{sample_id}_{timestamp}.h5ad"
    h5ad_path = os.path.join(output_dir, 'EXPORTS', h5ad_filename)
    log_step(f"Exporting to {h5ad_filename}")
    
    if 'umap' in data.tl.result:
        umap_result = data.tl.result['umap']
        log_step(f"UMAP result type: {type(umap_result)}")
        
        # Convert to DataFrame if it's numpy array
        if isinstance(umap_result, np.ndarray):
            data.tl.result['umap'] = pd.DataFrame(
                umap_result.astype(np.float32),
                columns=['UMAP_1', 'UMAP_2']
            )
        elif isinstance(umap_result, pd.DataFrame):
            # Ensure it's float32
            data.tl.result['umap'] = umap_result.astype(np.float32)
        else:
            log_step(f"WARNING: Unexpected UMAP type: {type(umap_result)}")
    
    # Get DataFrame from Cell object
    cells_df = data.cells.to_df()
    log_step(f"Retrieved cells DataFrame: {cells_df.shape}")
    
    # Track which columns were categorical
    categorical_cols = []
    for col in cells_df.columns:
        if hasattr(cells_df[col], 'cat'):
            categorical_cols.append(col)
            cells_df[col] = cells_df[col].astype(str)
    
    if categorical_cols:
        log_step(f"Converted {len(categorical_cols)} categorical columns: {', '.join(categorical_cols[:5])}")
    
    adata = st.io.stereo_to_anndata(data, flavor='scanpy')
    log_step(f"AnnData object created: {adata}")
    
    if hasattr(adata, 'obsm'):
        log_step(f"Post-processing {len(adata.obsm)} obsm entries...")
        for key in list(adata.obsm.keys()):
            try:
                # Convert to dense if sparse
                if issparse(adata.obsm[key]):
                    adata.obsm[key] = adata.obsm[key].toarray()
                
                # Ensure float32
                if not isinstance(adata.obsm[key], np.ndarray):
                    adata.obsm[key] = np.array(adata.obsm[key], dtype=np.float32)
                else:
                    adata.obsm[key] = adata.obsm[key].astype(np.float32)
                
                log_step(f"{key}: shape {adata.obsm[key].shape}, dtype {adata.obsm[key].dtype}")
                
            except Exception as e:
                log_step(f"Could not process obsm['{key}']: {e}")
                try:
                    del adata.obsm[key]
                except:
                    pass
    
    # Get clustering info from the DataFrame
    clustering_info = {}
    for method in ['louvain', 'leiden', 'spatial_leiden']:
        if method in cells_df.columns:
            n_clusters = int(cells_df[method].nunique())
            clustering_info[method] = n_clusters
    
    adata.uns['analysis_metadata'] = {
        'analysis_date': datetime.now().isoformat(),
        'stereopy_version': st.__version__,
        'data_source': os.path.basename(data_path),
        'bin_size': int(BIN_SIZE),
        'preprocessing': {
            'min_counts': int(MIN_COUNTS),
            'min_genes': int(MIN_GENES),
            'pct_counts_mt': int(PCT_COUNTS_MT),
            'normalization': 'total_count_None_log1p',
            'scaling': 'max10_no_zero_center',
            'spatial_filter': f'X({MIN_X}-{MAX_X}) Y({MIN_Y}-{MAX_Y})',
            'hvg_params': {
                'min_mean': $HVG_MIN_MEAN,
                'max_mean': $HVG_MAX_MEAN,
                'min_disp': $HVG_DISP,
                'n_top_genes': $HVG_TOP
            }
        },
        'clustering_methods': clustering_info
    }
    
    adata.write_h5ad(h5ad_path, compression='gzip')
    file_size_mb = os.path.getsize(h5ad_path) / 1e6
    log_step(f"AnnData export completed. File size: {file_size_mb:.1f} MB")

    try:
        log_step("Exporting metadata CSV")
        metadata_csv = os.path.join(output_dir, 'EXPORTS', 'CELL_METADATA.csv')
        
        # Start with cells DataFrame
        metadata_df = cells_df.copy()
        
        # Add cell names if available
        if hasattr(data.cells, 'cell_name'):
            metadata_df.insert(0, 'cell_name', data.cells.cell_name)
        
        # Add spatial coordinates
        if hasattr(data, 'position'):
            spatial_coords = data.position
            metadata_df['spatial_x'] = spatial_coords[:, 0]
            metadata_df['spatial_y'] = spatial_coords[:, 1]

        # Add UMAP coordinates
        if 'umap' in data.tl.result:
            umap_coords = data.tl.result['umap']
            metadata_df['umap_1'] = umap_coords.iloc[:, 0]
            metadata_df['umap_2'] = umap_coords.iloc[:, 1]
        
        # Add QC metrics if available
        if hasattr(data.cells, 'total_counts'):
            metadata_df['total_counts'] = data.cells.total_counts
        if hasattr(data.cells, 'n_genes_by_counts'):
            metadata_df['n_genes'] = data.cells.n_genes_by_counts
        if hasattr(data.cells, 'pct_counts_mt'):
            metadata_df['pct_mt'] = data.cells.pct_counts_mt
        
        metadata_df.to_csv(metadata_csv, index=False)
        log_step(f"Metadata CSV saved: {len(metadata_df)} cells, {len(metadata_df.columns)} columns")

        # Export gene metadata
        gene_metadata_csv = os.path.join(output_dir, 'EXPORTS', 'MARKERS_METADATA.csv')
        if hasattr(data.genes, 'gene_name'):
            gene_df = pd.DataFrame({
                'gene_name': data.genes.gene_name
            })
        
            # Add highly variable genes info if available
            if 'highly_variable_genes' in data.tl.result:
                hvg_result = data.tl.result['highly_variable_genes']
                if isinstance(hvg_result, pd.DataFrame):
                    gene_df['highly_variable'] = False
                    gene_df.loc[gene_df['gene_name'].isin(hvg_result.index), 'highly_variable'] = True
                    hvg_result_subset = hvg_result[['means', 'dispersions']]
                    gene_df = gene_df.merge(hvg_result_subset, left_on='gene_name', right_index=True, how='left')

            gene_df.to_csv(gene_metadata_csv, index=False)
            log_step(f"Gene metadata saved: {len(gene_df)} genes")

        # Export expression matrix
        if 'highly_variable_genes' in data.tl.result:
            # Get highly variable genes
            hvg_result = data.tl.result['highly_variable_genes']
            if isinstance(hvg_result, pd.DataFrame) and 'highly_variable' in hvg_result.columns:
                hvg_genes = hvg_result[hvg_result['highly_variable'] == True].index.tolist()
                log_step(f"Found {len(hvg_genes)} highly variable genes (HVGs).")
                log_step(f"Top 5: {hvg_genes[:5]}")
            
                if len(hvg_genes) == 0:
                    log_step("No highly variable genes (HVGs) found")
                else:
                    try:
                        gene_name_to_idx = {gene: idx for idx, gene in enumerate(data.genes.gene_name)}
                        gene_indices = []
                        valid_hvg_genes = []
                        
                        for gene in hvg_genes:
                            if gene in gene_name_to_idx:
                                gene_indices.append(gene_name_to_idx[gene])
                                valid_hvg_genes.append(gene)
                        
                        log_step(f"Found {len(gene_indices)} highly variable genes (HVGs) of {len(hvg_genes)}")
                    
                        hvg_sparse_matrix = data.exp_matrix[:, gene_indices]
                        hvg_matrix_dense = hvg_sparse_matrix.toarray()
                        
                        hvg_expression_matrix_df = pd.DataFrame(
                            hvg_matrix_dense,
                            index=data.cells.cell_name,
                            columns=valid_hvg_genes
                        )
                    
                        expression_file = os.path.join(output_dir, 'EXPORTS', 'HVG_EXPRESSION_MATRIX.csv')
                        
                        hvg_expression_matrix_df.to_csv(expression_file)
                        log_step(f"HVG matrix saved in {expression_file} with {hvg_expression_matrix_df.shape[0]} cells and {hvg_expression_matrix_df.shape[1]} genes")
                    
                    except Exception as e1:
                        log_step(f"ERROR exporting HVG matrix: {e1}")
                        try:
                            log_step("Trying fallback")
                            
                            # Find index
                            gene_indices = []
                            for gene in hvg_genes[:100]:
                                try:
                                    idx = data.genes.gene_name.tolist().index(gene)
                                    gene_indices.append(idx)
                                except ValueError:
                                    continue
                        
                            if len(gene_indices) > 0:
                                log_step(f"Exporting {len(gene_indices)} HVGs")
                                submatrix = data.exp_matrix[:, gene_indices]
                                hvg_df = pd.DataFrame(
                                    submatrix.toarray(),
                                    index=data.cells.cell_name,
                                    columns=[hvg_genes[hvg_genes.index(data.genes.gene_name.iloc[idx])] for idx in gene_indices if data.genes.gene_name.iloc[idx] in hvg_genes]
                                )
                                expression_file = os.path.join(output_dir, 'EXPORTS', 'HVG_EXPRESSION_MATRIX_DF.csv')
                                hvg_df.to_csv(expression_file)
                                log_step(f"Fallback success: File saved with {hvg_df.shape[0]} cells and {hvg_df.shape[1]} genes")
                            else:
                                log_step("No valid HVG found for exporting")
                            
                        except Exception as e2:
                            log_step(f"Fallback method failed: {e2}")
                            log_step("Skipping HVG export")
            else:
                log_step("ERROR: No valid HVG in the dataframe.")

    except Exception as e:
        log_step(f"WARNING: Error in CSV export: {e}")
        import traceback
        log_step(f"CSV export traceback: {traceback.format_exc()}")
    
except Exception as e:
    log_step(f"ERROR in data export: {e}")
    import traceback
    log_step(f"Full export traceback: {traceback.format_exc()}")

log_memory_usage("EXPORT_COMPLETE")

# Visualization

log_step("="*100)
log_step("CREATING DIRECT GENE INTEREST VISUALIZATION")
log_step("="*100)
direct_results = create_direct_gene_visualization()

if direct_results:
    log_step("SUCCESS: Direct gene visualization completed!")
else:
    log_step("FAILED: Could not create direct visualization")

# Individual analysis
log_step("="*100)
log_step("INDIVIDUAL GENES ANALYSIS")
log_step("="*100)

interest_genes_path = "$INTEREST_GENES_PATH"
expression_thr = $EXPRESSION_THR

if interest_genes_path and os.path.exists(interest_genes_path):
    # Parse gene list (format: category,gene1,gene2,gene3)
    log_step(f"Loading genes from {interest_genes_path}")
    interest_genes = []
    try:
        with open(interest_genes_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    # Split by comma and take all genes (skip first column which is category)
                    parts = line.split(',')
                    if len(parts) > 1:
                        genes = [g.strip() for g in parts[1:] if g.strip()]
                        interest_genes.extend(genes)
        
        # Remove duplicates while preserving order
        seen = set()
        interest_genes = [g for g in interest_genes if not (g in seen or seen.add(g))]
        
        log_step(f"Loaded {len(interest_genes)} unique genes: {interest_genes}")
        
        if not interest_genes:
            log_step("Warning: No genes found in file. Skipping individual analysis.")
        else:
            # Create output directories
            individual_analysis_dir = os.path.join(output_dir, 'INTEREST_ANALYSIS')
            os.makedirs(individual_analysis_dir, exist_ok=True)
            
            maps_dir = os.path.join(individual_analysis_dir, 'MAPS')
            os.makedirs(maps_dir, exist_ok=True)
            
            total_counts_dir = os.path.join(maps_dir, 'TOTAL_COUNTS')
            os.makedirs(total_counts_dir, exist_ok=True)
            
            for gene in interest_genes:
                gene_dir = os.path.join(maps_dir, gene)
                os.makedirs(gene_dir, exist_ok=True)
            
            # Get sample name
            sample_name = os.path.splitext(os.path.basename(data_path))[0]
            
            # Calculate total counts
            total_counts_array = data.exp_matrix.sum(axis=1).A1
            tissue_spots_mask = total_counts_array > 0
            
            x_tissue = data.position[tissue_spots_mask, 0]
            y_tissue = data.position[tissue_spots_mask, 1]
            
            sample_stats = {
                'sample_name': sample_name,
                'total_tissue_spots': int(np.sum(tissue_spots_mask)),
                'genes': {}
            }
            
            # PLOT: Total Counts Map
            log_step("Generating total counts map")
            fig, ax = plt.subplots(figsize=(10, 10), facecolor='black')
            ax.set_facecolor('black')
            
            if np.sum(tissue_spots_mask) > 0:
                total_counts_tissue = total_counts_array[tissue_spots_mask]
                norm_total = Normalize(vmin=0, vmax=np.max(total_counts_tissue))
                scatter_total = ax.scatter(x_tissue, y_tissue, 
                                           c=total_counts_tissue, cmap='viridis', 
                                           s=25, norm=norm_total, edgecolors='none')
                cbar = plt.colorbar(scatter_total, ax=ax, label='Total Expression')
                cbar.ax.yaxis.set_tick_params(color='white')
                plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
                cbar.set_label('Total Expression', color='white')
            
            ax.set_title(f"Total Counts - {sample_name}", fontsize=16, color='white')
            ax.set_xlabel('Spatial X (μm)', color='white')
            ax.set_ylabel('Spatial Y (μm)', color='white')
            ax.tick_params(colors='white')
            for spine in ax.spines.values():
                spine.set_color('white')
            
            path_counts = os.path.join(total_counts_dir, f"{sample_name}_total_counts_thr{expression_thr}.png")
            plt.savefig(path_counts, dpi=300, bbox_inches='tight', facecolor='black')
            plt.close(fig)
            log_step(f"Total counts map saved: {path_counts}")
            
            # Process individual genes
            for gene in interest_genes:
                log_step(f"Processing gene: {gene}")
                
                gene_indices = np.where(data.genes.gene_name == gene)[0]
                
                if len(gene_indices) > 0:
                    gene_index = gene_indices[0]
                    gene_expression_vector = data.exp_matrix[:, gene_index].toarray().flatten()
                    
                    fig, ax = plt.subplots(figsize=(12, 10), facecolor='black')
                    ax.set_facecolor('black')
                    
                    # Layer 1: Background (grayscale)
                    if np.sum(tissue_spots_mask) > 0:
                        total_counts_tissue = total_counts_array[tissue_spots_mask]
                        norm_background = Normalize(vmin=0, vmax=np.max(total_counts_tissue))
                        
                        scatter_bg = ax.scatter(x_tissue, y_tissue, 
                                               c=total_counts_tissue, 
                                               cmap='gray',
                                               s=15, 
                                               norm=norm_background, 
                                               alpha=0.5,
                                               edgecolors='none',
                                               label='Tissue Background')
                    
                    # Layer 2: Gene expression
                    expressed_gene_spots_mask = (gene_expression_vector > expression_thr) & tissue_spots_mask
                    
                    if np.sum(expressed_gene_spots_mask) > 0:
                        x_expressed = data.position[expressed_gene_spots_mask, 0]
                        y_expressed = data.position[expressed_gene_spots_mask, 1]
                        expression_values = gene_expression_vector[expressed_gene_spots_mask]
                        
                        norm_gene = Normalize(vmin=0, vmax=np.max(expression_values))
                        scatter_gene = ax.scatter(x_expressed, y_expressed,
                                                  c=expression_values, 
                                                  cmap='plasma',
                                                  s=35,
                                                  norm=norm_gene, 
                                                  edgecolors='white', 
                                                  linewidths=0.3, 
                                                  label=f'{gene}',
                                                  alpha=0.95)
                        
                        cbar = plt.colorbar(scatter_gene, ax=ax, 
                                           label=f'{gene} Expression (>{expression_thr})')
                        cbar.ax.yaxis.set_tick_params(color='white')
                        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
                        cbar.set_label(f'{gene} Expression (>{expression_thr})', color='white')
                        
                        expressed_spots = int(np.sum(expressed_gene_spots_mask))
                        max_expr = float(np.max(expression_values))
                        min_expr = float(np.min(expression_values))
                        mean_expr = float(np.mean(expression_values))
                        
                        log_step(f"{gene}: {expressed_spots} spots expressed, max: {max_expr:.2f}")
                        
                        sample_stats['genes'][gene] = {
                            'expressed_spots': expressed_spots,
                            'max_expression': max_expr,
                            'min_expression': min_expr,
                            'mean_expression': mean_expr,
                            'percentage_of_tissue': (expressed_spots / sample_stats['total_tissue_spots'] * 100) if sample_stats['total_tissue_spots'] > 0 else 0
                        }
                        
                        # Save individual stats
                        stats_file = os.path.join(maps_dir, gene, f"{sample_name}_{gene}_stats_thr{expression_thr}.txt")
                        with open(stats_file, 'w') as f:
                            f.write(f"Gene: {gene}\n")
                            f.write(f"Sample: {sample_name}\n")
                            f.write(f"Analysis date: {datetime.now().strftime('%Y/%m/%d %H:%M:%S')}\n")
                            f.write(f"Expression threshold: > {expression_thr}\n\n")
                            f.write(f"Total tissue spots: {sample_stats['total_tissue_spots']}\n")
                            f.write(f"Expressed spots: {expressed_spots}\n")
                            f.write(f"Percentage of tissue: {(expressed_spots / sample_stats['total_tissue_spots'] * 100):.2f}%\n")
                            f.write(f"Maximum expression: {max_expr:.2f}\n")
                            f.write(f"Minimum expression: {min_expr:.2f}\n")
                            f.write(f"Mean expression: {mean_expr:.2f}\n")
                        
                    else:
                        log_step(f"{gene} not expressed above threshold {expression_thr}")
                        sample_stats['genes'][gene] = {
                            'expressed_spots': 0,
                            'max_expression': 0,
                            'min_expression': 0,
                            'mean_expression': 0,
                            'percentage_of_tissue': 0
                        }
                        
                        stats_file = os.path.join(maps_dir, gene, f"{sample_name}_{gene}_stats_thr{expression_thr}.txt")
                        with open(stats_file, 'w') as f:
                            f.write(f"Gene: {gene}\n")
                            f.write(f"Sample: {sample_name}\n")
                            f.write(f"Status: Not expressed above threshold ({expression_thr})\n")
                        
                        ax.text(0.5, 0.95, f'{gene} not detected\n(expression ≤ {expression_thr})', 
                               transform=ax.transAxes, fontsize=12, 
                               bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7),
                               ha='center', va='top', color='black')
                    
                    ax.set_title(f"Spatial Expression: {gene} - {sample_name}\n(Threshold > {expression_thr})", 
                                fontsize=14, color='white', pad=20)
                    ax.set_xlabel('Spatial X (μm)', fontsize=12, color='white')
                    ax.set_ylabel('Spatial Y (μm)', fontsize=12, color='white')
                    ax.tick_params(colors='white')
                    for spine in ax.spines.values():
                        spine.set_color('white')
                    
                    if np.sum(expressed_gene_spots_mask) > 0:
                        legend = ax.legend(loc='upper left', bbox_to_anchor=(0, 1), 
                                          framealpha=0.8, facecolor='black', edgecolor='white')
                        for text in legend.get_texts():
                            text.set_color('white')
                    
                    ax.grid(False)
                    plt.tight_layout()
                    
                    plot_path = os.path.join(maps_dir, gene, f"{sample_name}_{gene}_thr{expression_thr}.png")
                    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='black')
                    plt.close(fig)
                    log_step(f"Plot saved: {plot_path}")
                    
                else:
                    log_step(f"{gene} not found in dataset")
                    sample_stats['genes'][gene] = 'Gene not found'
                    
                    stats_file = os.path.join(maps_dir, gene, f"{sample_name}_{gene}_stats_thr{expression_thr}.txt")
                    with open(stats_file, 'w') as f:
                        f.write(f"Gene: {gene}\n")
                        f.write(f"Sample: {sample_name}\n")
                        f.write(f"Status: Gene not found in dataset\n")
            
            # Save global stats
            stats_file = os.path.join(individual_analysis_dir, f"expression_stats_thr{expression_thr}.txt")
            with open(stats_file, 'w') as f:
                f.write("INDIVIDUAL GENE EXPRESSION ANALYSIS\n")
                f.write(f"Analysis date: {datetime.now().strftime('%Y/%m/%d %H:%M:%S')}\n")
                f.write(f"Sample: {sample_name}\n")
                f.write(f"Expression threshold: > {expression_thr}\n")
                f.write(f"Total tissue spots: {sample_stats['total_tissue_spots']}\n\n")
                f.write(f"Genes analyzed: {len(interest_genes)}\n")
                f.write(f"Gene list: {', '.join(interest_genes)}\n\n")
                f.write("="*80 + "\n\n")
                
                for gene, stats in sample_stats['genes'].items():
                    f.write(f"GENE: {gene}\n")
                    if isinstance(stats, dict):
                        f.write(f"  Expressed spots: {stats['expressed_spots']}\n")
                        f.write(f"  Percentage of tissue: {stats['percentage_of_tissue']:.2f}%\n")
                        if stats['expressed_spots'] > 0:
                            f.write(f"  Maximum expression: {stats['max_expression']:.2f}\n")
                            f.write(f"  Minimum expression: {stats['min_expression']:.2f}\n")
                            f.write(f"  Mean expression: {stats['mean_expression']:.2f}\n")
                        else:
                            f.write(f"  Expression: Below threshold\n")
                    else:
                        f.write(f"  Status: {stats}\n")
                    f.write("\n")
            
            log_step(f"Individual gene analysis completed. Results in {individual_analysis_dir}")
            
    except Exception as e:
        log_step(f"ERROR in individual gene analysis: {e}")
        import traceback
        traceback.print_exc()
        
else:
    log_step("No interest genes file provided or file doesn't exist. Skipping individual gene analysis.")

# DATA INTEGRATION

if direct_results:
    log_step("SUCCESS: Direct gene visualization completed!")
    clean_spatial_success = create_clean_spatial_visualization_from_direct_results(
        direct_results, data, output_dir
    )
    
    if clean_spatial_success:
        log_step("SUCCESS: Clean spatial visualization integrated!")
    else:
        log_step("WARNING: Clean spatial visualization failed")
else:
    log_step("FAILED: Could not create direct visualization")

# FINAL REPORT

log_step("="*100)
log_step("GENERATING FINAL REPORT")
log_step("="*100)

end_time = datetime.now()
total_duration = (end_time - start_time).total_seconds()

# Generate Final report
final_report_path = os.path.join(output_dir, 'ANALYSIS_REPORT.txt')
with open(final_report_path, 'w') as f:
    f.write("="*100 + "\n")
    f.write("SPATIAL TRANSCRIPTOMICS ANALYSIS REPORT\n")
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

    if 'leiden' in data.cells:
        cluster_counts = data.cells['leiden'].value_counts()
        f.write("LEIDEN CLUSTERING ANALYSIS RESULTS:\n")
        f.write(f"- Total clusters identified: {len(cluster_counts)}\n")
        f.write(f"- Total cells analyzed: {cluster_counts.sum()}\n")
        f.write(f"- Largest cluster: {cluster_counts.max()} cells (Cluster {cluster_counts.idxmax()})\n")
        f.write(f"- Smallest cluster: {cluster_counts.min()} cells (Cluster {cluster_counts.idxmin()})\n")
        f.write(f"- Average cluster size: {cluster_counts.mean():.1f} cells\n")
        f.write(f"- Total genes analyzed: {data.n_genes}\n\n")

    if 'spatial_leiden' in data.cells:
        cluster_counts = data.cells['spatial_leiden'].value_counts()
        f.write("SPATIAL_LEIDEN CLUSTERING ANALYSIS RESULTS:\n")
        f.write(f"- Total clusters identified: {len(cluster_counts)}\n")
        f.write(f"- Total cells analyzed: {cluster_counts.sum()}\n")
        f.write(f"- Largest cluster: {cluster_counts.max()} cells (Cluster {cluster_counts.idxmax()})\n")
        f.write(f"- Smallest cluster: {cluster_counts.min()} cells (Cluster {cluster_counts.idxmin()})\n")
        f.write(f"- Average cluster size: {cluster_counts.mean():.1f} cells\n")
        f.write(f"- Total genes analyzed: {data.n_genes}\n\n")
    
# Save processing summary (LOUVAIN)
processing_summary_file = os.path.join(output_dir, 'LOGS', 'LOUVAIN_SUMMARY.json')
processing_summary = {
    'analysis_type': 'NO_LIMITS',
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

log_step(f"Processing summary (LOUVAIN) saved to {processing_summary_file}")

# Save processing summary (LEIDEN)
processing_summary_file = os.path.join(output_dir, 'LOGS', 'LEIDEN_SUMMARY.json')
processing_summary = {
    'analysis_type': 'NO_LIMITS',
    'start_time': start_time.isoformat(),
    'end_time': end_time.isoformat(),
    'total_duration_seconds': total_duration,
    'total_duration_minutes': total_duration/60,
    'data_file': data_path,
    'n_clusters': len(data.cells['leiden'].value_counts()) if 'leiden' in data.cells else 0,
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

log_step(f"Processing summary (LEIDEN) saved to {processing_summary_file}")

# Save processing summary (SPATIAL_LEIDEN)
processing_summary_file = os.path.join(output_dir, 'LOGS', 'SPATIAL_LEIDEN_SUMMARY.json')
processing_summary = {
    'analysis_type': 'NO_LIMITS',
    'start_time': start_time.isoformat(),
    'end_time': end_time.isoformat(),
    'total_duration_seconds': total_duration,
    'total_duration_minutes': total_duration/60,
    'data_file': data_path,
    'n_clusters': len(data.cells['spatial_leiden'].value_counts()) if 'spatial_leiden' in data.cells else 0,
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

log_step(f"Processing summary (SPATIAL_LEIDEN) saved to {processing_summary_file}")

log_step(f"analysis report saved to {final_report_path}")
log_memory_usage("FINAL")

# Final cleanup
gc.collect()

log_step("="*100)
log_step("ANALYSIS COMPLETED SUCCESSFULLY!")
log_step(f"Total processing time: {total_duration/60:.1f} minutes")
log_step(f"Peak memory usage: {psutil.virtual_memory().percent:.1f}%")
log_step(f"Results location: {output_dir}")
log_step("="*100)

EOF
echo "==========================================="
echo "Primary Script created successfully!"
echo "==========================================="


# Generate analysis script
echo "==========================================="
echo "Creating Secondary Script!"
echo "==========================================="
cat > bin/SCRIPT_SECONDARY_ANALYSIS.py << EOF
#!/usr/bin/env python3
# Import dependencies
import os
from datetime import datetime

# Find the latest RESULTS folder
results_folders = [d for d in os.listdir('.') if d.startswith('RESULTS_') and os.path.isdir(d)]
if not results_folders:
    print("ERROR: No RESULTS folder found!")
    exit(1)

latest_results = sorted(results_folders)[-1]
print(f"Working on: {latest_results}")

# Generate Summary File
summary_file = os.path.join(latest_results, 'SECONDARY_ANALYSIS_SUMMARY.txt')
with open(summary_file, 'w') as f:
    f.write("=" * 50 + "\n")
    f.write("SECONDARY ANALYSIS\n")
    f.write("=" * 50 + "\n")
    f.write(f"Execution time: {datetime.now()}\n")
    f.write(f"Results folder: {latest_results}\n")
    f.write("\n")
    f.write("Contents of results folder:\n")
    
    # List all files in the results folder
    for item in sorted(os.listdir(latest_results)):
        item_path = os.path.join(latest_results, item)
        if os.path.isfile(item_path):
            size = os.path.getsize(item_path)
            f.write(f"  [FILE] {item} ({size} bytes)\n")
        elif os.path.isdir(item_path):
            f.write(f"  [DIR]  {item}/\n")

# Function 1

# Finish script
with open(summary_file, 'a') as f:
    f.write("\n" + "=" * 50 + "\n")
    f.write("Secondary Analysis completed successfully!\n")
    f.write("=" * 50 + "\n")
print("Secondary analysis test completed!")
EOF
echo "==========================================="
echo "Secondary Script created successfully!"
echo "==========================================="


echo "==========================================="
echo "Creating Converter Script..."
echo "==========================================="
cat > bin/SCRIPT_CONVERTER.py << EOF
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
INPUT_PATH = "$INPUT_PATH"
OUTPUT_DIR = "$OUTPUT_DIR"
BIN_SIZE = $BIN_SIZE

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
    
    log_step(f"Processing: {os.path.basename(path)}")
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
        log_step(f"[FINISH] Created {output_path} in {time.time() - start_time:.2f}s")
    except Exception as e:
        log_step(f"[ERROR] Failed to process {path}: {e}")

if __name__ == "__main__":
    # 1. Check if INPUT_PATH is empty
    if not INPUT_PATH:
        log_step("[ERROR] No INPUT_PATH provided.")
        sys.exit(0)

    # 2. Check if Path exists
    if not os.path.exists(INPUT_PATH):
        log_step(f"[ERROR] Path '{INPUT_PATH}' not found.")
        sys.exit(1)

    # 3. Process File or Directory
    if os.path.isfile(INPUT_PATH):
        convert_file(INPUT_PATH)
    elif os.path.isdir(INPUT_PATH):
        files = glob.glob(os.path.join(INPUT_PATH, "*.gem*")) + glob.glob(os.path.join(INPUT_PATH, "*.h5ad"))
        if not files:
            log_step(f"[ERROR] No compatible files (GEM/H5AD) found in {INPUT_PATH}")
        for f in files:
            convert_file(f)
EOF

echo "==========================================="
echo "Creating Converter Script..."
echo "==========================================="

# Generate NETWORK (R) analysis script
echo "==========================================="
echo "Creating NETWORK Analysis Script (R)..."
echo "==========================================="
cat > bin/SCRIPT_NETWORK_ANALYSIS.r << 'EOF'
lib_path <- paste0(getwd(), "/R_libs")
if(!dir.exists(lib_path)) dir.create(lib_path)
.libPaths(c(lib_path, .libPaths()))

# Check and install packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = lib_path, repos = "https://cloud.r-project.org")
  }
}

# Dependencies
needed_pkgs <- c("Seurat", "WGCNA", "tidyverse", "Matrix")
lapply(needed_pkgs, install_if_missing)

library(Seurat)
library(WGCNA)
library(tidyverse)
library(Matrix)

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 10)

# Find lastest FOLDER
all_dirs <- list.dirs("..", full.names = TRUE, recursive = FALSE)
results_dirs <- all_dirs[grepl("/RESULTS_", all_dirs)]

if (length(results_dirs) == 0) {
    stop("ERROR: No RESULTS folder found (RESULTS_*)")
}

latest_results <- sort(results_dirs, decreasing = TRUE)[1]
cat(paste0("Working on: ", latest_results, "\n"))

# Set working directory to the results folder
setwd(latest_results)

# Setup Network subdirectory
output_dir <- "NETWORK/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
project_name <- sub("RESULTS_", "", basename(latest_results))
cat(paste0("Project Name: ", project_name, "\n"))
checkpoint_file <- paste0(output_dir, "WGCNA_INTERMEDIATE_DATA.rds")

if (file.exists(checkpoint_file)) {
    cat("\n[CHECKPOINT] Loading previous progress\n")
    checkpoint_data <- readRDS(checkpoint_file)
    datExpr <- checkpoint_data$datExpr
    TOM <- checkpoint_data$TOM
    geneTree <- checkpoint_data$geneTree
    mergedColors <- checkpoint_data$mergedColors
    mergedMEs <- checkpoint_data$mergedMEs
} else {
    cat("\n[PROCESS] Starting full processing\n")
    
    # Load data from within the results folder
    hvg_matrix_path <- "EXPORTS/HVG_EXPRESSION_MATRIX.csv"

    if (!file.exists(hvg_matrix_path)) {
        stop(paste("ERROR: HVG Matrix not found at", hvg_matrix_path))
    }

    cat("\n[PROCESS] Loading HVG matrix from EXPORTS\n")
    
    counts_raw <- read.csv(hvg_matrix_path, row.names = 1, check.names = FALSE)
    counts_seurat <- t(as.matrix(counts_raw))
    seurat_obj <- CreateSeuratObject(counts = counts_seurat, project = project_name)
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

    # Metacells
    k_value <- 25
    pca_embeddings <- Embeddings(seurat_obj, reduction = "pca")[, 1:30]
    set.seed(12345)
    km_result <- kmeans(pca_embeddings, centers = k_value, iter.max = 100, nstart = 25)
    expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
    metacell_expr <- matrix(0, nrow = nrow(expr_matrix), ncol = k_value)
    rownames(metacell_expr) <- rownames(expr_matrix)
    colnames(metacell_expr) <- paste0("MC_", 1:k_value)
    for(i in 1:k_value) {
        cluster_cells <- which(km_result$cluster == i)
        if(length(cluster_cells) > 0) {
            metacell_expr[, i] <- Matrix::rowMeans(expr_matrix[, cluster_cells, drop = FALSE])
        }
    }

    # WGCNA
    datExpr <- as.data.frame(t(metacell_expr))
    vars <- apply(datExpr, 2, var)
    bad_genes <- names(vars[vars == 0 | is.na(vars)])
    if(length(bad_genes) > 0) {
        write.table(bad_genes, paste0(output_dir, "INITIAL_DISCARD_GSG.txt"), row.names=F, col.names=F, quote=F)
        datExpr <- datExpr[, !colnames(datExpr) %in% bad_genes]
    }

    # Thresholding and TOM
    powers <- c(seq(1, 10, by = 1), seq(12, 30, by = 2))
    sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "unsigned", verbose = 5)
    selected_power <- sft$powerEstimate
    if(is.na(selected_power)) selected_power <- 6

    adjacency <- adjacency(datExpr, power = selected_power, type = "unsigned")
    TOM <- TOMsimilarity(adjacency)
    dissTOM <- 1 - TOM
    geneTree <- hclust(as.dist(dissTOM), method = "average")
    dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
    dynamicColors <- labels2colors(dynamicMods)
    merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)
    mergedColors <- merge$colors
    mergedMEs <- merge$newMEs

    saveRDS(list(datExpr=datExpr, TOM=TOM, geneTree=geneTree, mergedColors=mergedColors, mergedMEs=mergedMEs), checkpoint_file)
}

# Export Results
kME <- cor(datExpr, mergedMEs, use = "p")
modules_df <- data.frame(gene_name = colnames(datExpr), module = mergedColors, color = mergedColors, stringsAsFactors = FALSE)
modules_df$kME <- sapply(1:nrow(modules_df), function(i) {
    mod <- modules_df$module[i]
    me_col <- paste0("ME", mod)
    if(me_col %in% colnames(kME)) return(kME[i, me_col]) else return(NA)
})

threshold <- 0.15
all_edges_list <- list()
for(mod in unique(mergedColors)) {
    mod_genes <- modules_df %>% filter(module == mod) %>% pull(gene_name)
    mod_idx <- which(colnames(datExpr) %in% mod_genes)
    if(length(mod_genes) < 2) next
    tom_sub <- TOM[mod_idx, mod_idx]
    rownames(tom_sub) <- colnames(tom_sub) <- mod_genes
    edges_mod <- as.data.frame(as.table(tom_sub)) %>%
        filter(Freq > threshold & Var1 != Var2) %>%
        rename(fromNode = Var1, toNode = Var2, weight = Freq)
    if(nrow(edges_mod) > 0) {
        edges_mod <- edges_mod[as.character(edges_mod$fromNode) < as.character(edges_mod$toNode), ]
        all_edges_list[[mod]] <- edges_mod
        write.table(edges_mod, paste0(output_dir, mod, "_EDGE.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    write.table(modules_df %>% filter(module == mod), paste0(output_dir, mod, "_NODE.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}
exports_dir <- "EXPORTS/"
if(!dir.exists(exports_dir)) dir.create(exports_dir, recursive = TRUE)
edge_filename <- paste0(exports_dir, project_name, "_FULL_EDGES.txt")
node_filename <- paste0(exports_dir, project_name, "_FULL_NODES.txt")
cat(paste0("Exporting Complete Edges and Nodes files for visualization: ", exports_dir, "\n"))
write.table(bind_rows(all_edges_list), edge_filename, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(modules_df, node_filename, sep = "\t", quote = FALSE, row.names = FALSE)
EOF

# Dynamic analysis selection
case $ANALYSIS in
    0)
        echo "==========================================="
        echo "Executing DATA CONVERSION (GEM/H5AD to GEF)..."
        echo "==========================================="
        singularity exec -B /g/data,/scratch $IMAGE python3 bin/SCRIPT_CONVERTER.py
        EXIT_CODE=$?
        
        if [ $EXIT_CODE -eq 0 ]; then
            echo "Conversion finished successfully."
        else
            echo "Conversion failed. Check logs."
        fi
        ;;
    1)
        echo "==========================================="
        echo "Executing PRIMARY ANALYSIS..."
        echo "==========================================="
        singularity exec -B /g/data,/scratch $IMAGE python3 bin/SCRIPT_PRIMARY_ANALYSIS.py
        EXIT_CODE=$?
        ;;
    2)
        echo "==========================================="
        echo "Executing SECONDARY ANALYSIS..."
        echo "==========================================="
        
        # Find and decompress the latest results folder
        LATEST_COMPRESSED=$(ls -t RESULTS_*.tar.gz 2>/dev/null | head -1)
        
        if [ -z "$LATEST_COMPRESSED" ]; then
            echo "ERROR: No compressed results folder found (RESULTS_*.tar.gz)"
            echo "Please run PRIMARY ANALYSIS (ANALYSIS=1) first"
            exit 1
        fi
        
        echo "Found compressed results: $LATEST_COMPRESSED"
        echo "Decompressing..."
        
        tar -xzf "$LATEST_COMPRESSED"
        
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to decompress $LATEST_COMPRESSED"
            exit 1
        fi
        
        # Get the decompressed folder name (remove .tar.gz extension)
        RESULTS_FOLDER="${LATEST_COMPRESSED%.tar.gz}"
        echo "Decompressed to: $RESULTS_FOLDER"
        echo "==========================================="
        
        # Execute secondary analysis
        singularity exec -B /g/data,/scratch $IMAGE python3 bin/SCRIPT_SECONDARY_ANALYSIS.py
        EXIT_CODE=$?
        ;;
     3)
        echo "==========================================="
        echo "Executing NETWORK ANALYSIS (hdWGCNA)"
        echo "==========================================="

        # Find and decompress the latest results folder
        LATEST_COMPRESSED=$(ls -t RESULTS_*.tar.gz 2>/dev/null | head -1)
        if [ -z "$LATEST_COMPRESSED" ]; then
            echo "ERROR: No compressed results folder found (RESULTS_*.tar.gz)"
            exit 1
        fi

        echo "Found compressed results: $LATEST_COMPRESSED"
        tar -xzf "$LATEST_COMPRESSED"
        RESULTS_FOLDER="${LATEST_COMPRESSED%.tar.gz}"
        
        # Enter folder and prepare environment
        cd "$RESULTS_FOLDER" || exit 1
        mkdir -p NETWORK
        
        echo "[$(date +%H:%M:%S)] Loading Environment for R..."
        module load R/4.3.1
        export R_LIBS_USER=~/.local/R/library
        ulimit -Sn 4000

        # Execute R script (pointing back to bin/)
        echo "[$(date +%H:%M:%S)] Running Rscript..."
        Rscript ../bin/SCRIPT_NETWORK_ANALYSIS.r 2>&1 | tee LOGS/NETWORK_ANALYSIS.log
        
        EXIT_CODE=${PIPESTATUS[0]}
        
        # 5. Return to root and cleanup
        cd ..
        
        if [ $EXIT_CODE -eq 0 ]; then
            echo "SUCCESS: Network analysis finished. Results in $RESULTS_FOLDER/NETWORK"
        else
            echo "ERROR: Network analysis failed."
        fi
        ;;
esac

# Check analysis execution status
if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "  Analysis $ANALYSIS completed successfully!"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "  ERROR: Analysis $ANALYSIS failed (exit code: $EXIT_CODE)"
    echo "=========================================="
    exit $EXIT_CODE
fi

echo "==========================================="
echo "Final system status:"
free -h
echo "==========================================="

# Finish: Move logs and compress results
if [ $EXIT_CODE -eq 0 ] && [ "$ANALYSIS" -ne 0 ]; then
    echo ""
    echo "==========================================="
    echo "ANALYSIS COMPLETED SUCCESSFULLY!"
    date
    
    if [ $ANALYSIS -eq 2 ]; then
        LATEST_RESULTS="$RESULTS_FOLDER"
    else
        LATEST_RESULTS=$(ls -td RESULTS_* 2>/dev/null | grep -v ".tar.gz" | head -1)
    fi
    
    if [ -n "$LATEST_RESULTS" ] && [ -d "$LATEST_RESULTS" ]; then
        echo "MOVING LOGS TO RESULTS FOLDER ($LATEST_RESULTS)"
        echo "==========================================="
        cp "SPATIAL_ANALYSYS_${$PBS_JOBID}.out" "$LATEST_RESULTS/" 2>/dev/null
        cp "SPATIAL_ANALYSIS_${$PBS_JOBID}.err" "$LATEST_RESULTS/" 2>/dev/null
        
        # Remove old compressed version if doing secondary analysis
        if [ $ANALYSIS -eq 2 ] && [ -f "${LATEST_RESULTS}.tar.gz" ]; then
            echo "Removing old compressed version..."
            rm -f "${LATEST_RESULTS}.tar.gz"
        fi
        
        echo "COMPRESSING RESULTS"
        tar -czf "${LATEST_RESULTS}.tar.gz" "$LATEST_RESULTS"
        
        if [ -f "${LATEST_RESULTS}.tar.gz" ]; then
            echo "COMPRESSION SUCCESSFUL. REMOVING FOLDER"
            rm -rf "$LATEST_RESULTS"
            (sleep 30 && rm -f "SPATIAL_ANALYSYS_${JOB_ID}.out" "SPATIAL_ANALYSIS_${JOB_ID}.err" 2>/dev/null) &
            echo "COMPRESSION FINISHED AND FOLDER CLEANED: ${LATEST_RESULTS}.tar.gz"
        else
            echo "ERROR: COMPRESSION FAILED."
        fi
    else
        echo "WARNING: No results folder found to compress"
    fi
    echo "==========================================="
    echo ""
else
    if [ "$ANALYSIS" -eq 0 ]; then
        echo "Conversion task finished. Results kept uncompressed for further analysis."
    else
        echo ""
        echo "==========================================="
        echo "ANALYSIS FAILED!"
        date
        echo "==========================================="
        echo "Check error logs for details"
        echo "==========================================="
        exit 1
    fi
fi

