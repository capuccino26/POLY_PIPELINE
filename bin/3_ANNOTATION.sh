#!/bin/bash
#$ -V
#$ -cwd
#$ -pe smp 16
#$ -l h_vmem=12G
#$ -N STEREOPY_ANNOTATION_ANALYSIS
#$ -o stereopy_annotation_$JOB_ID.out
#$ -e stereopy_annotation_$JOB_ID.err

echo "==========================================="
echo "STEREOPY ANNOTATION & EXPORT ANALYSIS"
echo "Cluster annotation and final data export"
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

# Load environment with explicit paths
echo "Loading miniconda3 module..."
module load miniconda3

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
    echo "   qsub -v ST_PYTHON=/home/user/.conda/envs/st/bin/python,MIN_COUNTS=50,MIN_GENES=5,PCT_COUNTS_MT=100,N_PCS=30 bin/3_ANNOTATION.sh" >&2
    
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
import anndata
print(f'AnnData version: {anndata.__version__}')
import scanpy as sc
print(f'Scanpy version: {sc.__version__}')
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
echo "  This step can be inproved after first run. Check the Elbow Plot (plots/qc/pca_elbow_enhanced.png) and insert the value of the elbow as N_PCS"
echo "You can alter the parameters inline:"
echo "  qsub -v ST_PYTHON="/home/user/.conda/envs/st/bin/python",MIN_COUNTS=50,MIN_GENES=5,PCT_COUNTS_MT=100,N_PCS=30 bin/3_ANNOTATION.sh"
echo ""

# Create the annotation analysis script
echo "Creating annotation analysis script"
cat > bin/stereopy_annotation_analysis.py << EOF
#!/usr/bin/env python3

# Import dependencies
import stereo as st
import warnings
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import psutil
import gc
from datetime import datetime
import sys
import json
from collections import defaultdict
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from scipy.sparse import issparse
import glob
import csv

# Set matplotlib backend for cluster
import matplotlib
matplotlib.use('Agg')

warnings.filterwarnings('ignore')

# Function for memory usage log
def log_memory_usage(step_name=""):
    try:
        memory = psutil.virtual_memory()
        swap = psutil.swap_memory()
        print(f"[MEMORY {step_name}] RAM: {memory.percent:.1f}% ({memory.used/1e9:.1f}GB/{memory.total/1e9:.1f}GB) | Swap: {swap.percent:.1f}%")
    except:
        print(f"[MEMORY {step_name}] Unable to get memory info")

# Function for correct logging
def log_step(message):
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()

# Function for checkpoints
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

# Function for analysis of genes of interest
def create_custom_gene_markers():
    file_path = 'INPUT/interest_genes.txt'
    gene_markers = {}

    if not os.path.exists(file_path):
        print(f"Warning: File of genes of interest not found in {file_path}. Proceeding without genes of interest.")
        return gene_markers

    try:
        with open(file_path, mode='r', newline='', encoding='utf-8') as file:
            reader = csv.reader(file)
            
            try:
                next(reader) 
            except StopIteration:
                print(f"Warning: The file {file_path} is empty. Proceeding without genes of interest")
                return gene_markers

            for row in reader:
                if not any(row):
                    continue

                marker_name = row[0].strip()
                loc_ids = [loc_id.strip() for loc_id in row[1:] if loc_id.strip()]
                
                if marker_name and loc_ids:
                    gene_markers[marker_name] = loc_ids
                elif marker_name and not loc_ids:
                    print(f"Warning: The gene '{marker_name}' was ignored because it has no gene IDs.")

    except Exception as e:
        print(f"ERROR: Processing {file_path} failed. Proceeding without genes of interest. Error: {e}")
        return {}
        
    return gene_markers

# Function for analyzing custom genes
def apply_gene_interest_annotation(data, gene_markers, cluster_res_key='louvain', threshold=1.2):
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

# Function for cluster visualization
def create_dynamic_gene_interest_visualization(data, annotation_details, gene_markers, output_dir):
    import matplotlib.pyplot as plt
    import numpy as np
    
    try:
        log_step("Creating dynamic visualization for gene interest clusters")
        
        if not hasattr(data.cells, 'spatial') or 'anno_custom_genes' not in data.cells:
            log_step("Warning: Missing spatial coordinates or custom gene annotations")
            return
        
        spatial_coords = data.cells.spatial
        annotations = data.cells['anno_custom_genes']
        
        # Dynamically identify enriched categories
        enriched_categories = []
        for category in gene_markers.keys():
            if gene_markers[category]:
                category_mask = annotations.str.contains(f'{category}_enriched', na=False)
                if category_mask.sum() > 0:
                    enriched_categories.append(category)
        
        log_step(f"Found {len(enriched_categories)} enriched categories: {enriched_categories}")
        
        if not enriched_categories:
            log_step("No enriched categories found for visualization")
            return
        
        # Generate colors dynamically
        colors = plt.cm.tab10(np.linspace(0, 1, len(enriched_categories)))
        
        # Create main visualization
        n_plots = len(enriched_categories) + 1  # +1 for combined plot
        n_cols = min(3, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
        if n_plots == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        axes = axes.flatten()
        
        # Combined plot (first subplot)
        ax = axes[0]
        
        # Background (non-enriched) in light gray
        enriched_mask = annotations.str.contains('enriched', na=False)
        non_enriched_mask = ~enriched_mask
        
        if non_enriched_mask.sum() > 0:
            ax.scatter(spatial_coords[non_enriched_mask, 0], 
                      spatial_coords[non_enriched_mask, 1],
                      c='lightgray', s=1, alpha=0.3, label='Other clusters')
        
        # Plot each enriched category with different color
        for i, category in enumerate(enriched_categories):
            category_mask = annotations.str.contains(f'{category}_enriched', na=False)
            if category_mask.sum() > 0:
                ax.scatter(spatial_coords[category_mask, 0], 
                          spatial_coords[category_mask, 1],
                          c=colors[i], s=2, alpha=0.8, 
                          label=f'{category.upper()} enriched')
        
        ax.set_title('All Gene Categories: Enriched Clusters', fontsize=14, fontweight='bold')
        ax.set_xlabel('Spatial X (μm)')
        ax.set_ylabel('Spatial Y (μm)')
        ax.legend(markerscale=3, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        
        # Individual plots for each category
        for i, category in enumerate(enriched_categories):
            ax = axes[i + 1]
            
            # Light gray background
            ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1], 
                      c='lightgray', s=0.5, alpha=0.3)
            
            # Highlight this category
            category_mask = annotations.str.contains(f'{category}_enriched', na=False)
            if category_mask.sum() > 0:
                ax.scatter(spatial_coords[category_mask, 0], 
                          spatial_coords[category_mask, 1],
                          c=colors[i], s=2, alpha=0.9)
                
                cell_count = category_mask.sum()
                total_cells = len(annotations)
                percentage = (cell_count / total_cells) * 100
                
                ax.set_title(f'{category.upper()} Enriched Only\n({cell_count:,} cells, {percentage:.1f}%)', 
                           fontweight='bold')
            else:
                ax.set_title(f'{category.upper()} Enriched\n(No cells found)', 
                           fontweight='bold')
            
            ax.set_xlabel('Spatial X')
            ax.set_ylabel('Spatial Y')
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for i in range(n_plots, len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        # Save the main plot
        main_plot_file = os.path.join(output_dir, 'plots', 'annotation', 'gene_interest_dynamic_spatial.png')
        plt.savefig(main_plot_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        # Create statistics plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Bar plot of cell counts
        category_names = []
        cell_counts = []
        category_colors = []
        
        for i, category in enumerate(enriched_categories):
            category_mask = annotations.str.contains(f'{category}_enriched', na=False)
            category_names.append(f'{category.upper()}\nEnriched')
            cell_counts.append(category_mask.sum())
            category_colors.append(colors[i])
        
        # Add "Other" category
        category_names.append('Other\nClusters')
        cell_counts.append(non_enriched_mask.sum())
        category_colors.append('lightgray')
        
        bars = ax1.bar(category_names, cell_counts, color=category_colors, alpha=0.8)
        ax1.set_title('Cell Distribution by Gene Interest Category', fontweight='bold')
        ax1.set_ylabel('Number of Cells')
        
        # Add value labels on bars
        total_cells = sum(cell_counts)
        for bar, count in zip(bars, cell_counts):
            ax1.text(bar.get_x() + bar.get_width()/2, 
                    bar.get_height() + max(cell_counts)*0.01,
                    f'{count:,}\n({count/total_cells*100:.1f}%)', 
                    ha='center', va='bottom', fontweight='bold')
        
        # Pie chart
        pie_colors = category_colors
        wedges, texts, autotexts = ax2.pie(cell_counts, labels=category_names, 
                                          colors=pie_colors, autopct='%1.1f%%',
                                          startangle=90)
        ax2.set_title('Proportion of Cells by Gene Interest', fontweight='bold')
        
        plt.tight_layout()
        
        # Save statistics plot
        stats_plot_file = os.path.join(output_dir, 'plots', 'annotation', 'gene_interest_statistics.png')
        plt.savefig(stats_plot_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        # Create publication-ready plot (clean version)
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        
        # Only show enriched clusters, no background
        for i, category in enumerate(enriched_categories):
            category_mask = annotations.str.contains(f'{category}_enriched', na=False)
            if category_mask.sum() > 0:
                ax.scatter(spatial_coords[category_mask, 0], 
                          spatial_coords[category_mask, 1],
                          c=colors[i], s=3, alpha=0.9, 
                          label=f'{category.upper()} enriched', rasterized=True)
        
        # Clean styling
        ax.set_aspect('equal')
        ax.set_xlabel('Spatial X (μm)', fontsize=12)
        ax.set_ylabel('Spatial Y (μm)', fontsize=12)
        ax.legend(frameon=False, fontsize=12, markerscale=2)
        
        # Remove spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5)
        
        plt.tight_layout()
        
        # Save publication version
        pub_plot_file = os.path.join(output_dir, 'plots', 'annotation', 'gene_interest_publication_ready.png')
        plt.savefig(pub_plot_file, dpi=300, bbox_inches='tight', facecolor='white')
        
        pub_plot_pdf = os.path.join(output_dir, 'plots', 'annotation', 'gene_interest_publication_ready.pdf')
        plt.savefig(pub_plot_pdf, bbox_inches='tight', facecolor='white')
        plt.close()
        
        # Log results
        log_step(f"Dynamic visualizations saved:")
        log_step(f"  - Main plot: gene_interest_dynamic_spatial.png")
        log_step(f"  - Statistics: gene_interest_statistics.png")
        log_step(f"  - Publication: gene_interest_publication_ready.png/.pdf")
        
        for category in enriched_categories:
            category_mask = annotations.str.contains(f'{category}_enriched', na=False)
            count = category_mask.sum()
            percentage = (count / len(annotations)) * 100
            log_step(f"  - {category.upper()}: {count:,} cells ({percentage:.1f}%)")
        
    except Exception as e:
        log_step(f"Warning: Error creating dynamic visualization: {e}")
        import traceback
        log_step(f"Traceback: {traceback.format_exc()}")

#Function for individual clusters
def create_individual_gene_plots(data, annotation_details, gene_markers, output_dir):
    try:
        log_step("Creating individual plots for each gene interest category")
        # Get spatial coordinates
        if not hasattr(data, 'position') or data.position is None or data.position.size == 0:
            log_step("Warning: Missing spatial coordinates. Cannot create plots.")
            return
        spatial_coords = data.position

        # Get annotations
        anno_key = 'anno_custom_genes'
        if not hasattr(data, 'cells') or anno_key not in data.cells:
            log_step(f"Warning: Missing annotation '{anno_key}'. Cannot create plots.")
            return
        annotations = data.cells[anno_key]
        total_cells = len(annotations)
        
        # Identify enriched categories
        enriched_categories = []
        for category in gene_markers.keys():
            category_mask = annotations.str.contains(f'{category}_enriched', na=False)
            if category_mask.sum() > 0:
                enriched_categories.append(category)

        if not enriched_categories:
            log_step("No enriched categories found to plot.")
            return

        # Generate colors
        colors = plt.cm.tab10(np.linspace(0, 1, len(enriched_categories)))

        # Create plots for each category
        for i, category in enumerate(enriched_categories):
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))

            non_enriched_mask = ~annotations.str.contains('enriched', na=False)
            if non_enriched_mask.sum() > 0:
                ax.scatter(spatial_coords[non_enriched_mask, 0],
                           spatial_coords[non_enriched_mask, 1],
                           c='lightgray', s=0.5, alpha=0.3, label='Other clusters')

            category_mask = annotations.str.contains(f'{category}_enriched', na=False)
            cell_count = category_mask.sum()
            percentage = (cell_count / total_cells) * 100

            ax.scatter(spatial_coords[category_mask, 0],
                       spatial_coords[category_mask, 1],
                       c=[colors[i]], s=3, alpha=0.9,
                       label=f'{category.upper()} enriched\n({cell_count:,} cells, {percentage:.1f}%)')

            ax.set_title(f'Spatial Distribution of {category.upper()} Enriched Clusters',
                         fontweight='bold', fontsize=14)
            ax.set_xlabel('Spatial X (μm)')
            ax.set_ylabel('Spatial Y (μm)')
            ax.set_aspect('equal', adjustable='box')
            ax.legend(markerscale=3, loc='upper right')
            ax.grid(True, alpha=0.3)
            plt.tight_layout()

            plot_dir = os.path.join(output_dir, 'plots', 'annotation')
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir, exist_ok=True)
            plot_file = os.path.join(plot_dir, f'spatial_enrichment_{category}.png')
            
            plt.savefig(plot_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()

            log_step(f"Plot for '{category.upper()}' saved: {plot_file}")
            print(f"File size: {os.path.getsize(plot_file) / 1e6:.1f} MB")

    except Exception as e:
        log_step(f"Warning: Error creating individual gene plots: {e}")
        import traceback
        log_step(f"Traceback: {traceback.format_exc()}")

# Function for individual clustering
def create_individual_cluster_plots_per_annotation(data, annotation_templates, output_dir):
    try:
        log_step("Creating individual cluster plots for each annotation category...")
        
        # Get spatial coordinates
        if hasattr(data, 'position') and data.position is not None:
            spatial_coords = data.position
            x_coords = spatial_coords[:, 0]
            y_coords = spatial_coords[:, 1]
        elif hasattr(data.cells, 'spatial'):
            spatial_coords = data.cells.spatial
            x_coords = spatial_coords[:, 0]
            y_coords = spatial_coords[:, 1]
        elif 'x' in data.cells.columns and 'y' in data.cells.columns:
            x_coords = data.cells['x'].values
            y_coords = data.cells['y'].values
        else:
            log_step("ERROR: Could not find spatial coordinates")
            return
        
        # Main individual clusters directory
        individual_clusters_base_dir = os.path.join(output_dir, 'plots', 'individual_clusters')
        os.makedirs(individual_clusters_base_dir, exist_ok=True)
        
        # Process each annotation template
        for template_name, annotation_dict in annotation_templates.items():
            if not annotation_dict:
                log_step(f"Skipping {template_name} - no annotations available")
                continue
                
            res_key = f'anno_{template_name}'
            if res_key not in data.cells:
                log_step(f"Skipping {template_name} - annotation not applied to data")
                continue
            
            log_step(f"Processing individual clusters for {template_name} annotation...")
            
            # Create subdirectory
            annotation_clusters_dir = os.path.join(individual_clusters_base_dir, template_name)
            os.makedirs(annotation_clusters_dir, exist_ok=True)
            
            # Get the annotation labels for all cells
            cell_annotations = data.cells[res_key]
            unique_annotations = sorted(cell_annotations.unique())
            n_annotations = len(unique_annotations)
            
            log_step(f"Found {n_annotations} unique annotations for {template_name}")
            
            # Create individual plots for each annotation
            for i, annotation_label in enumerate(unique_annotations):
                log_step(f"Processing {template_name} annotation '{annotation_label}' ({i+1}/{n_annotations})...")
                
                try:
                    # Create a copy of annotations for highlighting
                    annotation_highlight = cell_annotations.copy()
                    
                    # Handle Categorical data properly
                    if hasattr(annotation_highlight, 'cat'):
                        if 'Other' not in annotation_highlight.cat.categories:
                            annotation_highlight = annotation_highlight.cat.add_categories(['Other'])
                    
                    # Set all other annotations to 'Other'
                    annotation_highlight[annotation_highlight != annotation_label] = 'Other'
                    
                    # Create the plot
                    fig, ax = plt.subplots(figsize=(12, 10))
                    
                    # Plot background (other annotations) in gray
                    other_mask = annotation_highlight == 'Other'
                    if other_mask.any():
                        ax.scatter(x_coords[other_mask], y_coords[other_mask], 
                                  c='lightgray', s=0.8, alpha=0.4, rasterized=True, 
                                  label='Other annotations')
                    
                    # Plot the highlighted annotation in color
                    target_mask = annotation_highlight == annotation_label
                    if target_mask.any():
                        # Use a distinct color - cycle through color palette
                        colors = plt.cm.tab20(i % 20)
                        ax.scatter(x_coords[target_mask], y_coords[target_mask], 
                                  c=[colors], s=2.5, alpha=0.9, rasterized=True,
                                  label=f'{annotation_label}')
                    
                    # Formatting
                    ax.set_xlabel('Spatial X (μm)', fontsize=12)
                    ax.set_ylabel('Spatial Y (μm)', fontsize=12)
                    ax.set_title(f'{template_name.capitalize()} Annotation: {annotation_label}\nSpatial Distribution', 
                               fontsize=14, fontweight='bold')
                    
                    # Add legend
                    ax.legend(loc='upper right', fontsize=10, frameon=True, 
                             fancybox=True, shadow=True, markerscale=3)
                    
                    # Add scale bar (assuming coordinates are in micrometers)
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
                    
                    # Add statistics as text
                    n_cells_annotation = target_mask.sum()
                    total_cells = len(annotation_highlight)
                    percentage = (n_cells_annotation / total_cells) * 100
                    
                    stats_text = f'Cells with this annotation: {n_cells_annotation:,}\nPercentage: {percentage:.1f}%'
                    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                           verticalalignment='top', fontsize=11, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))
                    
                    # Set equal aspect ratio and adjust layout
                    ax.set_aspect('equal')
                    ax.grid(True, alpha=0.3)
                    plt.tight_layout()
                    
                    # Save the plot - create safe filename
                    safe_annotation_name = str(annotation_label).replace('/', '_').replace(' ', '_')
                    output_file = os.path.join(annotation_clusters_dir, f'{safe_annotation_name}_spatial.png')
                    plt.savefig(output_file, dpi=300, bbox_inches='tight', 
                               facecolor='white', edgecolor='none')
                    plt.close()
                    
                    log_step(f"Saved individual plot: {safe_annotation_name}_spatial.png")
                    
                except Exception as e:
                    log_step(f"ERROR generating plot for annotation '{annotation_label}': {e}")
                    continue
            
            # Generate HTML summary for this annotation type
            log_step(f"Generating HTML summary for {template_name} individual annotations...")
            html_content = f"""
            <!DOCTYPE html>
            <html>
            <head>
                <title>{template_name.capitalize()} Annotation Analysis</title>
                <style>
                    body {{ font-family: 'Arial', sans-serif; margin: 20px; background-color: #f8f9fa; }}
                    .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                              color: white; padding: 20px; border-radius: 10px; margin-bottom: 20px; }}
                    .annotation-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(450px, 1fr)); 
                                       gap: 20px; padding: 20px 0; }}
                    .annotation-item {{ 
                        background: white; border: 2px solid #e9ecef; padding: 15px; 
                        border-radius: 10px; box-shadow: 0 4px 6px rgba(0,0,0,0.1);
                        transition: transform 0.2s ease, box-shadow 0.2s ease;
                    }}
                    .annotation-item:hover {{ 
                        transform: translateY(-2px); 
                        box-shadow: 0 6px 12px rgba(0,0,0,0.15);
                    }}
                    .annotation-item img {{ width: 100%; height: auto; border-radius: 5px; }}
                    .annotation-title {{ font-weight: bold; margin-bottom: 10px; font-size: 16px; 
                                        color: #495057; border-bottom: 2px solid #007bff; padding-bottom: 5px; }}
                    h1 {{ color: white; text-align: center; margin: 0; font-size: 28px; }}
                    .summary {{ background: white; padding: 20px; border-radius: 10px; 
                               margin-bottom: 20px; border-left: 5px solid #007bff; }}
                    .summary h3 {{ color: #007bff; margin-top: 0; }}
                    .stats {{ display: flex; justify-content: space-around; margin: 10px 0; }}
                    .stat-item {{ text-align: center; }}
                    .stat-number {{ font-size: 24px; font-weight: bold; color: #007bff; }}
                    .stat-label {{ font-size: 14px; color: #6c757d; }}
                </style>
            </head>
            <body>
                <div class="header">
                    <h1>{template_name.capitalize()} Annotation - Individual Analysis</h1>
                </div>
                <div class="summary">
                    <h3>Analysis Summary</h3>
                    <div class="stats">
                        <div class="stat-item">
                            <div class="stat-number">{n_annotations}</div>
                            <div class="stat-label">Unique Annotations</div>
                        </div>
                        <div class="stat-item">
                            <div class="stat-number">{len(cell_annotations):,}</div>
                            <div class="stat-label">Total Cells</div>
                        </div>
                        <div class="stat-item">
                            <div class="stat-number">{template_name.capitalize()}</div>
                            <div class="stat-label">Annotation Type</div>
                        </div>
                    </div>
                    <p><strong>Analysis method:</strong> Individual annotation highlighting (target annotation in color, others in gray)</p>
                    <p><strong>Visualization:</strong> Spatial distribution with statistical information</p>
                </div>
                <div class="annotation-grid">
            """
            
            for annotation_label in unique_annotations:
                safe_annotation_name = str(annotation_label).replace('/', '_').replace(' ', '_')
                target_mask = cell_annotations == annotation_label
                count = target_mask.sum()
                percentage = (count / len(cell_annotations)) * 100
                
                html_content += f"""
                    <div class="annotation-item">
                        <div class="annotation-title">{annotation_label}</div>
                        <img src="{safe_annotation_name}_spatial.png" alt="{annotation_label}">
                        <div style="margin-top: 10px; font-size: 14px; color: #6c757d;">
                            <strong>Cells:</strong> {count:,} ({percentage:.1f}%)
                        </div>
                    </div>
                """
            
            html_content += """
                </div>
                <div style="text-align: center; margin-top: 40px; padding: 20px; background: white; border-radius: 10px;">
                    <p style="color: #6c757d; font-style: italic;">
                        Individual annotation analysis completed. Each plot shows the spatial distribution 
                        of cells with the specific annotation highlighted.
                    </p>
                </div>
            </body>
            </html>
            """
            
            html_file = os.path.join(annotation_clusters_dir, f'{template_name}_individual_annotations_summary.html')
            with open(html_file, 'w') as f:
                f.write(html_content)
            
            log_step(f"Generated HTML summary for {template_name}: {template_name}_individual_annotations_summary.html")
            log_step(f"Individual plots saved to: {annotation_clusters_dir}")
        
        # Create master HTML index for all annotation types
        log_step("Creating master HTML index for all annotation types...")
        master_html = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Complete Annotation Analysis - Individual Clusters</title>
            <style>
                body { font-family: 'Arial', sans-serif; margin: 20px; background-color: #f8f9fa; }
                .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                          color: white; padding: 30px; border-radius: 15px; margin-bottom: 30px; text-align: center; }
                .annotation-type-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); 
                                       gap: 25px; margin: 20px 0; }
                .annotation-type-card { 
                    background: white; border: 2px solid #e9ecef; padding: 25px; 
                    border-radius: 15px; box-shadow: 0 6px 12px rgba(0,0,0,0.1);
                    transition: all 0.3s ease;
                }
                .annotation-type-card:hover { 
                    transform: translateY(-5px); 
                    box-shadow: 0 10px 20px rgba(0,0,0,0.15);
                    border-color: #007bff;
                }
                .card-title { font-size: 24px; font-weight: bold; color: #007bff; margin-bottom: 15px; }
                .card-description { color: #6c757d; margin-bottom: 20px; line-height: 1.6; }
                .view-button { 
                    background: linear-gradient(135deg, #007bff 0%, #0056b3 100%);
                    color: white; padding: 12px 25px; border: none; border-radius: 25px;
                    font-weight: bold; cursor: pointer; text-decoration: none;
                    display: inline-block; transition: all 0.3s ease;
                }
                .view-button:hover { 
                    background: linear-gradient(135deg, #0056b3 0%, #004085 100%);
                    transform: scale(1.05);
                }
                h1 { margin: 0; font-size: 32px; }
                .subtitle { margin: 10px 0 0 0; font-size: 18px; opacity: 0.9; }
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Complete Annotation Analysis</h1>
                <p class="subtitle">Individual Cluster Visualization by Annotation Type</p>
            </div>
            
            <div class="annotation-type-grid">
        """
        
        # Add cards for each annotation type that has data
        annotation_descriptions = {
            'custom_genes': 'Gene interest-based annotations using predefined marker genes and expression thresholds',
            'alphabetical': 'Simple alphabetical labeling of clusters for easy reference and communication',
            'numerical': 'Numerical cluster naming with size-based descriptive labels',
            'biological': 'Biologically-informed annotations based on marker gene expression profiles'
        }
        
        for template_name in annotation_templates.keys():
            if annotation_templates[template_name] and f'anno_{template_name}' in data.cells:
                description = annotation_descriptions.get(template_name, 'Custom annotation strategy')
                n_unique = data.cells[f'anno_{template_name}'].nunique()
                
                master_html += f"""
                    <div class="annotation-type-card">
                        <div class="card-title">{template_name.capitalize()} Annotations</div>
                        <div class="card-description">
                            {description}
                            <br><br>
                            <strong>Unique annotations:</strong> {n_unique}
                        </div>
                        <a href="{template_name}/{template_name}_individual_annotations_summary.html" class="view-button">
                            View Individual Plots
                        </a>
                    </div>
                """
        
        master_html += """
            </div>
            
            <div style="text-align: center; margin-top: 40px; padding: 30px; background: white; border-radius: 15px;">
                <h3 style="color: #007bff;">Analysis Complete</h3>
                <p style="color: #6c757d; line-height: 1.6;">
                    Individual cluster visualization has been generated for all annotation strategies. 
                    Each annotation type has its own dedicated analysis page with spatial distribution plots 
                    and statistical summaries.
                </p>
            </div>
        </body>
        </html>
        """
        
        master_index_file = os.path.join(individual_clusters_base_dir, 'index.html')
        with open(master_index_file, 'w') as f:
            f.write(master_html)
        
        log_step(f"Master HTML index created: {master_index_file}")
        log_step(f"Individual cluster analysis completed for all annotation types")
        log_step(f"Results organized in: {individual_clusters_base_dir}")
        
        return True
        
    except Exception as e:
        log_step(f"ERROR in individual cluster plots creation: {e}")
        import traceback
        log_step(f"Full traceback: {traceback.format_exc()}")
        return False

# Function for direct gene visualization
def create_direct_gene_visualization():
    try:
        log_step("Extracting data directly for gene-level visualization...")
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
                       c='lightgray', s=0.5, alpha=0.3, label='Other clusters')
        
        for i, category in enumerate(enriched_categories_list):
            mask = cell_categories == category
            count = mask.sum()
            ax.scatter(spatial_coords[mask, 0], spatial_coords[mask, 1],
                       c=[colors[i]], s=3, alpha=0.8,
                       label=f'{category.upper()}: {count:,} cells')
        
        ax.set_title('Gene Interest Spatial Distribution\n(Dynamic Gene Categories)', fontweight='bold', fontsize=14)
        ax.set_xlabel('Spatial X')
        ax.set_ylabel('Spatial Y')
        ax.legend(markerscale=3, fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        plt.tight_layout(rect=[0, 0, 0.8, 1])
        main_plot = os.path.join(annotation_output_dir, 'plots', 'gene_interest_direct_spatial_combined.png')
        plt.savefig(main_plot, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        log_step(f"Main combined plot saved: {main_plot}")

        log_step("Creating individual plots for each gene of interest")
        individual_gene_plot_dir = os.path.join(annotation_output_dir, 'plots', 'gene_expression_individual')
        os.makedirs(individual_gene_plot_dir, exist_ok=True)
        
        all_genes_of_interest = []
        for category, genes in gene_markers.items():
            all_genes_of_interest.extend(genes)
        
        for gene_name in sorted(list(set(all_genes_of_interest))):
            if gene_name not in gene_names:
                log_step(f"Warning: Gene '{gene_name}' not found in expression data. Skipping individual plot.")
                continue
            
            gene_idx = gene_names.index(gene_name)
            gene_expression = exp_matrix[:, gene_idx]
            
            normalized_expression = gene_expression / (gene_expression.max() if gene_expression.max() > 0 else 1)
            
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            sc = ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                            c=gene_expression,
                            cmap='viridis',
                            s=3, alpha=0.8, vmin=0, vmax=gene_expression.max())
            
            ax.set_title(f'Spatial Expression of Gene: {gene_name}',
                         fontweight='bold', fontsize=14)
            ax.set_xlabel('Spatial X (μm)')
            ax.set_ylabel('Spatial Y (μm)')
            ax.set_aspect('equal', adjustable='box')
            
            cbar = fig.colorbar(sc, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
            cbar.set_label('Gene Expression Level', rotation=270, labelpad=15)
            
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            
            plot_file = os.path.join(individual_gene_plot_dir, f'spatial_expression_{gene_name}.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            
            log_step(f"Individual plot for gene '{gene_name}' saved: {plot_file}")
            print(f"File size: {os.path.getsize(plot_file) / 1e6:.1f} MB")
            
        report_file = os.path.join(annotation_output_dir, 'direct_gene_enrichment_report.txt')
        with open(report_file, 'w') as f:
            f.write("DIRECT GENE ENRICHMENT ANALYSIS REPORT\n")
            f.write("="*60 + "\n\n")
            f.write(f"Analysis date: {datetime.now()}\n")
            f.write(f"Total cells: {len(cell_categories):,}\n")
            f.write(f"Total clusters: {base_clusters.nunique()}\n")
            f.write(f"Enrichment threshold: 1.2\n\n")
            
            f.write("GENE CATEGORIES DEFINED:\n")
            for category, genes in gene_markers.items():
                f.write(f"\n{category.upper()}:\n")
                if genes:
                    for gene in genes:
                        f.write(f"  - {gene}\n")
                else:
                    f.write("  - (No genes defined)\n")
            
            f.write(f"\nENRICHMENT RESULTS BY CATEGORY:\n")
            f.write("-"*40 + "\n")
            
            for category in enriched_categories_list:
                count = cell_categories.value_counts().get(category, 0) # Obter count novamente
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
        
        log_step("Creating clean spatial visualization from direct analysis results...")
        
        spatial_coords = direct_results['spatial_coords']
        cell_categories = direct_results['cell_categories']
        enriched_clusters = direct_results['enriched_clusters']
        stereopy_spatial = None
        if hasattr(data.cells, 'obs'):
            obs_df = data.cells.obs
            numeric_cols = obs_df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) >= 2:
                stereopy_spatial = obs_df[numeric_cols[:2]].values
        print(f"Using direct analysis spatial coordinates: {spatial_coords.shape}")
        fig, ax = plt.subplots(1, 1, figsize=(14, 12))
        other_mask = cell_categories == 'Other'
        if other_mask.sum() > 0:
            ax.scatter(spatial_coords[other_mask, 0], spatial_coords[other_mask, 1],
                      c='whitesmoke', s=0.3, alpha=0.4, edgecolors='lightgray', 
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
                      c=category_colors[category], s=4, alpha=0.9,
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
        ax.set_aspect('equal', adjustable='box')
        plt.tight_layout()
        
        clean_plot = os.path.join(output_dir, 'plots', 'gene_interest_CLEAN_spatial.png')
        plt.savefig(clean_plot, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        log_step(f"Clean spatial plot saved: {clean_plot}")
        print(f"File size: {os.path.getsize(clean_plot) / 1e6:.1f} MB")
        
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
        ax.set_aspect('equal', adjustable='box')
        
        plt.tight_layout()
        
        genes_only_plot = os.path.join(output_dir, 'plots', 'gene_interest_ONLY_spatial.png')
        plt.savefig(genes_only_plot, dpi=300, bbox_inches='tight', facecolor='white')
        
        genes_only_pdf = os.path.join(output_dir, 'plots', 'gene_interest_ONLY_spatial.pdf')
        plt.savefig(genes_only_pdf, bbox_inches='tight', facecolor='white')
        plt.close()
        
        log_step(f"Genes-only plot saved: {genes_only_plot} and {genes_only_pdf}")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        if other_mask.sum() > 0:
            ax1.scatter(spatial_coords[other_mask, 0], spatial_coords[other_mask, 1],
                       c='lightgray', s=0.5, alpha=0.3, rasterized=True)
        
        for category in enriched_categories:
            mask = cell_categories == category
            ax1.scatter(spatial_coords[mask, 0], spatial_coords[mask, 1],
                       c=category_colors[category], s=3, alpha=0.9,
                       label=f'{category.upper()}', rasterized=True)
        
        ax1.set_title('With Context\n(All Clusters)', fontweight='bold', fontsize=14)
        ax1.set_xlabel('Spatial X (μm)')
        ax1.set_ylabel('Spatial Y (μm)')
        ax1.legend(frameon=False, fontsize=10, markerscale=2)
        ax1.grid(True, alpha=0.2)
        ax1.set_aspect('equal', adjustable='box')
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
        ax2.set_aspect('equal', adjustable='box')
        
        for ax in [ax1, ax2]:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        
        comparison_plot = os.path.join(output_dir, 'plots', 'gene_interest_COMPARISON_spatial.png')
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
        
        stats_plot = os.path.join(output_dir, 'plots', 'gene_interest_STATISTICS.png')
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

# INITIALIZATION
start_time = datetime.now()
log_step("Starting annotation and export analysis")
log_memory_usage("START")

# Check if previous analysis results exist
results_dir = 'results_ultimate'
if not os.path.exists(results_dir):
    log_step("ERROR: Previous analysis results not found. Please run bin/2_DOC_ANALYSIS.sh first")
    sys.exit(1)

# Set up directories for annotation analysis
annotation_output_dir = 'RESULTS/results_annotation'
directories = [
    annotation_output_dir,
    os.path.join(annotation_output_dir, 'plots'),
    os.path.join(annotation_output_dir, 'plots', 'annotation'),
    os.path.join(annotation_output_dir, 'annotations'),
    os.path.join(annotation_output_dir, 'annotations', 'automated'),
    os.path.join(annotation_output_dir, 'exports'),
    os.path.join(annotation_output_dir, 'exports', 'h5ad'),
    os.path.join(annotation_output_dir, 'exports', 'csv'),
    os.path.join(annotation_output_dir, 'exports', 'metadata'),
    os.path.join(annotation_output_dir, 'biological_interpretation'),
    os.path.join(annotation_output_dir, 'logs'),
    os.path.join(annotation_output_dir, 'checkpoints')
]

for directory in directories:
    os.makedirs(directory, exist_ok=True)

log_step(f"Created {len(directories)} folders")

# Filepaths
data_path = get_single_gef_file()
if not os.path.exists(data_path):
    log_step(f"ERROR: Data file not found at {data_path}")
    sys.exit(1)

# Filtering parameters from bash
MIN_COUNTS = $MIN_COUNTS
MIN_GENES = $MIN_GENES
PCT_COUNTS_MT = $PCT_COUNTS_MT
N_PCS = $N_PCS

log_step("Reloading data and restoring analysis state")
try:
    data = st.io.read_gef(file_path=data_path, bin_size=100)
    log_step(f"Data reloaded: {data}")
    
    log_step("Applying preprocessing pipeline")
    data.tl.filter_cells(min_counts=MIN_COUNTS, min_genes=MIN_GENES, pct_counts_mt=PCT_COUNTS_MT, inplace=True)
    data.tl.normalize_total(target_sum=10000)
    data.tl.log1p()
    data.tl.raw_checkpoint()
    
    data.tl.highly_variable_genes(min_mean=0.0125, max_mean=3, min_disp=0.5, 
                                  n_top_genes=2000, res_key='highly_variable_genes')
    data.tl.scale(max_value=10, zero_center=False)
    data.tl.pca(use_highly_genes=True, n_pcs=N_PCS, res_key='pca')
    data.tl.neighbors(pca_res_key='pca', n_pcs=N_PCS, res_key='neighbors')
    data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap')
    data.tl.louvain(neighbors_res_key='neighbors', res_key='louvain')
    
    data.tl.find_marker_genes(cluster_res_key='louvain', method='t_test', 
                             use_highly_genes=False, use_raw=True, 
                             res_key='marker_genes')
    
    log_step("Analysis state restored successfully")
    log_memory_usage("DATA_RESTORED")
    
except Exception as e:
    log_step(f"ERROR restoring analysis state: {e}")
    sys.exit(1)

# CLUSTER ANNOTATION

log_step("="*80)
log_step("CLUSTER ANNOTATION")
log_step("="*80)

annotation_log_file = os.path.join(annotation_output_dir, 'logs', 'annotation_log.txt')
with open(annotation_log_file, 'w') as f:
    f.write("CLUSTER ANNOTATION LOG\n")
    f.write("="*60 + "\n\n")
    f.write(f"Analysis start: {start_time}\n")
    f.write(f"Clustering method: Louvain\n")
    f.write(f"Total clusters: {len(data.cells['louvain'].unique()) if 'louvain' in data.cells else 'Unknown'}\n\n")

# Strategy 1: Apply custom gene interest annotation
log_step("Applying custom gene interest annotation")

custom_gene_markers = create_custom_gene_markers()
print(custom_gene_markers)
custom_annotations, custom_scores, annotation_details = apply_gene_interest_annotation(
    data, custom_gene_markers, threshold=1.2
)

# Create detailed report for custom gene annotation
custom_report_file = os.path.join(annotation_output_dir, 'biological_interpretation', 'custom_gene_interest_report.txt')
with open(custom_report_file, 'w') as f:
    f.write("CUSTOM GENES OF INTEREST ANNOTATION REPORT\n")
    f.write("="*70 + "\n\n")
    
    f.write("GENES OF INTEREST DEFINED:\n")
    for category, genes in custom_gene_markers.items():
        f.write(f"\n{category.upper()}:\n")
        if genes:
            for gene in genes:
                f.write(f"  - {gene}\n")
        else:
            f.write("  - (No genes defined yet)\n")
    
    f.write(f"\nANNOTATION RESULTS (threshold: 1.2):\n")
    f.write("-" * 50 + "\n")
    
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

# Strategy 2: Load marker genes for biological interpretation
marker_genes_file = os.path.join(results_dir, 'marker_genes', 'top50_per_cluster_for_viz.csv')
biological_annotations = {}

if os.path.exists(marker_genes_file):
    try:
        log_step("Loading marker genes for biological annotation...")
        marker_df = pd.read_csv(marker_genes_file)
        log_step(f"Loaded {len(marker_df)} marker gene records")
        
        # Create biological annotations based on marker genes
        biological_annotations = create_biological_annotation_dict(marker_df, n_top_genes=5)
        log_step(f"Created biological annotations for {len(biological_annotations)} clusters")
        
        # Save biological annotation mapping
        bio_anno_file = os.path.join(annotation_output_dir, 'annotations', 'automated', 'biological_annotations.json')
        with open(bio_anno_file, 'w') as f:
            json.dump(biological_annotations, f, indent=2)
        
        # Create detailed marker-based annotation report
        marker_report_file = os.path.join(annotation_output_dir, 'biological_interpretation', 'marker_based_annotation_report.txt')
        with open(marker_report_file, 'w') as f:
            f.write("MARKER-BASED BIOLOGICAL ANNOTATION REPORT\n")
            f.write("="*70 + "\n\n")
            
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
        log_step(f"Warning: Error loading marker genes for annotation: {e}")

# Strategy 3: Manual annotation
log_step("Creating manual annotation")

# Basic alphabetical annotation
alphabetical_annotation = {}
if 'louvain' in data.cells:
    unique_clusters = sorted(data.cells['louvain'].unique(), key=lambda x: int(x))
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    
    for i, cluster in enumerate(unique_clusters):
        if i < len(alphabet):
            alphabetical_annotation[str(cluster)] = alphabet[i]
        else:
            # Handle more than 26 clusters
            alphabetical_annotation[str(cluster)] = f"cluster_{cluster}"

# Numerical annotation with descriptive names
numerical_annotation = {}
if 'louvain' in data.cells:
    cluster_sizes = data.cells['louvain'].value_counts().sort_index()
    
    for cluster in cluster_sizes.index:
        size = cluster_sizes[cluster]
        if size > 1000:
            size_desc = "large"
        elif size > 100:
            size_desc = "medium"
        else:
            size_desc = "small"
        
        numerical_annotation[str(cluster)] = f"Cluster_{cluster}_{size_desc}"

# Compile all annotation templates
annotation_templates = {
    'custom_genes': custom_annotations,
    'alphabetical': alphabetical_annotation,
    'numerical': numerical_annotation,
    'biological': biological_annotations
}

# Save annotation templates
for template_name, template_dict in annotation_templates.items():
    if template_dict:  # Only save non-empty templates
        template_file = os.path.join(annotation_output_dir, 'annotations', 'automated', f'{template_name}_annotation.json')
        with open(template_file, 'w') as f:
            if template_name == 'custom_genes':
                # Save additional details for custom genes
                save_data = {
                    'annotations': template_dict,
                    'scores': custom_scores,
                    'details': annotation_details,
                    'gene_markers': custom_gene_markers,
                    'threshold_used': 1.2
                }
                json.dump(save_data, f, indent=2)
            else:
                json.dump(template_dict, f, indent=2)

log_step(f"Created {len([t for t in annotation_templates.values() if t])} annotation templates")

# Strategy 4: Apply annotations and create visualizations
for template_name, annotation_dict in annotation_templates.items():
    if not annotation_dict:
        continue
        
    try:
        log_step(f"Applying {template_name} annotation")
        
        # Apply annotation
        res_key = f'anno_{template_name}'
        data.tl.annotation(
            annotation_information=annotation_dict,
            cluster_res_key='louvain',
            res_key=res_key
        )
        
        # Create visualization
        plot_file = os.path.join(annotation_output_dir, 'plots', 'annotation', f'{template_name}_annotation_spatial.png')
        data.plt.cluster_scatter(res_key=res_key, out_path=plot_file)
        
        # Create UMAP visualization if available
        if 'umap' in data.tl.result:
            umap_plot_file = os.path.join(annotation_output_dir, 'plots', 'annotation', f'{template_name}_annotation_umap.png')
            data.plt.umap(res_key='umap', cluster_key=res_key, out_path=umap_plot_file)
        
        log_step(f"{template_name.capitalize()} annotation applied and visualized")
        
        # Save annotation assignments
        if res_key in data.cells:
            annotation_df = pd.DataFrame({
                'cell_name': data.cells.cell_name,
                'cluster_id': data.cells['louvain'],
                'annotation': data.cells[res_key]
            })
            
            # Add spatial coordinates if available
            if hasattr(data.cells, 'spatial'):
                spatial_coords = data.cells.spatial
                annotation_df['x_coord'] = spatial_coords[:, 0]
                annotation_df['y_coord'] = spatial_coords[:, 1]
            
            # Add UMAP coordinates with merge approach
            if 'umap' in data.tl.result:
                umap_coords = data.tl.result['umap'].copy()
                umap_coords['cell_name'] = data.cells.cell_name
                
                # Merge DataFrames for correct alignment
                annotation_df = pd.merge(
                    annotation_df,
                    umap_coords[['cell_name', 0, 1]], 
                    on='cell_name',
                    how='left'
                )
                
                # Rename UMAP columns
                annotation_df.rename(columns={0: 'umap1', 1: 'umap2'}, inplace=True)
            
            # Add custom gene scores if this is custom_genes annotation
            if template_name == 'custom_genes':
                cluster_score_map = {int(k): v for k, v in custom_scores.items()}
                annotation_df['gene_interest_score'] = annotation_df['cluster_id'].map(cluster_score_map)
                
                # Add best category for each cluster
                cluster_category_map = {}
                for cluster_str, details in annotation_details.items():
                    cluster_category_map[int(cluster_str)] = details['best_category'] or 'none'
                annotation_df['best_gene_category'] = annotation_df['cluster_id'].map(cluster_category_map)
            
            # Save CSV file
            annotation_csv = os.path.join(annotation_output_dir, 'exports', 'csv', f'{template_name}_annotations.csv')
            annotation_df.to_csv(annotation_csv, index=False)
            log_step(f"{template_name.capitalize()} annotation assignments saved")
            
    except Exception as e:
        log_step(f"Warning: Error applying {template_name} annotation: {e}")
        continue

# Strategy 5: Create individual cluster plots for each annotation type
log_step("="*80)
log_step("CREATING INDIVIDUAL CLUSTER PLOTS FOR EACH ANNOTATION")
log_step("="*80)

try:
    individual_plots_success = create_individual_cluster_plots_per_annotation(
        data, annotation_templates, annotation_output_dir
    )
    
    if individual_plots_success:
        log_step("SUCCESS: Individual cluster plots created for all annotation types")
        log_step("Directory structure created:")
        log_step(f"  - Main directory: plots/individual_clusters/")
        
        for template_name in annotation_templates.keys():
            if annotation_templates[template_name]:
                subdir_path = os.path.join(annotation_output_dir, 'plots', 'individual_clusters', template_name)
                if os.path.exists(subdir_path):
                    file_count = len([f for f in os.listdir(subdir_path) if f.endswith('.png')])
                    log_step(f"  - {template_name}/: {file_count} individual plots + HTML summary")
        
        log_step("Access via: plots/individual_clusters/index.html")
    else:
        log_step("WARNING: Individual cluster plots creation failed")

except Exception as e:
    log_step(f"ERROR in individual cluster plots creation: {e}")
    import traceback
    log_step(f"Traceback: {traceback.format_exc()}")

# Save individual plots checkpoint
save_progress_checkpoint(data, annotation_output_dir, 'individual_plots_complete')
log_memory_usage("INDIVIDUAL_PLOTS_COMPLETE")

# DATA EXPORT

log_step("="*80)
log_step("DATA EXPORT")
log_step("="*80)

export_start_time = datetime.now()

# Export to AnnData format
try:
    log_step("Converting to AnnData format")
    
    # Choose the best annotation for primary export
    primary_annotation = 'biological' if biological_annotations else 'alphabetical'
    
    # Create comprehensive output filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    h5ad_filename = f"D02266B1_stereopy_complete_analysis_{timestamp}.h5ad"
    h5ad_path = os.path.join(annotation_output_dir, 'exports', 'h5ad', h5ad_filename)
    
    log_step(f"Exporting to {h5ad_filename}...")
    
    adata = st.io.stereo_to_anndata(data, output=h5ad_path)
    
    log_step(f"AnnData export completed. File size: {os.path.getsize(h5ad_path) / 1e6:.1f} MB")
    
    # Add additional metadata to the exported AnnData object
    log_step("Adding comprehensive metadata to AnnData object")
    
    # Load the saved AnnData to add metadata
    adata = ad.read_h5ad(h5ad_path)
    
    # Add analysis metadata
    adata.uns['analysis_metadata'] = {
        'analysis_date': datetime.now().isoformat(),
        'stereopy_version': st.__version__,
        'data_source': data_path,
        'bin_size': 100,
        'preprocessing': {
            'min_counts': MIN_COUNTS,
            'min_genes': MIN_GENES,
            'pct_counts_mt': PCT_COUNTS_MT,
            'normalization': 'total_count_10k_log1p',
            'scaling': 'max10_no_zero_center',
            'hvg_params': {'min_mean': 0.0125, 'max_mean': 3, 'min_disp': 0.5, 'n_top_genes': 2000}
        },
        'clustering': {
            'method': 'louvain',
            'neighbors_params': {'n_pcs': N_PCS},
            'n_clusters': len(data.cells['louvain'].unique()) if 'louvain' in data.cells else 0
        },
        'annotations_available': list(annotation_templates.keys())
    }
    
    # Add cluster statistics
    if 'louvain' in adata.obs:
        cluster_stats = adata.obs['louvain'].value_counts().sort_index()
        adata.uns['cluster_statistics'] = {
            'cluster_sizes': cluster_stats.to_dict(),
            'total_clusters': len(cluster_stats),
            'largest_cluster': int(cluster_stats.max()),
            'smallest_cluster': int(cluster_stats.min()),
            'mean_cluster_size': float(cluster_stats.mean()),
            'median_cluster_size': float(cluster_stats.median())
        }
    
    # Save updated AnnData
    adata.write(h5ad_path)
    
    # Create summary report
    export_summary_file = os.path.join(annotation_output_dir, 'exports', 'export_summary.txt')
    with open(export_summary_file, 'w') as f:
        f.write("DATA EXPORT SUMMARY\n")
        f.write("="*60 + "\n\n")
        f.write(f"Export completed: {datetime.now()}\n")
        f.write(f"Primary export file: {h5ad_filename}\n")
        f.write(f"File size: {os.path.getsize(h5ad_path) / 1e6:.1f} MB\n\n")
        
        f.write("ANNDATA OBJECT STRUCTURE:\n")
        f.write(f"Observations (cells): {adata.n_obs}\n")
        f.write(f"Variables (genes): {adata.n_vars}\n")
        f.write(f"Observation keys: {', '.join(adata.obs.keys())}\n")
        f.write(f"Variable keys: {', '.join(adata.var.keys())}\n")
        f.write(f"Unstructured keys: {', '.join(adata.uns.keys())}\n")
        f.write(f"Obsm keys: {', '.join(adata.obsm.keys())}\n")
        f.write(f"Obsp keys: {', '.join(adata.obsp.keys())}\n\n")
        
        if 'cluster_statistics' in adata.uns:
            stats = adata.uns['cluster_statistics']
            f.write("CLUSTER STATISTICS:\n")
            f.write(f"Total clusters: {stats['total_clusters']}\n")
            f.write(f"Largest cluster: {stats['largest_cluster']} cells\n")
            f.write(f"Smallest cluster: {stats['smallest_cluster']} cells\n")
            f.write(f"Mean cluster size: {stats['mean_cluster_size']:.1f} cells\n")
            f.write(f"Median cluster size: {stats['median_cluster_size']:.1f} cells\n\n")
        
        f.write("AVAILABLE ANNOTATIONS:\n")
        for template_name in annotation_templates.keys():
            if f'anno_{template_name}' in adata.obs:
                unique_annotations = adata.obs[f'anno_{template_name}'].nunique()
                f.write(f"- {template_name}: {unique_annotations} unique labels\n")
    
    log_step("Export summary saved")
    
except Exception as e:
    log_step(f"ERROR in AnnData export: {e}")
    import traceback
    log_step(f"Traceback: {traceback.format_exc()}")

# Export cluster assignments and annotations to CSV
try:
    log_step("Exporting detailed CSV files...")
    
    # Comprehensive metadata export
    metadata_csv = os.path.join(annotation_output_dir, 'exports', 'metadata', 'comprehensive_cell_metadata.csv')
    
    if 'louvain' in data.cells:
        metadata_df = pd.DataFrame({
            'cell_name': data.cells.cell_name,
            'cluster_louvain': data.cells['louvain']
        })
        
        # Add all available annotations
        for template_name in annotation_templates.keys():
            res_key = f'anno_{template_name}'
            if res_key in data.cells:
                metadata_df[f'annotation_{template_name}'] = data.cells[res_key]
        
        # Add spatial coordinates
        if hasattr(data.cells, 'spatial'):
            spatial_coords = data.cells.spatial
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
        log_step(f"Comprehensive metadata saved: {len(metadata_df)} cells")

    # Export gene metadata
    gene_metadata_csv = os.path.join(annotation_output_dir, 'exports', 'metadata', 'gene_metadata.csv')
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

    # Export expression matrix (top variable genes for downstream analysis)
    if 'highly_variable_genes' in data.tl.result:
        # Get highly variable genes
        hvg_result = data.tl.result['highly_variable_genes']
        if isinstance(hvg_result, pd.DataFrame) and 'highly_variable' in hvg_result.columns:
            hvg_genes = hvg_result[hvg_result['highly_variable'] == True].index.tolist()
            print(f"Found {len(hvg_genes)} highly variable genes (HVGs).")
            print(f"Top 5: {hvg_genes[:5]}")
            
            if len(hvg_genes) == 0:
                print("No highly variable genes (HVGs) found")
            else:
                try:
                    gene_name_to_idx = {gene: idx for idx, gene in enumerate(data.genes.gene_name)}
                    gene_indices = []
                    valid_hvg_genes = []
                    
                    for gene in hvg_genes:
                        if gene in gene_name_to_idx:
                            gene_indices.append(gene_name_to_idx[gene])
                            valid_hvg_genes.append(gene)
                    
                    print(f"Found {len(gene_indices)} highly variable genes (HVGs) of {len(hvg_genes)}")
                    
                    hvg_sparse_matrix = data.exp_matrix[:, gene_indices]
                    hvg_matrix_dense = hvg_sparse_matrix.toarray()
                    
                    hvg_expression_matrix_df = pd.DataFrame(
                        hvg_matrix_dense,
                        index=data.cells.cell_name,
                        columns=valid_hvg_genes
                    )
                    
                    expression_file = os.path.join(annotation_output_dir, 'exports', 'csv', 'hvg_expression_matrix.csv')
                    
                    hvg_expression_matrix_df.to_csv(expression_file)
                    print(f"HVG matrix saved in {expression_file} with {hvg_expression_matrix_df.shape[0]} cells and {hvg_expression_matrix_df.shape[1]} genes")
                    
                except Exception as e1:
                    print(f"Error exporting HVG matrix: {e1}")
                    try:
                        print("Trying fallback")
                        
                        # Encontrar índices numericamente
                        gene_indices = []
                        for gene in hvg_genes[:100]:
                            try:
                                idx = data.genes.gene_name.tolist().index(gene)
                                gene_indices.append(idx)
                            except ValueError:
                                continue
                        
                        if len(gene_indices) > 0:
                            print(f"Exporting {len(gene_indices)} HVGs...")
                            submatrix = data.exp_matrix[:, gene_indices]
                            hvg_df = pd.DataFrame(
                                submatrix.toarray(),
                                index=data.cells.cell_name,
                                columns=[hvg_genes[hvg_genes.index(data.genes.gene_name.iloc[idx])] for idx in gene_indices if data.genes.gene_name.iloc[idx] in hvg_genes]
                            )
                            expression_file = os.path.join(annotation_output_dir, 'exports', 'csv', 'hvg_expression_matrix.csv')
                            hvg_df.to_csv(expression_file)
                            print(f"Fallback success: File saved with {hvg_df.shape[0]} cells and {hvg_df.shape[1]} genes")
                        else:
                            print("No valid HVG found for exporting")
                            
                    except Exception as e2:
                        print(f"Fallback method failed: {e2}")
                        print("Skipping HVG export")
        else:
            print("ERROR: No valid HVG in the dataframe.")

except Exception as e:
    log_step(f"Warning: Error in CSV export: {e}")
    print(f"Error in exporting: {str(e)}")

log_memory_usage("EXPORT_COMPLETE")

# Visualization

log_step("="*70)
log_step("CREATING DIRECT GENE INTEREST VISUALIZATION")
log_step("="*70)
direct_results = create_direct_gene_visualization()

if direct_results:
    log_step("SUCCESS: Direct gene visualization completed!")
    log_step("Check the plots/ directory for generated images")
else:
    log_step("FAILED: Could not create direct visualization")

# DATA INTEGRATION

if direct_results:
    log_step("SUCCESS: Direct gene visualization completed!")
    clean_spatial_success = create_clean_spatial_visualization_from_direct_results(
        direct_results, data, annotation_output_dir
    )
    
    if clean_spatial_success:
        log_step("SUCCESS: Clean spatial visualization integrated!")
    else:
        log_step("WARNING: Clean spatial visualization failed")
    
    log_step("Check the plots/ directory for all generated images")
else:
    log_step("FAILED: Could not create direct visualization")


# FINAL REPORT AND SUMMARY

log_step("="*70)
log_step("GENERATING FINAL ANNOTATION REPORT")
log_step("="*70)

end_time = datetime.now()
total_duration = (end_time - start_time).total_seconds()
export_duration = (end_time - export_start_time).total_seconds()

final_report_path = os.path.join(annotation_output_dir, 'ANNOTATION_ANALYSIS_REPORT.txt')
with open(final_report_path, 'w') as f:
    f.write("="*100 + "\n")
    f.write("STEREOPY ANNOTATION & EXPORT ANALYSIS REPORT\n")
    f.write("="*100 + "\n\n")
    f.write(f"Analysis completed: {end_time}\n")
    f.write(f"Total processing time: {total_duration/60:.1f} minutes ({total_duration:.1f} seconds)\n")
    f.write(f"Export processing time: {export_duration:.1f} seconds\n")
    f.write(f"Data source: {data_path}\n")
    f.write(f"Output directory: {annotation_output_dir}\n\n")
    
    # System utilization
    memory = psutil.virtual_memory()
    f.write("SYSTEM RESOURCES UTILIZED:\n")
    f.write(f"- Total RAM available: {memory.total/1e9:.1f} GB\n")
    f.write(f"- Peak memory usage: {memory.percent:.1f}% ({memory.used/1e9:.1f} GB)\n")
    f.write(f"- CPU cores utilized: {psutil.cpu_count()}\n\n")
    
    # Annotation results
    if 'louvain' in data.cells:
        cluster_counts = data.cells['louvain'].value_counts()
        f.write("ANNOTATION ANALYSIS RESULTS:\n")
        f.write(f"- Base clustering method: Louvain\n")
        f.write(f"- Total clusters: {len(cluster_counts)}\n")
        f.write(f"- Total cells: {cluster_counts.sum()}\n")
        f.write(f"- Annotation strategies applied: {len(annotation_templates)}\n")
        
        for template_name in annotation_templates.keys():
            if annotation_templates[template_name]:
                f.write(f"  * {template_name}: {len(annotation_templates[template_name])} cluster labels\n")
        f.write("\n")
    
    # Export summary
    f.write("DATA EXPORT SUMMARY:\n")
    f.write(f"- Primary format: AnnData (H5AD)\n")
    if os.path.exists(h5ad_path):
        f.write(f"- Main export file: {h5ad_filename} ({os.path.getsize(h5ad_path) / 1e6:.1f} MB)\n")
    
    # Directory structure
    f.write("OUTPUT DIRECTORY STRUCTURE:\n")
    f.write("- plots/: All visualization outputs\n")
    f.write("  * annotation/: Annotation-specific plots\n")
    f.write("  * comparison/: Comparison between methods\n")
    f.write("- annotations/: Annotation definitions and mappings\n")
    f.write("  * manual/: Manual annotation templates\n") 
    f.write("  * automated/: Algorithm-generated annotations\n")
    f.write("- exports/: Data export files\n")
    f.write("  * h5ad/: AnnData format exports\n")
    f.write("  * csv/: Tabular data exports\n")
    f.write("  * metadata/: Cell and gene metadata\n")
    f.write("- downstream_analysis/: Templates for further analysis\n")
    f.write("- biological_interpretation/: Marker-based annotations\n\n")

# Save processing summary
processing_summary_file = os.path.join(annotation_output_dir, 'logs', 'annotation_processing_summary.json')
processing_summary = {
    'analysis_type': 'ANNOTATION_AND_EXPORT',
    'start_time': start_time.isoformat(),
    'end_time': end_time.isoformat(),
    'total_duration_seconds': total_duration,
    'export_duration_seconds': export_duration,
    'annotation_strategies': len(annotation_templates),
    'successful_annotations': sum(1 for templates in annotation_templates.values() if templates),
    'export_file_size_mb': os.path.getsize(h5ad_path) / 1e6 if os.path.exists(h5ad_path) else 0,
    'peak_memory_percent': psutil.virtual_memory().percent,
    'peak_memory_gb': psutil.virtual_memory().used/1e9,
    'system_info': {
        'ram_gb': psutil.virtual_memory().total/1e9,
        'cpu_cores': psutil.cpu_count(),
        'python_version': sys.version.split()[0],
        'stereopy_version': st.__version__
    }
}

with open(processing_summary_file, 'w') as f:
    json.dump(processing_summary, f, indent=2)

log_step(f"Final annotation report saved to {final_report_path}")
log_step(f"Processing summary saved to {processing_summary_file}")

# Final cleanup
gc.collect()
log_memory_usage("FINAL")

log_step("="*80)
log_step("ANNOTATION & EXPORT ANALYSIS COMPLETED SUCCESSFULLY!")
log_step("="*80)
log_step(f"Total processing time: {total_duration/60:.1f} minutes")
log_step(f"Export file size: {os.path.getsize(h5ad_path) / 1e6:.1f} MB" if os.path.exists(h5ad_path) else "Export file not found")
log_step(f"Results location: {annotation_output_dir}")
log_step("")

EOF

echo ""
echo "Annotation script created successfully!"
echo ""

# Execute the annotation analysis
echo "Starting annotation and export analysis"
echo ""

$ST_PYTHON bin/stereopy_annotation_analysis.py

if [ $? -eq 0 ]; then
    echo ""
    echo "==========================================="
    echo "ANNOTATION ANALYSIS COMPLETED SUCCESSFULLY!"
    echo "==========================================="
    echo ""
    echo "Complete results saved to: results_annotation/"
    echo "==========================================="
else
    echo ""
    echo "==========================================="
    echo "ANNOTATION ANALYSIS FAILED!"
    echo "==========================================="
    echo "Check error logs for details"
    echo "==========================================="
    exit 1
fi

echo "==========================================="
echo ""
echo "Final system status:"
free -h
echo ""
echo "ANNOTATION ANALYSIS FINISHED:"
date
echo "==========================================="
