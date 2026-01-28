import stereo as st
import os
import glob
import warnings
import numpy as np
import time
import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import argparse
import sys

# Function for logging stamps
def log_step(message):
    timestamp = datetime.datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()

process_stt = time.time()

warnings.filterwarnings('ignore')

# Parse command line arguments
parser = argparse.ArgumentParser(description='Plot gene expression from GEF files')
parser.add_argument('-o', '--output', type=str, default='', 
                    help='Output subdirectory name (e.g., CLUSTER10). If not specified, uses default directory.')
parser.add_argument('-t', '--thr', type=float, default=1.0, help='Minimum Threshold (e.g., 10.5) for Differential Expression')
parser.add_argument('--min_x', type=float, default=None, help='Minimum X coordinate for filtering')
parser.add_argument('--max_x', type=float, default=None, help='Maximum X coordinate for filtering')
parser.add_argument('--min_y', type=float, default=None, help='Minimum Y coordinate for filtering')
parser.add_argument('--max_y', type=float, default=None, help='Maximum Y coordinate for filtering')
parser.add_argument('-i', '--interest', type=str, required=True,
                    help='Path to the file containing the list of genes of interest (one gene per line)')
args = parser.parse_args()

# Directories

input_dir = "INPUT/datasets"
if args.output:
    main_output_dir = os.path.join("RESULTS/output_images", args.output)
else:
    main_output_dir = "RESULTS/output_images/"

# Read genes from file
gene_file = args.interest
try:
    with open(gene_file, 'r') as f:
        int_genes = [line.strip() for line in f if line.strip()]
    log_step(f"Loaded {len(int_genes)} genes from {gene_file}")
    log_step(f"Genes: {int_genes}")
    if not int_genes:
        log_step(f"Warning: {gene_file} is empty. Check example for correct usage.")
        sys.exit(1)

    log_step(f"Loaded {len(int_genes)} genes from {gene_file}")
    log_step(f"Genes: {int_genes}")
except FileNotFoundError:
    log_step(f"ERROR: The interest gene file '{gene_file}' was not found.")
    log_step("The script requires a valid gene file to run. Check example for correct usage.")
    sys.exit(1)

except Exception as e:
    log_step(f"Error reading {gene_file}: {e}")
    log_step("The script requires a valid gene file to run. Check example for correct usage.")
    sys.exit(1)

# EXPRESSION THRESHOLD
EXPRESSION_THRESHOLD = args.thr

# Create output directories
os.makedirs(main_output_dir, exist_ok=True)
total_counts_dir = os.path.join(main_output_dir, "maps_counts")
os.makedirs(total_counts_dir, exist_ok=True)
for gene in int_genes:
    gene_output_dir = os.path.join(main_output_dir, gene)
    os.makedirs(gene_output_dir, exist_ok=True)

gef_files = glob.glob(os.path.join(input_dir, '*.gef'))
log_step(f"Found {len(gef_files)} .gef files to process.")
log_step(f"Using expression threshold: > {EXPRESSION_THRESHOLD}")

# Initialize statistics storage
all_sample_stats = []

for gef_path in gef_files:
    sample_name = os.path.splitext(os.path.basename(gef_path))[0]
    log_step(f"\nProcessing {sample_name}")

    try:
        # Load the data
        data = st.io.read_gef(file_path=gef_path)
        
        # Apply coordinate filtering if specified
        if any([args.min_x, args.max_x, args.min_y, args.max_y]):
            data.tl.filter_coordinates(
                min_x=args.min_x if args.min_x is not None else data.position[:, 0].min(),
                max_x=args.max_x if args.max_x is not None else data.position[:, 0].max(),
                min_y=args.min_y if args.min_y is not None else data.position[:, 1].min(),
                max_y=args.max_y if args.max_y is not None else data.position[:, 1].max()
            )
        
        # Calculate total counts
        total_counts_array = data.exp_matrix.sum(axis=1).A1
        
        # Identify spots with total_counts > 0
        tissue_spots_mask = total_counts_array > 0
        
        # Coordinates for all tissue spots (background)
        x_tissue = data.position[tissue_spots_mask, 0]
        y_tissue = data.position[tissue_spots_mask, 1]

        # Sample statistics
        sample_stats = {
            'sample_name': sample_name,
            'total_tissue_spots': int(np.sum(tissue_spots_mask)),
            'genes': {}
        }

        # PLOT Total Counts Map
        log_step("Generating total counts map")
        fig, ax = plt.subplots(figsize=(10, 10), facecolor='black')
        ax.set_facecolor('black')
        
        # Plot only tissue spots, colored by total counts
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
        else:
            ax.set_title("No tissue spots found with expression.", color='white')

        ax.set_title(f"Total Counts - {sample_name}", fontsize=16, color='white')
        ax.set_xlabel('Spt X', color='white')
        ax.set_ylabel('Spt Y', color='white')
        ax.tick_params(colors='white')
        ax.spines['bottom'].set_color('white')
        ax.spines['top'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.spines['right'].set_color('white')
        ax.set_aspect('equal', adjustable='box')
        
        path_counts = os.path.join(total_counts_dir, f"{sample_name}_thr_{EXPRESSION_THRESHOLD}.png")
        plt.savefig(path_counts, dpi=300, bbox_inches='tight', facecolor='black')
        plt.close(fig)
        log_step("Total counts map saved.")

        # Plot Individual Genes (dark tissue background) ---
        for gene in int_genes:
            log_step(f"\nGenerating map for {gene}")
            path_gene = os.path.join(main_output_dir, gene, f"{sample_name}_{gene}_thr_{EXPRESSION_THRESHOLD}.png")
            
            gene_indices = np.where(data.genes.gene_name == gene)[0]
            
            if len(gene_indices) > 0:
                gene_index = gene_indices[0]
                gene_expression_vector = data.exp_matrix[:, gene_index].toarray().flatten()

                fig, ax = plt.subplots(figsize=(12, 10), facecolor='black')
                ax.set_facecolor('black')

                # Layer 1: Background in grayscale
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
                
                # Layer 2: Expressed spots
                expressed_gene_spots_mask = (gene_expression_vector > EXPRESSION_THRESHOLD) & tissue_spots_mask
                
                if np.sum(expressed_gene_spots_mask) > 0:
                    x_expressed_gene = data.position[expressed_gene_spots_mask, 0]
                    y_expressed_gene = data.position[expressed_gene_spots_mask, 1]
                    expression_values_gene = gene_expression_vector[expressed_gene_spots_mask]
                    norm_gene = Normalize(vmin=0, vmax=np.max(expression_values_gene))
                    scatter_gene = ax.scatter(x_expressed_gene, y_expressed_gene,
                                              c=expression_values_gene, 
                                              cmap='plasma',
                                              s=35,
                                              norm=norm_gene, 
                                              edgecolors='white', 
                                              linewidths=0.3, 
                                              label=f'{gene}',
                                              alpha=0.95)

                    # Colorbar
                    cbar = plt.colorbar(scatter_gene, ax=ax, label=f'{gene} Expression (>{EXPRESSION_THRESHOLD})')
                    cbar.ax.yaxis.set_tick_params(color='white')
                    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
                    cbar.set_label(f'{gene} Expression (>{EXPRESSION_THRESHOLD})', color='white')
                    expressed_spots = int(np.sum(expressed_gene_spots_mask))
                    max_expr = float(np.max(expression_values_gene))
                    min_expr = float(np.min(expression_values_gene))
                    mean_expr = float(np.mean(expression_values_gene))
                    
                    log_step(f"{gene}: {expressed_spots} spots expressed (>{EXPRESSION_THRESHOLD}), max expression: {max_expr:.2f}")
                    sample_stats['genes'][gene] = {
                        'expressed_spots': expressed_spots,
                        'max_expression': max_expr,
                        'min_expression': min_expr,
                        'mean_expression': mean_expr,
                        'percentage_of_tissue': (expressed_spots / sample_stats['total_tissue_spots'] * 100) if sample_stats['total_tissue_spots'] > 0 else 0
                    }
                    individual_stats_file = os.path.join(main_output_dir, gene, f"{sample_name}_{gene}_stats_thr{EXPRESSION_THRESHOLD}.txt")
                    with open(individual_stats_file, 'w') as f_ind:
                        f_ind.write(f"{gene}\n")
                        f_ind.write(f"Sample: {sample_name}\n")
                        f_ind.write(f"Analysis date: {datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')}\n")
                        f_ind.write(f"Expression threshold used: > {EXPRESSION_THRESHOLD}\n\n")
                        f_ind.write(f"Total tissue spots: {sample_stats['total_tissue_spots']}\n\n")
                        f_ind.write(f"Expressed spots: {expressed_spots}\n")
                        f_ind.write(f"Percentage of tissue expressing gene: {(expressed_spots / sample_stats['total_tissue_spots'] * 100):.2f}%\n")
                        f_ind.write(f"Maximum expression level: {max_expr:.2f}\n")
                        f_ind.write(f"Minimum expression level: {min_expr:.2f}\n")
                        f_ind.write(f"Mean expression level: {mean_expr:.2f}\n")
                    
                else:
                    log_step(f"{gene} found but not expressed in any tissue spot (threshold > {EXPRESSION_THRESHOLD}).")
                    sample_stats['genes'][gene] = {
                        'expressed_spots': 0,
                        'max_expression': 0,
                        'min_expression': 0,
                        'mean_expression': 0,
                        'percentage_of_tissue': 0
                    }
                    individual_stats_file = os.path.join(main_output_dir, gene, f"{sample_name}_{gene}_stats_thr{EXPRESSION_THRESHOLD}.txt")
                    with open(individual_stats_file, 'w') as f_ind:
                        f_ind.write(f"{gene} but not expressed above threshold\n")
                        f_ind.write(f"Sample: {sample_name}\n")
                        f_ind.write(f"Analysis date: {datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')}\n")
                        f_ind.write(f"Expression threshold used: > {EXPRESSION_THRESHOLD}\n\n")
                        f_ind.write(f"Total tissue spots: {sample_stats['total_tissue_spots']}\n\n")
                        f_ind.write(f"Expressed spots: 0\n")
                        f_ind.write(f"Percentage of tissue expressing gene: 0.00%\n")
                        f_ind.write(f"Expression levels: Below threshold ({EXPRESSION_THRESHOLD})\n")
                    ax.text(0.5, 0.95, f'{gene} not detected\n(expression ≤ {EXPRESSION_THRESHOLD})', 
                           transform=ax.transAxes, fontsize=12, 
                           bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7),
                           ha='center', va='top', color='black')
                ax.set_title(f"Spatial Expression of {gene} - {sample_name}\n(Expression threshold > {EXPRESSION_THRESHOLD})", 
                            fontsize=14, color='white', pad=20)
                ax.set_xlabel('Spt X', fontsize=12, color='white')
                ax.set_ylabel('Spt Y', fontsize=12, color='white')
                ax.tick_params(colors='white')
                ax.spines['bottom'].set_color('white')
                ax.spines['top'].set_color('white')
                ax.spines['left'].set_color('white')
                ax.spines['right'].set_color('white')
                ax.set_aspect('equal', adjustable='box')
                if np.sum(expressed_gene_spots_mask) > 0:
                    legend = ax.legend(loc='upper left', bbox_to_anchor=(0, 1), 
                                      framealpha=0.8, facecolor='black', edgecolor='white')
                    for text in legend.get_texts():
                        text.set_color('white')
                ax.grid(False)
                plt.tight_layout()
                plt.savefig(path_gene, dpi=300, bbox_inches='tight', facecolor='black')
                plt.close(fig)
                log_step(f"Plot saved for gene {gene}.")
            else:
                log_step(f"{gene} not found in sample {sample_name}.")
                sample_stats['genes'][gene] = {
                    'expressed_spots': 'Gene not found',
                    'max_expression': 'Gene not found',
                    'min_expression': 'Gene not found',
                    'mean_expression': 'Gene not found',
                    'percentage_of_tissue': 'Gene not found'
                }
                individual_stats_file = os.path.join(main_output_dir, gene, f"{sample_name}_{gene}_stats_thr{EXPRESSION_THRESHOLD}.txt")
                with open(individual_stats_file, 'w') as f_ind:
                    f_ind.write(f"{gene}\n")
                    f_ind.write(f"Sample: {sample_name}\n")
                    f_ind.write(f"Analysis date: {datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')}\n")
                    f_ind.write(f"Expression threshold used: > {EXPRESSION_THRESHOLD}\n\n")
                    f_ind.write(f"Total tissue spots: {sample_stats['total_tissue_spots']}\n\n")
                    f_ind.write(f"Status: GENE NOT FOUND in this sample\n")

        # Add sample statistics to overall list
        all_sample_stats.append(sample_stats)

    except Exception as e:
        log_step(f"Error processing sample {sample_name}: {e}")
        
process_fit = time.time()
process_tot = process_fit - process_stt
timestamp_final = datetime.datetime.now().strftime("[%Y/%m/%d %H:%M:%S]")

# Save stats
stats_file = os.path.join(main_output_dir, f"expression_stats_thr{EXPRESSION_THRESHOLD}.txt")
with open(stats_file, 'w') as f:
    f.write("SPATIAL EXPRESSION ANALYSIS STATS\n")
    f.write(f"Analysis completed: {timestamp_final}\n")
    f.write(f"Processing: {process_tot:.2f} seconds\n")
    f.write(f"Expression threshold: > {EXPRESSION_THRESHOLD}\n")
    f.write(f"Genes analyzed: {', '.join(int_genes)}\n")
    f.write(f"Total samples processed: {len(all_sample_stats)}\n\n")
    
    for sample_stat in all_sample_stats:
        f.write(f"SAMPLE: {sample_stat['sample_name']}\n")
        f.write(f"Total tissue spots: {sample_stat['total_tissue_spots']}\n\n")
        
        for gene, stats in sample_stat['genes'].items():
            f.write(f"{gene}:\n")
            if isinstance(stats['expressed_spots'], int):
                f.write(f"  Expressed spots: {stats['expressed_spots']}\n")
                f.write(f"  Percentage of tissue: {stats['percentage_of_tissue']:.2f}%\n")
                f.write(f"  Max expression: {stats['max_expression']:.2f}\n")
                f.write(f"  Min expression: {stats['min_expression']:.2f}\n")
                f.write(f"  Mean expression: {stats['mean_expression']:.2f}\n")
            else:
                f.write(f"  Status: {stats['expressed_spots']}\n")
            f.write("\n")
        f.write("-" * 50 + "\n\n")

log_step(f"\n{timestamp_final}Finished processing all files! Duration: {process_tot:.2f}s")
