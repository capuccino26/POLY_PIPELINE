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
