lib_path <- paste0(getwd(), "/R_libs")
if (!dir.exists(lib_path)) dir.create(lib_path, recursive = TRUE)
.libPaths(c(lib_path, .libPaths()))

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = lib_path, repos = "https://cloud.r-project.org")
  }
}

needed_pkgs <- c("Seurat", "WGCNA", "tidyverse", "Matrix", "hdWGCNA")
lapply(needed_pkgs, install_if_missing)

library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(Matrix)

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 10)

args <- commandArgs(trailingOnly = TRUE)
network_gene_set <- if (length(args) >= 1 && nzchar(args[1])) tolower(args[1]) else tolower(Sys.getenv("NETWORK_GENE_SET", "hvg"))
custom_gene_list_path <- if (length(args) >= 2 && nzchar(args[2])) args[2] else Sys.getenv("CUSTOM_GENE_LIST", "INPUT/interest_genes.txt")

if (!network_gene_set %in% c("hvg", "all", "custom")) {
  warning(paste0("Invalid NETWORK_GENE_SET value '", network_gene_set, "'. Falling back to 'hvg'."))
  network_gene_set <- "hvg"
}

resolve_gene_list_path <- function(path) {
  if (file.exists(path)) {
    return(path)
  }
  alt_path <- file.path("..", path)
  if (file.exists(alt_path)) {
    return(alt_path)
  }
  alt_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (file.exists(alt_path)) {
    return(alt_path)
  }
  path
}

load_custom_gene_list <- function(path) {
  resolved_path <- resolve_gene_list_path(path)
  if (!file.exists(resolved_path)) {
    warning(paste0("Custom gene list not found at '", path, "'. Falling back to HVGs."))
    return(character())
  }

  raw_genes <- readLines(resolved_path, warn = FALSE)
  raw_genes <- trimws(raw_genes)
  raw_genes <- raw_genes[nzchar(raw_genes)]
  if (length(raw_genes) == 0) {
    return(character())
  }
  if (any(grepl("[,\t]", raw_genes))) {
    raw_genes <- unlist(strsplit(raw_genes, "[,\t]"))
    raw_genes <- trimws(raw_genes)
    raw_genes <- raw_genes[nzchar(raw_genes)]
  }
  unique(raw_genes)
}

build_modules_df <- function(mergedColors, module_table = NULL, kME_matrix = NULL) {
  modules_df <- data.frame(
    gene_name = names(mergedColors),
    module = as.character(mergedColors),
    color = as.character(mergedColors),
    stringsAsFactors = FALSE
  )

  modules_df$kME <- NA_real_

  if (!is.null(module_table)) {
    module_table <- as.data.frame(module_table, stringsAsFactors = FALSE)
    if ("gene_name" %in% colnames(module_table)) {
      rownames(module_table) <- module_table$gene_name
      for (i in seq_len(nrow(modules_df))) {
        gene_name <- modules_df$gene_name[i]
        mod <- modules_df$module[i]
        kme_col <- paste0("kME_", mod)
        if (gene_name %in% rownames(module_table) && kme_col %in% colnames(module_table)) {
          modules_df$kME[i] <- module_table[gene_name, kme_col]
        }
      }
      return(modules_df)
    }
  }

  if (!is.null(kME_matrix)) {
    for (i in seq_len(nrow(modules_df))) {
      mod <- modules_df$module[i]
      me_col <- paste0("ME", mod)
      if (me_col %in% colnames(kME_matrix) && i <= nrow(kME_matrix)) {
        modules_df$kME[i] <- kME_matrix[i, me_col]
      }
    }
  }

  modules_df
}

compute_kme_matrix <- function(datExpr, mergedMEs) {
  datExpr_df <- as.data.frame(datExpr)
  mergedMEs_df <- as.data.frame(mergedMEs)

  if (nrow(datExpr_df) == nrow(mergedMEs_df)) {
    return(cor(datExpr_df, mergedMEs_df, use = "p"))
  }

  if (ncol(datExpr_df) == nrow(mergedMEs_df)) {
    return(cor(t(as.matrix(datExpr_df)), as.matrix(mergedMEs_df), use = "p"))
  }

  stop(
    paste0(
      "Unable to align datExpr and mergedMEs for kME calculation. datExpr dims: ",
      nrow(datExpr_df), "x", ncol(datExpr_df),
      "; mergedMEs dims: ",
      nrow(mergedMEs_df), "x", ncol(mergedMEs_df)
    )
  )
}

build_network_fallback <- function(seurat_obj, k_value = 25) {
  cat("\n[PROCESS] Using fallback k-means metacells + WGCNA pipeline\n")
  pca_embeddings <- Embeddings(seurat_obj, reduction = "pca")[, 1:30]
  set.seed(12345)
  km_result <- kmeans(pca_embeddings, centers = k_value, iter.max = 100, nstart = 25)

  expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  metacell_expr <- matrix(0, nrow = nrow(expr_matrix), ncol = k_value)
  rownames(metacell_expr) <- rownames(expr_matrix)
  colnames(metacell_expr) <- paste0("MC_", seq_len(k_value))

  for (i in seq_len(k_value)) {
    cluster_cells <- which(km_result$cluster == i)
    if (length(cluster_cells) > 0) {
      metacell_expr[, i] <- Matrix::rowMeans(expr_matrix[, cluster_cells, drop = FALSE])
    }
  }

  datExpr <- as.data.frame(t(metacell_expr))
  vars <- apply(datExpr, 2, var)
  bad_genes <- names(vars[vars == 0 | is.na(vars)])
  if (length(bad_genes) > 0) {
    write.table(bad_genes, paste0(output_dir, "INITIAL_DISCARD_GSG.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    datExpr <- datExpr[, !colnames(datExpr) %in% bad_genes, drop = FALSE]
  }

  powers <- c(seq(1, 10, by = 1), seq(12, 30, by = 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "unsigned", verbose = 5)
  selected_power <- sft$powerEstimate
  if (is.na(selected_power) || !is.finite(selected_power) || selected_power < 1 || selected_power > 50) {
    selected_power <- 6
  }

  adjacency <- adjacency(datExpr, power = selected_power, type = "unsigned")
  TOM <- TOMsimilarity(adjacency)
  dissTOM <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
  dynamicColors <- labels2colors(dynamicMods)
  merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)
  mergedColors <- merge$colors
  if (is.null(names(mergedColors))) {
    names(mergedColors) <- colnames(datExpr)
  }

  kME_matrix <- compute_kme_matrix(datExpr, merge$newMEs)
  modules_df <- build_modules_df(mergedColors, kME_matrix = kME_matrix)

  list(
    datExpr = datExpr,
    TOM = TOM,
    geneTree = geneTree,
    mergedColors = mergedColors,
    mergedMEs = merge$newMEs,
    modules_df = modules_df
  )
}

build_network_hdwgcna <- function(seurat_obj, gene_set, custom_gene_list_path, k_value = 25) {
  cat("\n[PROCESS] Using hdWGCNA metacell workflow\n")

  selected_features <- NULL
  if (gene_set == "custom") {
    selected_features <- load_custom_gene_list(custom_gene_list_path)
    if (length(selected_features) == 0) {
      warning("Custom gene list is empty or unavailable. Falling back to HVGs.")
      gene_set <- "hvg"
    }
  }

  if (gene_set == "custom") {
    seurat_obj <- SetupForWGCNA(
      seurat_obj,
      wgcna_name = "poly_pipeline",
      features = selected_features
    )
  } else if (gene_set == "all") {
    seurat_obj <- SetupForWGCNA(
      seurat_obj,
      wgcna_name = "poly_pipeline",
      gene_select = "all"
    )
  } else {
    seurat_obj <- SetupForWGCNA(
      seurat_obj,
      wgcna_name = "poly_pipeline",
      gene_select = "variable"
    )
  }

  seurat_obj$all_cells <- "all"
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = "all_cells",
    ident.group = "all_cells",
    k = k_value,
    target_metacells = 250,
    min_cells = 0,
    max_shared = 10,
    reduction = "pca",
    assay = "RNA",
    slot = "data",
    layer = "data"
  )
  seurat_obj <- NormalizeMetacells(seurat_obj)
  seurat_obj <- SetDatExpr(seurat_obj, group_name = "all", group.by = "all_cells", use_metacells = TRUE)
  seurat_obj <- TestSoftPowers(seurat_obj, networkType = "unsigned")

  soft_power <- NA_real_
  power_table <- tryCatch(as.data.frame(GetPowerTable(seurat_obj)), error = function(e) NULL)
  if (!is.null(power_table) && all(c("Power", "SFT.R.sq") %in% colnames(power_table))) {
    candidate_rows <- power_table %>%
      filter(is.finite(Power), is.finite(SFT.R.sq), Power >= 1, Power <= 50, SFT.R.sq >= 0.8)
    if (nrow(candidate_rows) > 0) {
      soft_power <- min(candidate_rows$Power)
    } else {
      finite_rows <- power_table %>%
        filter(is.finite(Power), is.finite(SFT.R.sq), Power >= 1, Power <= 50)
      if (nrow(finite_rows) > 0) {
        soft_power <- finite_rows$Power[which.max(finite_rows$SFT.R.sq)]
      }
    }
  }

  if (is.na(soft_power) || !is.finite(soft_power) || soft_power < 1 || soft_power > 50) {
    soft_power <- 6
    cat(paste0("[PROCESS] No valid hdWGCNA soft power found; falling back to soft_power = ", soft_power, "\n"))
  } else {
    cat(paste0("[PROCESS] Using hdWGCNA soft_power = ", soft_power, "\n"))
  }

  seurat_obj <- ConstructNetwork(seurat_obj, soft_power = soft_power, tom_name = "poly_pipeline", overwrite_tom = TRUE)

  if ("ModuleEigengenes" %in% getNamespaceExports("hdWGCNA")) {
    seurat_obj <- ModuleEigengenes(seurat_obj)
  }

  datExpr <- GetDatExpr(seurat_obj)
  TOM <- GetTOM(seurat_obj)
  modules <- as.data.frame(GetModules(seurat_obj), stringsAsFactors = FALSE)
  mergedColors <- modules$color
  names(mergedColors) <- modules$gene_name
  modules_df <- build_modules_df(mergedColors, module_table = modules)

  mergedMEs <- tryCatch(
    as.data.frame(GetMEs(seurat_obj, harmonized = FALSE)),
    error = function(e) NULL
  )
  if (is.null(mergedMEs)) {
    mergedMEs <- tryCatch(
      as.data.frame(GetMEs(seurat_obj, harmonized = TRUE)),
      error = function(e) NULL
    )
  }
  if (is.null(mergedMEs)) {
    stop("hdWGCNA completed but module eigengenes were not available.")
  }

  dissTOM <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")

  list(
    seurat_obj = seurat_obj,
    datExpr = datExpr,
    TOM = TOM,
    geneTree = geneTree,
    mergedColors = mergedColors,
    mergedMEs = mergedMEs,
    modules_df = modules_df
  )
}

use_hdWGCNA <- all(c("SetupForWGCNA", "MetacellsByGroups", "NormalizeMetacells", "SetDatExpr", "TestSoftPowers", "ConstructNetwork", "GetDatExpr", "GetTOM", "GetModules", "GetMEs") %in% getNamespaceExports("hdWGCNA"))

all_dirs <- list.dirs("..", full.names = TRUE, recursive = FALSE)
results_dirs <- all_dirs[grepl("/RESULTS_", all_dirs)]

if (length(results_dirs) == 0) {
  stop("ERROR: No RESULTS folder found (RESULTS_*)")
}

latest_results <- sort(results_dirs, decreasing = TRUE)[1]
cat(paste0("Working on: ", latest_results, "\n"))

setwd(latest_results)

output_dir <- "NETWORK/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
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
  modules_df <- checkpoint_data$modules_df
} else {
  cat("\n[PROCESS] Starting full processing\n")

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

  k_value <- 25
  network_state <- NULL

  if (use_hdWGCNA) {
    network_state <- tryCatch(
      build_network_hdwgcna(seurat_obj, network_gene_set, custom_gene_list_path, k_value = k_value),
      error = function(e) {
        warning(paste0("hdWGCNA workflow failed, falling back to legacy pipeline: ", conditionMessage(e)))
        NULL
      }
    )
  }

  if (is.null(network_state)) {
    network_state <- build_network_fallback(seurat_obj, k_value = k_value)
  }

  datExpr <- network_state$datExpr
  TOM <- network_state$TOM
  geneTree <- network_state$geneTree
  mergedColors <- network_state$mergedColors
  mergedMEs <- network_state$mergedMEs
  modules_df <- network_state$modules_df

  saveRDS(list(datExpr = datExpr, TOM = TOM, geneTree = geneTree, mergedColors = mergedColors, mergedMEs = mergedMEs, modules_df = modules_df), checkpoint_file)
}

if (is.null(modules_df)) {
  kME <- tryCatch(
    compute_kme_matrix(datExpr, mergedMEs),
    error = function(e) {
      warning(paste0("Unable to reconstruct kME values from checkpoint data: ", conditionMessage(e), ". Exporting nodes with kME = NA."))
      NULL
    }
  )
  modules_df <- build_modules_df(mergedColors, kME_matrix = kME)
}

threshold <- 0.15
all_edges_list <- list()

for (mod in unique(mergedColors)) {
  mod_genes <- modules_df %>% filter(module == mod) %>% pull(gene_name)
  mod_idx <- which(colnames(datExpr) %in% mod_genes)
  if (length(mod_genes) < 2) next
  tom_sub <- TOM[mod_idx, mod_idx]
  rownames(tom_sub) <- colnames(tom_sub) <- mod_genes

  edges_mod <- as.data.frame(as.table(tom_sub), stringsAsFactors = FALSE)
  colnames(edges_mod) <- c("fromNode", "toNode", "weight")

  edges_mod <- edges_mod %>%
    filter(weight > threshold & fromNode != toNode)

  if (nrow(edges_mod) > 0) {
    edges_mod <- edges_mod[as.character(edges_mod$fromNode) < as.character(edges_mod$toNode), ]
    all_edges_list[[mod]] <- edges_mod
    write.table(edges_mod, paste0(output_dir, mod, "_EDGE.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  write.table(modules_df %>% filter(module == mod), paste0(output_dir, mod, "_NODE.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

exports_dir <- "EXPORTS/"
if (!dir.exists(exports_dir)) dir.create(exports_dir, recursive = TRUE)
edge_filename <- paste0(exports_dir, project_name, "_FULL_EDGES.txt")
node_filename <- paste0(exports_dir, project_name, "_FULL_NODES.txt")
cat(paste0("Exporting Complete Edges and Nodes files for visualization: ", exports_dir, "\n"))

if (length(all_edges_list) > 0) {
  write.table(bind_rows(all_edges_list), edge_filename, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(data.frame(), edge_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

write.table(modules_df, node_filename, sep = "\t", quote = FALSE, row.names = FALSE)