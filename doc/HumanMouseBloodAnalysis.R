## ----setup--------------------------------------------------------------------
# load SBF package
library(SBF)

## -----------------------------------------------------------------------------
# install packages
pkgs <- c("data.table", "dplyr", "matrixStats")
require_install <- pkgs[!(pkgs %in% row.names(installed.packages()))]
if (length(require_install))
  install.packages(require_install)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(matrixStats)
})

## -----------------------------------------------------------------------------
species <- c("Homo_sapiens", "Mus_musculus")
species_short <- sapply(species, getSpeciesShortName)
species_short
common_celltypes <- c("hsc", "clp", "cmp", "nkcell", "cd8tcell", "cd4tcell", 
                      "bcell", "cd14monocytes", "eosinophil", "neutrophil")

## -----------------------------------------------------------------------------
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
counts_list <- metadata_list <- list()
for (sp in species) {
  # read blood logTPM counts for each species
  counts <- read.table(paste0(path, "human_mouse_blood_counts/", sp,
                              "_blood_logtpm.tsv"), header = T, sep = "\t",
                       row.names = 1)
  info <- tstrsplit(colnames(counts), "_")
  metadata <- data.frame(project = info[[1]],
        species = info[[2]],
        tissue = info[[3]],
        gsm = info[[4]],
        name = colnames(counts),
        stringsAsFactors = F)
  metadata$ref <- seq_len(nrow(metadata))
  metadata$key <- paste0(metadata$species, "_", metadata$tissue)
  metadata$tissue_factor <- factor(metadata$tissue)
  counts_list[[sp]] <- counts
  metadata_list[[sp]] <- metadata
}
sapply(counts_list, dim)

## -----------------------------------------------------------------------------
avg_counts <- list()
for (sp in species) {
  avg_counts[[sp]] <- calcAvgCounts(counts_list[[sp]], metadata_list[[sp]])
}

## -----------------------------------------------------------------------------
# check tissue columns are matching in each species
c_tissues <- as.data.frame(sapply(avg_counts, function(x) tstrsplit(colnames(x),"_")[[2]]))
if (!all(apply(c_tissues, 1, function(x) all(x == x[1])))) {
        stop("Error! columns not matching")
}

## -----------------------------------------------------------------------------
sapply(avg_counts, dim)

## -----------------------------------------------------------------------------
# remove empty rows
removeZeros <- function(df) {
    return(df[rowSums(df) > 0, ])
}
avg_counts <- lapply(avg_counts, removeZeros)
sapply(avg_counts, dim)
# update counts_list
counts_list_sub <- list()
for (sp in names(avg_counts)) {
  counts_list_sub[[sp]] <- counts_list[[sp]][row.names(avg_counts[[sp]]), ,
                                                 drop = F]
}

## -----------------------------------------------------------------------------
# first lets compute OSBF without updating the initial estimate of V. U and Delta
# are updated in this case
cat(format(Sys.time(), "%a %b %d %X %Y"),"\n")
osbf_noVupdate <- SBF(avg_counts, orthogonal = TRUE, transform_matrix = TRUE, minimizeError = TRUE,
                      optimizeV = FALSE, tol = 1e-4)
cat(format(Sys.time(), "%a %b %d %X %Y"),"\n")
# Now lets compute OSBF updating all three factors (U, Delta, V)
cat("=======================================================\n")
cat("optimizing V = TRUE\n")
cat(format(Sys.time(), "%a %b %d %X %Y"),"\n")
osbf <- SBF(avg_counts, orthogonal = TRUE, transform_matrix = TRUE, minimizeError = TRUE, 
            optimizeV = TRUE, tol = 1e-4)
cat(format(Sys.time(), "%a %b %d %X %Y"),"\n")

## -----------------------------------------------------------------------------
cat("\n", sprintf("%-27s:", "Final error [No V update]"), sprintf("%16.2f",
                                                        osbf_noVupdate$error))
cat("\n", sprintf("%-27s:", "Final error [With V update]"), sprintf("%16.2f",
                                                                    osbf$error))
cat("\n", sprintf("%-27s:", "# of update [No V update]"), sprintf("%16d",
                                                      osbf_noVupdate$error_pos))
cat("\n", sprintf("%-27s:", "# of update [With V update]"), sprintf("%16d",
                                                                osbf$error_pos))

## -----------------------------------------------------------------------------
osbf_noVupdate$error / osbf$error

## -----------------------------------------------------------------------------
par(mfrow = c(1, 3))
plot(x = seq_len(length(osbf$error_vec)), y = osbf$error_vec,
     xlab = "# of updates (step 1 -> )",
     ylab = "Factorization error", col = "red")
plot(x = 10:length(osbf$error_vec),
     y = osbf$error_vec[10:length(osbf$error_vec)],
     xlab = "# of updates (step 10 -> )",
     ylab = "Factorization error", col = "red")
plot(x = 50:length(osbf$error_vec),
     y = osbf$error_vec[50:length(osbf$error_vec)],
     xlab = "# of updates (step 50 -> )",
     ylab = "Factorization error", col = "red")

## -----------------------------------------------------------------------------
zapsmall(osbf_noVupdate$v %*% t(osbf_noVupdate$v))

zapsmall(osbf$v %*% t(osbf$v))

## -----------------------------------------------------------------------------
percentInfo_noVupdate <- lapply(osbf_noVupdate$d, function(x) {
  round(x^2 / sum(x^2) * 100, 2)})
cat("\nPercentage for each delta [No V update]:")
for (i in names(osbf_noVupdate$d)) {
  cat("\n", sprintf("%-25s:", i), sprintf("%8.2f", percentInfo_noVupdate[[i]]))
}


percentInfo <- lapply(osbf$d, function(x) {
  round(x^2 / sum(x^2) * 100, 2) })
cat("\nPercentage for each delta:")
for (i in names(osbf$d)) {
  cat("\n", sprintf("%-25s:", i), sprintf("%8.2f", percentInfo[[i]]))
}

## -----------------------------------------------------------------------------
d_inv_NoVupdate <- list()
for (sp in names(avg_counts)) {
  if (length(osbf_noVupdate$d[[sp]]) == 1) {
    d_inv_NoVupdate[[sp]] <- as.matrix(diag(as.matrix(1 / osbf_noVupdate$d[[sp]])))
    } else {
      d_inv_NoVupdate[[sp]] <- as.matrix(diag(1 / osbf_noVupdate$d[[sp]]))
    }
}
d_inv <- list()
for (sp in names(avg_counts)) {
  if (length(osbf$d[[sp]]) == 1) {
    d_inv[[sp]] <- as.matrix(diag(as.matrix(1 / osbf$d[[sp]])))
    } else {
      d_inv[[sp]] <- as.matrix(diag(1 / osbf$d[[sp]]))
    }
}

## -----------------------------------------------------------------------------
# project using no V update estimates
# we can project both mean expression profiles as well as individual expression
# profiles
avgcounts_projected_noVupdate  <- counts_projected_noVupdate <- list()
n_noVupdate <- n1_noVupdate <- c()
for (sp in names(avg_counts)) {
  avgcounts_projected_noVupdate[[sp]] <- as.matrix(t(avg_counts[[sp]])) %*%
    as.matrix(osbf_noVupdate$u[[sp]]) %*% d_inv_NoVupdate[[sp]]
  n_noVupdate <- c(n_noVupdate, row.names(avgcounts_projected_noVupdate[[sp]]))
  counts_projected_noVupdate[[sp]] <- as.matrix(t(counts_list_sub[[sp]])) %*%
    as.matrix(osbf_noVupdate$u[[sp]]) %*% d_inv_NoVupdate[[sp]]
  n1_noVupdate <- c(n1_noVupdate, row.names(counts_projected_noVupdate[[sp]]))
}

# combine projected profiles
df_proj_avg_noVupdate <- as.data.frame(do.call(rbind, avgcounts_projected_noVupdate))
rownames(df_proj_avg_noVupdate) <- n_noVupdate
colnames(df_proj_avg_noVupdate) <- paste0("Dim", 1:ncol(df_proj_avg_noVupdate))
meta <- tstrsplit(row.names(df_proj_avg_noVupdate), "_")
df_proj_avg_noVupdate$tissue <- factor(meta[[2]])
df_proj_avg_noVupdate$species <- factor(meta[[1]])
df_proj_avg_noVupdate <- df_proj_avg_noVupdate %>% mutate(species = factor(species,
                            levels = species_short))


df_proj_noVupdate <- as.data.frame(do.call(rbind, counts_projected_noVupdate))
rownames(df_proj_noVupdate) <- n1_noVupdate
colnames(df_proj_noVupdate) <- paste0("Dim", 1:ncol(df_proj_noVupdate))
meta1 <- data.table::tstrsplit(row.names(df_proj_noVupdate), "_")
df_proj_noVupdate$tissue <- factor(meta1[[3]])
df_proj_noVupdate$species <- factor(meta1[[2]])
df_proj_noVupdate <- df_proj_noVupdate %>% mutate(species = factor(species,
                            levels = species_short))

# project using V update estimates
avgcounts_projected  <- counts_projected <- list()
n <- n1 <- c()
for (sp in names(avg_counts)) {
  avgcounts_projected[[sp]] <- as.matrix(t(avg_counts[[sp]])) %*%
    as.matrix(osbf$u[[sp]]) %*% d_inv[[sp]]
  n <- c(n, row.names(avgcounts_projected[[sp]]))
  counts_projected[[sp]] <- as.matrix(t(counts_list_sub[[sp]])) %*%
    as.matrix(osbf$u[[sp]]) %*% d_inv[[sp]]
  n1 <- c(n1, row.names(counts_projected[[sp]]))
}

# combine projected profiles
df_proj_avg <- as.data.frame(do.call(rbind, avgcounts_projected))
rownames(df_proj_avg) <- n
colnames(df_proj_avg) <- paste0("Dim", 1:ncol(df_proj_avg))
meta <- data.table::tstrsplit(row.names(df_proj_avg), "_")
df_proj_avg$tissue <- factor(meta[[2]])
df_proj_avg$species <- factor(meta[[1]])
df_proj_avg <- df_proj_avg %>% mutate(species = factor(species,
                            levels = species_short))


df_proj <- as.data.frame(do.call(rbind, counts_projected))
rownames(df_proj) <- n1
colnames(df_proj) <- paste0("Dim", 1:ncol(df_proj))
meta1 <- tstrsplit(row.names(df_proj), "_")
df_proj$tissue <- factor(meta1[[3]])
df_proj$species <- factor(meta1[[2]])
df_proj <- df_proj %>% mutate(species = factor(species,
                            levels = species_short))


## -----------------------------------------------------------------------------
# install packages
pkgs <- c("grid", "ggthemes", "ggplot2")
require_install <- pkgs[!(pkgs %in% row.names(installed.packages()))]
if (length(require_install))
  install.packages(require_install)
suppressPackageStartupMessages({
  library(grid)
  library(ggthemes)
  library(ggplot2)
})

## -----------------------------------------------------------------------------
# custom theme function for ggplot2
customTheme <- function(base_size = 10, base_family = "helvetica") {
   require(grid)
   require(ggthemes)
   (ggthemes::theme_foundation(base_size = base_size)
   + ggplot2::theme(plot.title = element_text(face = "bold",
                                   size = rel(1.2), hjust = 0.5),
         text = element_text(),
         panel.background = element_rect(colour = NA),
         plot.background = element_rect(colour = NA),
         panel.border = element_rect(colour = NA),
         axis.title = element_text(size = rel(1)),
         axis.title.y = element_text(angle = 90, vjust = 2),
         axis.title.x = element_text(vjust = -0.2),
         axis.text = element_text(),
         axis.line = element_line(colour = "black"),
         axis.ticks = element_line(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.key = element_rect(colour = NA),
         legend.position = "top",
         legend.direction = "horizontal",
         legend.key.size = unit(0.2, "cm"),
         legend.spacing = unit(0, "cm"),
         legend.title = element_text(face = "italic"),
         plot.margin = unit(c(10, 5, 5, 5), "mm"),
         strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
         strip.text = element_text(face = "bold")))
}

## -----------------------------------------------------------------------------
sel_colors <- c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "mediumturquoise" ,
                "#E6AB02","darkmagenta", "#666666", "black","darkolivegreen2")
i <- 1
j <- 2
ggplot2::ggplot(df_proj_noVupdate, aes(x = df_proj_noVupdate[, i],
                              y = df_proj_noVupdate[, j], col = tissue,
                              shape = species, fill = tissue)) +
  xlab(paste("OSBF Dim", i)) + ylab(paste("OSBF Dim", j)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) +
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme() +
  theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# 2D plot for Dim1 and Dim2 [With V update]
sel_colors <- c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "mediumturquoise",
                "#E6AB02","darkmagenta", "#666666", "black","darkolivegreen2")
i <- 1
j <- 2
ggplot2::ggplot(df_proj, aes(x = df_proj[, i],
                              y = df_proj[, j], col = tissue,
                              shape = species, fill = tissue)) +
  xlab(paste("OSBF Dim", i)) + ylab(paste("OSBF Dim", j)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) + 
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme() +
  theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# 2D plot for Dim2 and Dim3 [No V update]
sel_colors <- c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "mediumturquoise",
                "#E6AB02","darkmagenta", "#666666", "black","darkolivegreen2")
i <- 2
j <- 3
ggplot2::ggplot(df_proj_noVupdate, aes(x = df_proj_noVupdate[, i],
                              y = df_proj_noVupdate[, j], col = tissue,
                              shape = species, fill = tissue)) +
  xlab(paste("OSBF Dim", i)) + ylab(paste("OSBF Dim", j)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) + 
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme() +
  theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# 2D plot for Dim2 and Dim3 [With V update]
sel_colors <- c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "mediumturquoise",
                "#E6AB02","darkmagenta", "#666666", "black","darkolivegreen2")
i <- 2
j <- 3
ggplot2::ggplot(df_proj, aes(x = df_proj[, i], y = df_proj[, j], col = tissue,
                              shape = species, fill = tissue)) +  
  xlab(paste("OSBF Dim", i)) + ylab(paste("OSBF Dim", j)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) + 
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme() +
  theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# create output directories. Change this path accordingly
# outdir <- "~/Dropbox/0.Analysis/9.bloodAnalysis/"
# dir.create(file.path(outdir), showWarnings = FALSE)
# subDir <- "2Dplots"
# dir.create(file.path(outdir, subDir), showWarnings = FALSE)

sel_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "mediumturquoise",
                "#E6AB02", "darkmagenta", "#666666", "black", "darkolivegreen2")
outputname <- "human_mouse_blood"
finished <- c()
for (i in 1:(ncol(df_proj)-2)) {
  for (j in 1:(ncol(df_proj)-2)) {
    if (i == j) next
    if (j %in% finished) next
    ggplot(df_proj, aes(x = df_proj[, i], y = df_proj[, j],
                            col = tissue, shape = species, fill = tissue)) +  
      xlab(paste0("OSBF Dim ", i)) + 
      ylab(paste0("OSBF Dim ", j)) +
      geom_point(size = 1.5) +
      scale_color_manual(values = sel_colors) + 
      scale_shape_manual(values = c(21:25,3:7)) +
      scale_fill_manual(values = sel_colors) +
      customTheme(base_size = 12) + 
      theme(legend.title = element_blank())
    #ggsave(filename = paste0(outdir, "2Dplots/opt_2Dplot_Dim_", i, "-", j, "_",
    #                         outputname, ".pdf"), device = "pdf",
    #       width = 7, height = 7, useDingbats = FALSE)
  }
  finished <- c(finished,i)
}

## -----------------------------------------------------------------------------
# install packages
pkgs <- c("RColorBrewer")
require_install <- pkgs[!(pkgs %in% row.names(installed.packages()))]
if (length(require_install))
  install.packages(require_install)
pkgs <- c("ComplexHeatmap")
require_install <- pkgs[!(pkgs %in% row.names(installed.packages()))]
if (length(require_install)) {
  if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(RColorBrewer)
})

## -----------------------------------------------------------------------------
data <- df_proj
data$tissue <- NULL
data$species <- NULL
data <- as.matrix(data)
data_dist <- as.matrix(dist(data, method = "euclidean"))
meta <- data.table::tstrsplit(colnames(data_dist), "_")
ht <- ComplexHeatmap::HeatmapAnnotation(celltype = meta[[3]], species = meta[[2]],
                      col = list(species = c("hsapiens" = "#66C2A5",
                                             "mmusculus" = "#E78AC3"),
                                 celltype = c("bcell" = "#1B9E77",
                                            "cd14monocytes" = "#D95F02",
                                            "cd4tcell" ="#7570B3",
                                            "cd8tcell" = "#E7298A",
                                            "clp" = "mediumturquoise",
                                            "cmp" = "#E6AB02",
                                            "eosinophil" = "darkmagenta",
                                            "hsc" = "#666666",
                                            "neutrophil" = "black",
                                            "nkcell" = "darkolivegreen2")),
                                        annotation_name_side = "left")
mypalette <- RColorBrewer::brewer.pal(9, "Blues")
morecolors <- colorRampPalette(mypalette)

myheatmap <- ComplexHeatmap::Heatmap(as.matrix(data_dist), cluster_rows = T,
                     clustering_method_rows = "centroid",
                     cluster_columns = T,
                     clustering_method_columns = "centroid",
                     top_annotation = ht, col = morecolors(50),
                     show_row_names = F, show_column_names = F,
                     name = "distance")
myheatmap

## -----------------------------------------------------------------------------
# 2D plot for Dim2 and Dim3 [With V update]
sel_colors <- c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "mediumturquoise",
                "#E6AB02","darkmagenta", "#666666", "black","darkolivegreen2")
i <- 4
j <- 7
ggplot2::ggplot(df_proj, aes(x = df_proj[, i], y = df_proj[, j], col = tissue,
                              shape = species, fill = tissue)) +  
  xlab(paste("OSBF Dim", i)) + ylab(paste("OSBF Dim", j)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) + 
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme() +
  theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# function to compute Tau
calc_tissue_specificity <- function(a) {
    a <- as.matrix(a)
    b <- a / matrixStats::rowMaxs(a)
    return(rowSums(1 - b) / (ncol(b) - 1))
}
    
Tau <- lapply(avg_counts, function(x) { calc_tissue_specificity(x)})
avg_counts_scaled <- lapply(avg_counts, function(x) { t(scale(t(x)))})

combine_expr <- list()
for (sp in names(avg_counts_scaled)) {
  x <- as.data.frame(avg_counts_scaled[[sp]])
  x[["Tau"]] <- Tau[[sp]]
  combine_expr[[sp]] <- x
}

## -----------------------------------------------------------------------------
sel_dim <- 4
sel_tissue <- "bcell"
species <- "Homo_sapiens"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = T]
expr1 <- expr[, c(paste0(getSpeciesShortName(species), "_", sel_tissue),
                   "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")
head(expr1)

## -----------------------------------------------------------------------------
# plot scatter
mid <- 0
p1 <- ggplot2::ggplot(expr1, aes(x = Tau, y = coef, col = tissue_zscore)) + theme_bw() +
  geom_point(size = 0.5) + xlab("Expression specificity") +
  ylab(paste0("Dim", sel_dim, " Coefficient")) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red", space = "Lab") +
  scale_y_continuous(limits = c(-1 * max(abs(expr1$coef)), max(abs(expr1$coef))),
                     breaks = seq(-1 * round(max(abs(expr$coef)), 2),
                                  round(max(abs(expr$coef)), 2), by = 0.01)) +
  customTheme() +  theme(legend.position = "right",
                         legend.direction = "vertical") +
  labs(title = getScientificName(species), color = "Z-score") +
  theme(legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(face = "italic"))
p1

## -----------------------------------------------------------------------------
sel_dim <- 7
sel_tissue <- "nkcell"
species <- "Homo_sapiens"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = T]
expr1 <- expr[, c(paste0(getSpeciesShortName(species), "_", sel_tissue),
                   "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")

# plot scatter
mid <- 0
p1 <- ggplot2::ggplot(expr1, aes(x = Tau, y = coef, col = tissue_zscore)) + theme_bw() +
  geom_point(size = 0.5) + xlab("Expression specificity") +
  ylab(paste0("Dim", sel_dim, " Coefficient")) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red", space = "Lab") +
  scale_y_continuous(limits = c(-1 * max(abs(expr1$coef)), max(abs(expr1$coef))),
                     breaks = seq(-1 * round(max(abs(expr$coef)), 2),
                                  round(max(abs(expr$coef)), 2), by = 0.01)) +
  customTheme() +  theme(legend.position = "right",
                         legend.direction = "vertical") +
  labs(title = getScientificName(species), color = "Z-score") +
  theme(legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(face = "italic"))
p1

## -----------------------------------------------------------------------------
# install packages
pkgs <- c("goseq")
require_install <- pkgs[!(pkgs %in% row.names(installed.packages()))]
if (length(require_install)) {
  if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("goseq")
}
suppressPackageStartupMessages({
  library(goseq)
})

## -----------------------------------------------------------------------------
sel_dim <- 4
sel_tissue <- "bcell"
top_genes <- 100
# axis positive (pos) or negative (neg)
sel_sign <- "pos"

## -----------------------------------------------------------------------------
species <- "Homo_sapiens"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = T]
expr1 <- expr[, c(paste0(getSpeciesShortName(species), "_", sel_tissue),
                         "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")
if (sel_sign == "neg") {
  cat("\n selecting negative loadings")
  expr1_selsign <- expr1[expr1$coef < 0, ]
  expr1_bgsign <- expr1[expr1$coef >= 0, ]
} else {
  cat("\n selecting positive loadings")
  expr1_selsign <- expr1[expr1$coef >= 0, ]
  expr1_bgsign <- expr1[expr1$coef < 0, ]
}
expr1_selsign$score <- expr1_selsign$Tau * abs(expr1_selsign$coef)
expr1_selsign$rank <- rank(-1 * expr1_selsign$score)
expr1_selsign <- expr1_selsign[order(expr1_selsign$rank), ]
# gene list of interest
genes_fg <- row.names(expr1_selsign[expr1_selsign$rank <= top_genes, ])
# background genes
# For GO analysis, we will use genes with opposite sign loadings as
# the background.
genes_bg <- row.names(expr1_bgsign)

genes_bg <- genes_bg[!genes_bg %in% genes_fg]
genome <- "hg38"
total_genes <- unique(c(genes_fg, genes_bg))
up_genes <- as.integer(total_genes %in%  genes_fg)
names(up_genes) <- total_genes

## ---- echo = FALSE, warning = FALSE-------------------------------------------
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
# load("hg38 length data")
load(paste0(path, "GOKeggFiles/hg38_length.EnsemblV94.RData"))
lengthData.up <- lengthData[names(up_genes)]
# load("hg38 EnsembleID to GO data")
load(paste0(path, "GOKeggFiles/hg38_geneID2GO.EnsemblID2GO.EnsembleV94.Robj"))
pwf <- goseq::nullp(up_genes, bias.data = lengthData.up, plot.fit = F)
go <- goseq::goseq(pwf, "hg38", "ensGene", gene2cat = geneID2GO, test.cats = c("GO:BP"))
go.sub <- go[go$ontology == "BP", ]
go.sub$padj <- p.adjust(go.sub$over_represented_pvalue, method = "BH")
go.sub[["ratio"]] <- round(go.sub[["numDEInCat"]] / go.sub[["numInCat"]], 4)
go.sub <- go.sub[with(go.sub, order(padj, decreasing = c(FALSE))), ]
go.sub$over_represented_pvalue <- NULL
go.sub$under_represented_pvalue <- NULL

## -----------------------------------------------------------------------------
head(go.sub)

## -----------------------------------------------------------------------------
# GO enrichment plot for human
go_out <- head(go.sub, n = 8)
go_out$padj <- as.numeric(go_out$padj)
go_out$term <- factor(go_out$term, levels = go_out$term)
breaks <- round(c(0, 1 / 4, 2 / 4, 3 / 4, 1) * max(go_out[["ratio"]]), 2)
go_plot <- ggplot2::ggplot(go_out, aes(x = term, y = ratio, fill = padj)) + geom_col() +
  scale_y_continuous(expand = c(0, 0), breaks = breaks,
                     limits = c(0, max(go_out[["ratio"]] + 0.05))) +
  scale_x_discrete() + coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  ylab(paste0("Ratio of genes in GO category (",
              species, ")")) +
  xlab("") + customTheme() + theme(legend.position = "right",
                                   legend.direction = "vertical",
                                   plot.margin = unit(c(10, 5, 5, 5), "mm"))
go_plot

## -----------------------------------------------------------------------------
sel_dim <- 7
sel_tissue <- "nkcell"
top_genes <- 100
# axis positive (pos) or negative (neg)
sel_sign <- "pos"

species <- "Homo_sapiens"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = T]
expr1 <- expr[, c(paste0(getSpeciesShortName(species), "_", sel_tissue), "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")
if (sel_sign == "neg") {
  cat("\n selecting negative loadings")
  expr1_selsign <- expr1[expr1$coef < 0, ]
  expr1_bgsign <- expr1[expr1$coef >= 0, ]
} else {
  cat("\n selecting positive loadings")
  expr1_selsign <- expr1[expr1$coef >= 0, ]
  expr1_bgsign <- expr1[expr1$coef < 0, ]
}
expr1_selsign$score <- expr1_selsign$Tau * abs(expr1_selsign$coef)
expr1_selsign$rank <- rank(-1 * expr1_selsign$score)
expr1_selsign <- expr1_selsign[order(expr1_selsign$rank), ]
# gene list of interest
genes_fg <- row.names(expr1_selsign[expr1_selsign$rank <= top_genes, ])
# background genes
# For GO analysis, we will use genes with opposite sign loadings as
# the background.
genes_bg <- row.names(expr1_bgsign)

genes_bg <- genes_bg[!genes_bg %in% genes_fg]
genome <- "hg38"
total_genes <- unique(c(genes_fg, genes_bg))
up_genes <- as.integer(total_genes %in%  genes_fg)
names(up_genes) <- total_genes

## ---- echo = FALSE, warning = FALSE-------------------------------------------
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
# load("hg38 length data")
load(paste0(path, "GOKeggFiles/hg38_length.EnsemblV94.RData"))
lengthData.up <- lengthData[names(up_genes)]
# load("hg38 EnsembleID to GO data")
load(paste0(path, "GOKeggFiles/hg38_geneID2GO.EnsemblID2GO.EnsembleV94.Robj"))
pwf <- goseq::nullp(up_genes, bias.data = lengthData.up, plot.fit = F)
go <- goseq::goseq(pwf, "hg38", "ensGene", gene2cat = geneID2GO, test.cats = c("GO:BP"))
go.sub <- go[go$ontology == "BP", ]
go.sub$padj <- p.adjust(go.sub$over_represented_pvalue, method = "BH")
go.sub[["ratio"]] <- round(go.sub[["numDEInCat"]] / go.sub[["numInCat"]], 4)
go.sub <- go.sub[with(go.sub, order(padj, decreasing = c(FALSE))), ]
go.sub$over_represented_pvalue <- NULL
go.sub$under_represented_pvalue <- NULL

## -----------------------------------------------------------------------------
head(go.sub)

## -----------------------------------------------------------------------------
# GO enrichment plot for human
go_out <- head(go.sub, n = 8)
go_out$padj <- as.numeric(go_out$padj)
go_out$term <- factor(go_out$term, levels = go_out$term)
breaks <- round(c(0, 1 / 4, 2 / 4, 3 / 4, 1) * max(go_out[["ratio"]]), 2)
go_plot <- ggplot2::ggplot(go_out, aes(x = term, y = ratio, fill = padj)) + geom_col() +
  scale_y_continuous(expand = c(0, 0), breaks = breaks,
                     limits = c(0, max(go_out[["ratio"]] + 0.05))) +
  scale_x_discrete() + coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  ylab(paste0("Ratio of genes in GO category (",
              species, ")")) +
  xlab("") + customTheme() + theme(legend.position = "right",
                                   legend.direction = "vertical",
                                   plot.margin = unit(c(10, 5, 5, 5), "mm"))
go_plot

## -----------------------------------------------------------------------------
# set seed
s1 <- 32
s2 <- 135
species <- c("Homo_sapiens", "Mus_musculus")
species_short <- sapply(species, getSpeciesShortName)
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
counts_list_shuff <- metadata_list_shuff <- avg_counts_shuff <- list()
for (sp in species) {
  # reading raw counts
  counts <- read.table(paste0(path, "human_mouse_blood_counts/", sp,
                              "_blood_rawcounts.tsv"), header = T,
                           sep = "\t", row.names = 1)
  info <- tstrsplit(colnames(counts), "_")
  metadata <- data.frame(project = info[[1]],
        species = info[[2]],
        tissue = info[[3]],
        gsm = info[[4]],
        name = colnames(counts),
        stringsAsFactors = F)
  metadata$ref <- seq_len(nrow(metadata))
  metadata$key <- paste0(metadata$species, "_", metadata$tissue)
  metadata$tissue_factor <- factor(metadata$tissue)
  
  counts_avg <- calcAvgCounts(counts, metadata)
  cnames <- colnames(counts_avg)
  rnames <- row.names(counts_avg)
  set.seed(s1)
  counts_avg <- as.data.frame(apply(counts_avg, 2, sample))
  set.seed(s2)
  counts_avg <- as.data.frame(t(apply(counts_avg, 1, sample)))
  colnames(counts_avg) <- cnames
  row.names(counts_avg) <- rnames
  # normalize the shuffled counts to log TPM
  # set the path to the working directory. Change this accordingly
  path <- "~/Dropbox/0.Analysis/0.paper/"
  gene_length <- read.table(paste0(path, "ensembl94_annotation/", sp, "_genelength.tsv"),
                        sep = "\t",header = T, row.names = 1, stringsAsFactors = F)
  if (!all(row.names(counts_avg) %in% row.names(gene_length))) stop("Error")
  gene_length$Length <- gene_length$Length/1e3
  gene_length <- gene_length[row.names(counts_avg), , drop=T]
  names(gene_length) <- row.names(counts_avg)
  counts_tpm <- normalizeTPM(rawCounts = counts_avg, gene_len = gene_length)
  min_tpm <- 1
  counts_tpm[counts_tpm < min_tpm] <- 1
  counts_tpm <- log2(counts_tpm)

  info <- tstrsplit(colnames(counts_tpm),"_")
  metadata <- data.frame(
        species = info[[1]],
        tissue = info[[2]],
        name = colnames(counts_tpm),
        stringsAsFactors = F)
  metadata$key <- paste0(metadata$species,"_",metadata$tissue)
  avg_counts_shuff[[sp]] <- calcAvgCounts(counts_tpm, metadata)
  # write.table(avg_counts_shuff[[sp]],file=paste0(mainDir, "shuffledavgCounts/avg_logTPM_",
  #                               species,"_SHUFF_counts_seed",z,"_",z1,".tsv"),
  #           sep="\t",quote=F,col.names=NA)
  counts_list_shuff[[sp]] <- counts_tpm
  metadata_list_shuff[[sp]] <- metadata
}
# dims
sapply(counts_list_shuff, dim)
# remove zero counts
avg_counts_shuff <- lapply(avg_counts_shuff, removeZeros)
sapply(avg_counts_shuff, dim)

## -----------------------------------------------------------------------------
cat(format(Sys.time(), "%a %b %d %X %Y"),"\n")
osbf_shuf <- SBF(avg_counts_shuff, transform_matrix = TRUE, orthogonal = TRUE,
                 tol = 1e-2)
cat(format(Sys.time(), "%a %b %d %X %Y"),"\n")
osbf_shuf$error

## -----------------------------------------------------------------------------
Tau_null <- lapply(avg_counts_shuff, function(x) {calc_tissue_specificity(x)})
avg_counts_shuff_scaled <- lapply(avg_counts_shuff, function(x){ t(scale(t(x)))})
combine_expr_null <- list()
for (sp in names(avg_counts_shuff_scaled)) {
  x <- as.data.frame(avg_counts_shuff_scaled[[sp]])
  x[["Tau"]] <- Tau_null[[sp]]
  combine_expr_null[[sp]] <- x
}

## -----------------------------------------------------------------------------
sel_dim <- 4
sel_tissue <- "bcell"
# axis positive (pos) or negative (neg)
sel_sign <- "pos"
species <- "Homo_sapiens"
species_short <- "hsapiens"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = T]
expr1 <- expr[, c(paste0(species_short, "_", sel_tissue),
                   "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")

# null loadings for the same dimensions
expr_null <- combine_expr_null[[species]]
null_u <- osbf_shuf$u[[species]]
expr_null[["coef"]] <- null_u[, sel_dim, drop = T]
expr1_null <- expr_null[, c(paste0(species_short,"_",sel_tissue),"Tau","coef")]
if (sel_sign == "pos") {
  expr1 <- expr1[expr1$coef >= 0, ]
  expr1_null <- expr1_null[expr1_null$coef >= 0, ]
} else if (sel_dim == "neg") {
  expr1 <- expr1[expr1$coef < 0, ]
  expr1_null <- expr1_null[expr1_null$coef < 0, ]
}
expr1$score <- expr1$Tau * abs(expr1$coef)
expr1$rank <- rank(-1 * expr1$score)
expr1 <- expr1[order(expr1$rank), ]

expr1_null$score <- expr1_null$Tau * abs(expr1_null$coef)
expr1$pvalue <- sapply(expr1$score, function(x) { sum(as.integer(expr1_null$score > x))/length(expr1_null$score)})
head(expr1)

## -----------------------------------------------------------------------------
# cut off for the p-value
alpha <- 1e-3
summary(expr1$pvalue <= alpha)

## -----------------------------------------------------------------------------
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
gene_info <- read.table(paste0(path, "ensembl94_annotation/", species_short, "_genes_completeinfo.tsv"),
                               sep = "\t", header = T, quote = "\"")
gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id), ]
gene_info <- gene_info[gene_info$ensembl_gene_id %in% row.names(expr1), ]
row.names(gene_info) <- gene_info$ensembl_gene_id
gene_info <- gene_info[row.names(expr1), ]
expr1$gene_name <- gene_info$external_gene_name
expr1$biotype <- gene_info$gene_biotype
head(expr1, n = 10)

## -----------------------------------------------------------------------------
sel_dim <- 7
sel_tissue <- "nkcell"
# axis positive (pos) or negative (neg)
sel_sign <- "pos"
species <- "Mus_musculus"
species_short <- "mmusculus"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = T]
expr1 <- expr[, c(paste0(species_short, "_", sel_tissue),
                   "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")

# null loadings for the same dimensions
expr_null <- combine_expr_null[[species]]
null_u <- osbf_shuf$u[[species]]
expr_null[["coef"]] <- null_u[, sel_dim, drop = T]
expr1_null <- expr_null[, c(paste0(species_short,"_",sel_tissue),"Tau","coef")]
if (sel_sign == "pos") {
  expr1 <- expr1[expr1$coef >= 0, ]
  expr1_null <- expr1_null[expr1_null$coef >= 0, ]
} else if (sel_dim == "neg") {
  expr1 <- expr1[expr1$coef < 0, ]
  expr1_null <- expr1_null[expr1_null$coef < 0, ]
}
expr1$score <- expr1$Tau * abs(expr1$coef)
expr1$rank <- rank(-1 * expr1$score)
expr1 <- expr1[order(expr1$rank), ]

expr1_null$score <- expr1_null$Tau * abs(expr1_null$coef)
expr1$pvalue <- sapply(expr1$score, function(x) { sum(as.integer(expr1_null$score > x))/length(expr1_null$score)})
#head(expr1)

## -----------------------------------------------------------------------------
# cut off for the p-value
alpha <- 1e-3
summary(expr1$pvalue <= alpha)

## -----------------------------------------------------------------------------
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
gene_info <- read.table(paste0(path, "ensembl94_annotation/", species_short,
                               "_genes_completeinfo.tsv"),
                               sep = "\t", header = T, quote = "\"")
gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id), ]
gene_info <- gene_info[gene_info$ensembl_gene_id %in% row.names(expr1), ]
row.names(gene_info) <- gene_info$ensembl_gene_id
gene_info <- gene_info[row.names(expr1), ]
expr1$gene_name <- gene_info$external_gene_name
expr1$biotype <- gene_info$gene_biotype
head(expr1, n = 10)

## -----------------------------------------------------------------------------
sessionInfo()

