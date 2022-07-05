## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
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
species <- c("Homo_sapiens", "Pan_troglodytes", "Macaca_mulatta",
             "Mus_musculus", "Rattus_norvegicus", "Bos_taurus",
             "Sus_scrofa", "Gallus_gallus")
species_short <- sapply(species, getSpeciesShortName)
species_short
# common tissues present in all 8 species
tissues <-  c("brain", "heart", "kidney", "liver", "lung", "testis")

## -----------------------------------------------------------------------------
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
counts_list <- metadata_list <- list()
require(data.table)
for (sp in species) {
  # read logTPM counts for each species
  counts <- data.table::fread(paste0(path, "counts/", sp, "_logTPM.tsv"),
                              sep = "\t", header = TRUE, data.table = FALSE,
                              nThread = 4)
  row.names(counts) <- counts$V1
  counts$V1 <- NULL
  col_fields <- data.table::tstrsplit(colnames(counts), "_")
  metadata <- data.frame(
      project = col_fields[[1]],
      species = col_fields[[2]],
      tissue = col_fields[[3]],
      gsm = col_fields[[4]],
      name = colnames(counts),
      stringsAsFactors = FALSE)
  metadata_sel <- metadata[metadata$tissue %in% tissues, , drop = FALSE]
  counts_sel <- counts[, colnames(counts) %in% metadata_sel$name, drop = FALSE]
  metadata_sel$ref <- seq_len(nrow(metadata_sel))
  metadata_sel$key <- paste0(metadata_sel$species, "_", metadata_sel$tissue)

  counts_list[[sp]] <- counts_sel
  metadata_list[[sp]] <- metadata_sel
}

## -----------------------------------------------------------------------------
sapply(counts_list, dim)

## -----------------------------------------------------------------------------
avg_counts <- list()
for (sp in species) {
    avg_counts[[sp]] <- calcAvgCounts(counts_list[[sp]],
                                      metadata_list[[sp]])
}
# check tissue columns are matching in each species
c_tissues <- as.data.frame(sapply(avg_counts, function(x) {
            data.table::tstrsplit(colnames(x), "_")[[2]]
  }))
if (!all(apply(c_tissues, 1, function(x) all(x == x[1])))) {
        stop("Error! tissues not matching")
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
                                                 drop = FALSE]
}

## -----------------------------------------------------------------------------
t1 <- proc.time()
# decrease tol to minimize error
osbf <- SBF(avg_counts, transform_matrix = TRUE, orthogonal = TRUE,
                     optimizeV = TRUE, tol = 1e-3)#, tol = 1e-10)
t2 <- proc.time()
cat("Time taken:\n")
t2 - t1

## -----------------------------------------------------------------------------
hogsvd <- HOGSVD(avg_counts)

## -----------------------------------------------------------------------------
names(hogsvd)

## -----------------------------------------------------------------------------
calcDecompError(avg_counts, osbf$u, osbf$delta, osbf$v)

## -----------------------------------------------------------------------------
calcDecompError(avg_counts, hogsvd$u, hogsvd$delta, hogsvd$v)

## -----------------------------------------------------------------------------
zapsmall(osbf$v %*% t(osbf$v))

## -----------------------------------------------------------------------------
zapsmall(t(hogsvd$v) %*% hogsvd$v)

## -----------------------------------------------------------------------------
# get u factor for fist species
osbf_u <- osbf$u[[names(osbf$u)[1]]]
zapsmall(t(osbf_u) %*% osbf_u)

## -----------------------------------------------------------------------------
# get u factor for fist species
hogsvd_u <- hogsvd$u[[names(hogsvd$u)[1]]]
zapsmall(t(hogsvd_u) %*% hogsvd_u)

## -----------------------------------------------------------------------------
percentInfo_osbf <- calcPercentInfo(osbf)
for (i in names(osbf$delta)) {
  cat("\n", sprintf("%-25s:", i), sprintf("%8.2f", percentInfo_osbf[[i]]))
}

## -----------------------------------------------------------------------------
percentInfo_hogsvd <- calcPercentInfo(hogsvd)
for (i in names(hogsvd$delta)) {
  cat("\n", sprintf("%-25s:", i), sprintf("%8.2f", percentInfo_hogsvd[[i]]))
}

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
sel_dim <- 1
species <- "Homo_sapiens"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = TRUE]

# plot scatter
mid <- 0.5
p1 <- ggplot2::ggplot(expr, aes(x = Tau, y = coef, col = Tau)) +
  theme_bw() +
  geom_point(size = 0.5) + xlab("Expression specificity") +
  ylab(paste0("OSBF Dim", sel_dim, " Coefficient")) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red", space = "Lab") +
  customTheme() +  theme(legend.position = "right",
                         legend.direction = "vertical") +
  labs(title = getScientificName(species), color = "Tau") +
  theme(legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(face = "italic"))
p1

## -----------------------------------------------------------------------------
sel_dim <- 3
species <- "Homo_sapiens"
expr <- combine_expr[[species]]
hogsvd_coef <- hogsvd$u[[species]]
expr[["coef"]] <- hogsvd_coef[, sel_dim, drop = TRUE]

# plot scatter
mid <- 0.5
p1 <- ggplot2::ggplot(expr, aes(x = Tau, y = coef, col = Tau)) +
  theme_bw() +
  geom_point(size = 0.5) + xlab("Expression specificity") +
  ylab(paste0("HO GSVD Dim", sel_dim, " Coefficient")) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red", space = "Lab") +
  customTheme() +  theme(legend.position = "right",
                         legend.direction = "vertical") +
  labs(title = getScientificName(species), color = "Tau") +
  theme(legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(face = "italic"))
p1

## -----------------------------------------------------------------------------
df_proj <- projectCounts(counts_list_sub, osbf)
meta1 <- data.table::tstrsplit(row.names(df_proj), "_")
df_proj$tissue <- factor(meta1[[3]])
df_proj$species <- factor(meta1[[2]])
df_proj <- df_proj %>% mutate(species = factor(species,
                            levels = species_short))

## -----------------------------------------------------------------------------
ho_proj <- projectCounts(counts_list_sub, hogsvd)
meta1 <- data.table::tstrsplit(row.names(ho_proj), "_")
ho_proj$tissue <- factor(meta1[[3]])
ho_proj$species <- factor(meta1[[2]])
ho_proj <- ho_proj %>% mutate(species = factor(species,
                            levels = species_short))

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
ht <- ComplexHeatmap::HeatmapAnnotation(tissue = meta[[3]], species = meta[[2]],
                        col = list(tissue = c("brain" = "#1B9E77",
                                              "heart" = "#D95F02",
                                              "kidney" = "#7570B3",
                                              "liver" = "#E7298A",
                                              "lung" = "#66A61E",
                                              "testis" = "#E6AB02"),
                                   species = c("hsapiens" = "#66C2A5",
                                               "ptroglodytes" = "#FC8D62",
                                               "mmulatta" = "#8DA0CB",
                                               "mmusculus" = "#E78AC3",
                                               "rnorvegicus" = "#A6D854",
                                               "btaurus" = "#FFD92F",
                                               "sscrofa" = "#E5C494",
                                               "ggallus" = "#B3B3B3")),
                        #show_annotation_name = F,
                        annotation_name_side = "left")
mypalette <- RColorBrewer::brewer.pal(9, "Blues")
morecolors <- colorRampPalette(mypalette)

myheatmap <- ComplexHeatmap::Heatmap(as.matrix(data_dist), cluster_rows = TRUE,
                                     clustering_method_rows = "centroid",
                                     cluster_columns = TRUE,
                                     clustering_method_columns = "centroid",
                                     top_annotation = ht, col = morecolors(50),
                                     show_row_names = FALSE,
                                     show_column_names = FALSE,
                                     name = "distance")
myheatmap

## -----------------------------------------------------------------------------
data <- ho_proj
data$tissue <- NULL
data$species <- NULL
data <- as.matrix(data)
data_dist <- as.matrix(dist(data, method = "euclidean"))
meta <- data.table::tstrsplit(colnames(data_dist), "_")
ht <- ComplexHeatmap::HeatmapAnnotation(tissue = meta[[3]], species = meta[[2]],
                        col = list(tissue = c("brain" = "#1B9E77",
                                              "heart" = "#D95F02",
                                              "kidney" = "#7570B3",
                                              "liver" = "#E7298A",
                                              "lung" = "#66A61E",
                                              "testis" = "#E6AB02"),
                                   species = c("hsapiens" = "#66C2A5",
                                               "ptroglodytes" = "#FC8D62",
                                               "mmulatta" = "#8DA0CB",
                                               "mmusculus" = "#E78AC3",
                                               "rnorvegicus" = "#A6D854",
                                               "btaurus" = "#FFD92F",
                                               "sscrofa" = "#E5C494",
                                               "ggallus" = "#B3B3B3")),
                        #show_annotation_name = F,
                        annotation_name_side = "left")
myheatmap <- ComplexHeatmap::Heatmap(as.matrix(data_dist), cluster_rows = TRUE,
                                     clustering_method_rows = "centroid",
                                     cluster_columns = TRUE,
                                     clustering_method_columns = "centroid",
                                     top_annotation = ht, col = morecolors(50),
                                     show_row_names = FALSE,
                                     show_column_names = FALSE,
                                     name = "distance")
myheatmap

## -----------------------------------------------------------------------------
# install packages
pkgs <- c("gridExtra")
require_install <- pkgs[!(pkgs %in% row.names(installed.packages()))]
if (length(require_install))
  install.packages(require_install)
suppressPackageStartupMessages({
  library(gridExtra)
})

## -----------------------------------------------------------------------------
sel_colors <- RColorBrewer::brewer.pal(8, "Dark2")[seq_len(length(unique(df_proj$tissue)))]

p_list <- lapply(c(1, 3:6), function(k) {
  ggplot2::ggplot(df_proj, aes(x = df_proj[, 2], y = df_proj[, k], col = tissue,
                             shape = species, fill = tissue)) +
  xlab(paste("OSBF Dim", 2)) + ylab(paste("OSBF Dim", k)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) +
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme(base_size = 6) +
  theme(legend.title = element_blank())
})
h_list <- lapply(c(1, 3:6), function(k) {
  ggplot2::ggplot(ho_proj, aes(x = ho_proj[, 2], y = ho_proj[, k], col = tissue,
                             shape = species, fill = tissue)) +
  xlab(paste("HO GSVD Dim", 2)) + ylab(paste("HO GSVD Dim", k)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) +
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme(base_size = 6) +
  theme(legend.title = element_blank())
})

## -----------------------------------------------------------------------------
grid.arrange(p_list[[1]], h_list[[1]], ncol = 2)

## -----------------------------------------------------------------------------
grid.arrange(p_list[[2]], h_list[[2]], p_list[[3]], h_list[[3]], ncol = 2)

## -----------------------------------------------------------------------------
grid.arrange(p_list[[4]], h_list[[4]], p_list[[5]], h_list[[5]], ncol = 2)

## -----------------------------------------------------------------------------
sessionInfo()

