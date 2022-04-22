## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

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
species <- c("Homo_sapiens", "Pan_troglodytes", "Macaca_mulatta",
             "Mus_musculus", "Rattus_norvegicus", "Bos_taurus",
             "Sus_scrofa", "Gallus_gallus")
# function to create short name for species
getSpeciesShortName <- function(species_list) {
    require(data.table)
    shortname <- c()
    for (sp in species_list) {
       if (!grepl("_", sp)) {
          cat("\nspecies:", sp, "\n")
          stop("Invalid species full name")
       }
       shortname <- c(shortname,
                       paste0(tolower(substr(tstrsplit(sp, "_")[[1]], 1, 1)),
                       tstrsplit(sp, "_")[[2]]))
    }
    return(shortname)
}
species_short <- sapply(species, getSpeciesShortName)
species_short
# common tissues present in all 8 species
tissues <-  c("brain", "heart", "kidney", "liver", "lung", "testis")

## -----------------------------------------------------------------------------
# set the path to the counts. Change it accordingly
counts_path <- "~/Dropbox/0.Analysis/0.paper/counts/"
counts_list <- metadata_list <- list()
require(data.table)
for (sp in species) {
  counts <- fread(paste0(counts_path, sp, "_logTPM.tsv"), sep = "\t",
                  header = T, data.table = FALSE, nThread = 4)
  row.names(counts) <- counts$V1
  counts$V1 <- NULL
  col_fields <- tstrsplit(colnames(counts), "_")
  metadata <- data.frame(
      project = col_fields[[1]],
      species = col_fields[[2]],
      tissue = col_fields[[3]],
      gsm = col_fields[[4]],
      name = colnames(counts),
      stringsAsFactors = F)
  metadata_sel <- metadata[metadata$tissue %in% tissues, , drop = F]
  counts_sel <- counts[, colnames(counts) %in% metadata_sel$name, drop = F]
  metadata_sel$ref <- seq_len(nrow(metadata_sel))
  metadata_sel$key <- paste0(metadata_sel$species, "_", metadata_sel$tissue)

  counts_list[[sp]] <- counts_sel
  metadata_list[[sp]] <- metadata_sel
}

## -----------------------------------------------------------------------------
sapply(counts_list, dim)

## -----------------------------------------------------------------------------
#function to calculate mean expression values for each tissue/group
calcAvgCounts <- function(counts, metadata, ndecimal = 4) {
   require(dplyr)
   counts.avg <- list()
   if ("key" %in% colnames(metadata)) {
     for (species.organ in sort(unique(metadata$key))) {
        no.of.libs <- length(colnames(counts)[grepl(species.organ,
                                                    colnames(counts))])
        if (no.of.libs > 1)
           counts.avg[[species.organ]] <- round(rowMeans(counts[, colnames(counts)[
                                   grepl(species.organ, colnames(counts))]],
                                   dim = 1), ndecimal)
        else if (no.of.libs == 1)
           counts.avg[[species.organ]] <- round(counts[, colnames(counts)[grepl(
                              species.organ, colnames(counts))], drop = T], ndecimal)
     }
     avg.counts <- as.data.frame(counts.avg %>% bind_cols())
     row.names(avg.counts) <- row.names(counts)
     return(avg.counts)
   } else {
     stop("'key' column not found in the metadata!Exiting!")
   }
}

avg_counts <- list()
for (sp in species) {
    avg_counts[[sp]] <- calcAvgCounts(counts_list[[sp]],
                                      metadata_list[[sp]])
}
# check tissue columns are matching in each species
c_tissues <- as.data.frame(sapply(avg_counts, function(x) {
            tstrsplit(colnames(x), "_")[[2]]
  }))
if (!all(apply(c_tissues, 1, function(x) all(x == x[1])))) {
        stop("Error! tissues not matching")
}

## -----------------------------------------------------------------------------
sapply(avg_counts, dim)

## -----------------------------------------------------------------------------
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
cat("\n=======================================================")
cat("\nASBF with optimizing V = FALSE started\n")
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n")
t1 <- proc.time()
sbf_noVupdate <- SBF(avg_counts, transform_matrix = TRUE, approximate = TRUE,
                     optimizeV = FALSE, tol = 1e-3)#, tol = 1e-10)
t2 <- proc.time()
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n")
cat("ASBF with optimizing V = FALSE finished\n")
cat("Time taken:\n")
t2 - t1

## -----------------------------------------------------------------------------
cat("\n=======================================================")
cat("\nASBF with optimizing V=TRUE started\n")
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n")
t1 <- proc.time()
sbf <- SBF(avg_counts, transform_matrix = TRUE, approximate = TRUE,
                     optimizeV = TRUE, tol = 1e-3)#, tol = 1e-10)
t2 <- proc.time()
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n")
cat("ASBF with optimizing V=TRUE finished\n")
cat("Time taken:\n")
t2 - t1

## -----------------------------------------------------------------------------
cat("\n", sprintf("%-27s:", "Final error [No V update]"), sprintf("%16.2f",
                                                        sbf_noVupdate$error))
cat("\n", sprintf("%-27s:", "Final error [With V update]"), sprintf("%16.2f",
                                                                    sbf$error))
cat("\n", sprintf("%-27s:", "# of update [No V update]"), sprintf("%16d",
                                                      sbf_noVupdate$error_pos))
cat("\n", sprintf("%-27s:", "# of update [With V update]"), sprintf("%16d",
                                                                sbf$error_pos))

## -----------------------------------------------------------------------------
sbf_noVupdate$error / sbf$error

## -----------------------------------------------------------------------------
par(mfrow = c(1, 3))
plot(x = seq_len(length(sbf$error_vec)), y = sbf$error_vec,
     xlab = "# of updates (step 1 -> )",
     ylab = "Factorization error", col = "red")
plot(x = 10:length(sbf$error_vec),
     y = sbf$error_vec[10:length(sbf$error_vec)],
     xlab = "# of updates (step 10 -> )",
     ylab = "Factorization error", col = "red")
plot(x = 50:length(sbf$error_vec),
     y = sbf$error_vec[50:length(sbf$error_vec)],
     xlab = "# of updates (step 50 -> )",
     ylab = "Factorization error", col = "red")

## -----------------------------------------------------------------------------
percentInfo_noVupdate <- lapply(sbf_noVupdate$d, function(x) {
  round(x^2 / sum(x^2) * 100, 2)})
cat("\nPercentage for each delta [No V update]:")
for (i in names(sbf_noVupdate$d)) {
  cat("\n", sprintf("%-25s:", i), sprintf("%8.2f", percentInfo_noVupdate[[i]]))
}
cat("\n-------------------------------------------------------------------------")


percentInfo <- lapply(sbf$d, function(x) {
  round(x^2 / sum(x^2) * 100, 2) })
cat("\nPercentage for each delta:")
for (i in names(sbf$d)) {
  cat("\n", sprintf("%-25s:", i), sprintf("%8.2f", percentInfo[[i]]))
}

## -----------------------------------------------------------------------------
d_inv_NoVupdate <- list()
for (sp in names(avg_counts)) {
  if (length(sbf_noVupdate$d[[sp]]) == 1) {
    d_inv_NoVupdate[[sp]] <- as.matrix(diag(as.matrix(1 / sbf_noVupdate$d[[sp]])))
    } else {
      d_inv_NoVupdate[[sp]] <- as.matrix(diag(1 / sbf_noVupdate$d[[sp]]))
    }
}
d_inv <- list()
for (sp in names(avg_counts)) {
  if (length(sbf$d[[sp]]) == 1) {
    d_inv[[sp]] <- as.matrix(diag(as.matrix(1 / sbf$d[[sp]])))
    } else {
      d_inv[[sp]] <- as.matrix(diag(1 / sbf$d[[sp]]))
    }
}

## -----------------------------------------------------------------------------
# project using no V update estimates
avgcounts_projected_noVupdate  <- counts_projected_noVupdate <- list()
n_noVupdate <- n1_noVupdate <- c()
for (sp in names(avg_counts)) {
  avgcounts_projected_noVupdate[[sp]] <- as.matrix(t(avg_counts[[sp]])) %*%
    as.matrix(sbf_noVupdate$u[[sp]]) %*% d_inv_NoVupdate[[sp]]
  n_noVupdate <- c(n_noVupdate, row.names(avgcounts_projected_noVupdate[[sp]]))
  counts_projected_noVupdate[[sp]] <- as.matrix(t(counts_list_sub[[sp]])) %*%
    as.matrix(sbf_noVupdate$u[[sp]]) %*% d_inv_NoVupdate[[sp]]
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
meta1 <- tstrsplit(row.names(df_proj_noVupdate), "_")
df_proj_noVupdate$tissue <- factor(meta1[[3]])
df_proj_noVupdate$species <- factor(meta1[[2]])
df_proj_noVupdate <- df_proj_noVupdate %>% mutate(species = factor(species,
                            levels = species_short))

# project using V update estimates
avgcounts_projected  <- counts_projected <- list()
n <- n1 <- c()
for (sp in names(avg_counts)) {
  avgcounts_projected[[sp]] <- as.matrix(t(avg_counts[[sp]])) %*%
    as.matrix(sbf$u[[sp]]) %*% d_inv[[sp]]
  n <- c(n, row.names(avgcounts_projected[[sp]]))
  counts_projected[[sp]] <- as.matrix(t(counts_list_sub[[sp]])) %*%
    as.matrix(sbf$u[[sp]]) %*% d_inv[[sp]]
  n1 <- c(n1, row.names(counts_projected[[sp]]))
}

# combine projected profiles
df_proj_avg <- as.data.frame(do.call(rbind, avgcounts_projected))
rownames(df_proj_avg) <- n
colnames(df_proj_avg) <- paste0("Dim", 1:ncol(df_proj_avg))
meta <- tstrsplit(row.names(df_proj_avg), "_")
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
pkgs <- c("ComplexHeatmap", "RColorBrewer")
require_install <- pkgs[!(pkgs %in% row.names(installed.packages()))]
if (length(require_install))
  install.packages(require_install)
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(RColorBrewer)
})

## -----------------------------------------------------------------------------
data_noVupdate <- df_proj_noVupdate
data_noVupdate$tissue <- NULL
data_noVupdate$species <- NULL
data_noVupdate <- as.matrix(data_noVupdate)
data_noVupdate_dist <- as.matrix(dist(data_noVupdate, method = "euclidean"))
meta <- tstrsplit(colnames(data_noVupdate_dist), "_")
ht <- HeatmapAnnotation(tissue = meta[[3]], species = meta[[2]],
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
mypalette <- brewer.pal(9, "Blues")
morecolors <- colorRampPalette(mypalette)

myheatmap1 <- Heatmap(as.matrix(data_noVupdate_dist), cluster_rows = T,
                     clustering_method_rows = "centroid",
                     cluster_columns = T,
                     clustering_method_columns = "centroid",
                     top_annotation = ht, col = morecolors(50),
                     show_row_names = F, show_column_names = F,
                     name = "distance")
myheatmap1

## -----------------------------------------------------------------------------
data <- df_proj
data$tissue <- NULL
data$species <- NULL
data <- as.matrix(data)
data_dist <- as.matrix(dist(data, method = "euclidean"))
meta <- tstrsplit(colnames(data_dist), "_")
ht <- HeatmapAnnotation(tissue = meta[[3]], species = meta[[2]],
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
mypalette <- brewer.pal(9, "Blues")
morecolors <- colorRampPalette(mypalette)

myheatmap <- Heatmap(as.matrix(data_dist), cluster_rows = T,
                     clustering_method_rows = "centroid",
                     cluster_columns = T,
                     clustering_method_columns = "centroid",
                     top_annotation = ht, col = morecolors(50),
                     show_row_names = F, show_column_names = F,
                     name = "distance")
myheatmap

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
customTheme <- function(base_size = 10, base_family="helvetica") {
   require(grid)
   require(ggthemes)
   (theme_foundation(base_size = base_size)
   + theme(plot.title = element_text(face = "bold",
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
# 2D plot for Dim1 and Dim2 [No V update]
sel_colors <- brewer.pal(8, "Dark2")[seq_len(length(unique(df_proj_noVupdate$tissue)))]
i <- 1
j <- 2
ggplot(df_proj_noVupdate, aes(x = df_proj_noVupdate[, i],
                              y = df_proj_noVupdate[, j], col = tissue,
                              shape = species, fill = tissue)) +
  xlab(paste("A-SBF Dim", i)) + ylab(paste("A-SBF Dim", j)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) +
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme() +
  theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# 2D plot for Dim1 and Dim2 [With V update]
sel_colors <- brewer.pal(8, "Dark2")[seq_len(length(unique(df_proj$tissue)))]
i <- 1
j <- 2
ggplot(df_proj, aes(x = df_proj[, i],
                              y = df_proj[, j], col = tissue,
                              shape = species, fill = tissue)) +
  xlab(paste("A-SBF Dim", i)) + ylab(paste("A-SBF Dim", j)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) + 
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme() +
  theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# 2D plot for Dim1 and Dim2 [No V update]
sel_colors <- brewer.pal(8, "Dark2")[seq_len(length(unique(df_proj_noVupdate$tissue)))]
i <- 2
j <- 3
ggplot(df_proj_noVupdate, aes(x = df_proj_noVupdate[, i],
                              y = df_proj_noVupdate[, j], col = tissue,
                              shape = species, fill = tissue)) +
  xlab(paste("A-SBF Dim", i)) + ylab(paste("A-SBF Dim", j)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) + 
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme() +
  theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# 2D plot for Dim1 and Dim2 [With V update]
sel_colors <- brewer.pal(8, "Dark2")[seq_len(length(unique(df_proj$tissue)))]
i <- 2
j <- 3
ggplot(df_proj, aes(x = df_proj[, i], y = df_proj[, j], col = tissue,
                              shape = species, fill = tissue)) +  
  xlab(paste("A-SBF Dim", i)) + ylab(paste("A-SBF Dim", j)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) + 
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme() +
  theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# function to compute Tau
calc_tissue_specificity <- function(a) {
    a <- as.matrix(a)
    b <- a / rowMaxs(a)
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
sel_dim <- 2
sel_tissue <- "testis"
species <- "Homo_sapiens"
expr <- combine_expr[[species]]
asbf_coef <- sbf$u[[species]]
expr[["coef"]] <- asbf_coef[, sel_dim, drop = T]
expr1 <- expr[, c(paste0(getSpeciesShortName(species), "_", sel_tissue),
                   "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")
head(expr1)

## -----------------------------------------------------------------------------
# function to return scientific name
species_scientific <- function(x) {
  parts <- tstrsplit(x, "_")
  return(paste0(substr(x, 1, 1), ".", parts[[2]]))
}

# plot scatter
mid <- 0
p1 <- ggplot(expr1, aes(x = Tau, y = coef, col = tissue_zscore)) + theme_bw() +
  geom_point(size = 0.5) + xlab("Expression specificity") +
  ylab(paste0("Dim", sel_dim, " Coefficient")) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red", space = "Lab") +
  scale_y_continuous(limits = c(-1 * max(abs(expr1$coef)), max(abs(expr1$coef))),
                     breaks = seq(-1 * round(max(abs(expr$coef)), 2),
                                  round(max(abs(expr$coef)), 2), by = 0.01)) +
  customTheme() +  theme(legend.position = "right",
                         legend.direction = "vertical") +
  labs(title = paste0(species_scientific(species)), color = "Z-score") +
  theme(legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(face = "italic"))
p1

## -----------------------------------------------------------------------------
sel_dim <- 2
sel_tissue <- "testis"
species <- "Sus_scrofa"
expr <- combine_expr[[species]]
asbf_coef <- sbf$u[[species]]
expr[["coef"]] <- asbf_coef[, sel_dim, drop = T]
expr1 <- expr[, c(paste0(getSpeciesShortName(species), "_", sel_tissue),
                   "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")

# plot scatter
mid <- 0
p2 <- ggplot(expr1, aes(x = Tau, y = coef, col = tissue_zscore)) + theme_bw() +
  geom_point(size = 0.5) + xlab("Expression specificity") +
  ylab(paste0("Dim", sel_dim, " Coefficient")) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red", space = "Lab") +
  scale_y_continuous(limits = c(-1 * max(abs(expr1$coef)), max(abs(expr1$coef))),
                     breaks = seq(-1 * round(max(abs(expr$coef)), 2),
                                  round(max(abs(expr$coef)), 2), by = 0.01)) +
  customTheme() +  theme(legend.position = "right",
                         legend.direction = "vertical") +
  labs(title = paste0(species_scientific(species)), color = "Z-score") +
  theme(legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(face = "italic"))
p2

## -----------------------------------------------------------------------------
# load GO analysis files
# go_path <- "path to shared 1-to-1 orthologs"
go_path <- "~/Dropbox/0.Analysis/0.paper/GOKeggFiles/"
file <- "allwayOrthologs_hsapiens-ptroglodytes-mmulatta-mmusculus-rnorvegicus-btaurus-sscrofa-ggallus_ens94.tsv"
orthologs <- read.table(file.path(go_path, file), header = T, sep = "\t")
head(orthologs)

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
sel_dim <- 2
sel_tissue <- "testis"
top_genes <- 100
sel_sign <- "neg"

## -----------------------------------------------------------------------------
species <- "Homo_sapiens"
expr <- combine_expr[[species]]
asbf_coef <- sbf$u[[species]]
expr[["coef"]] <- asbf_coef[, sel_dim, drop = T]
expr1 <- expr[, c(paste0(getSpeciesShortName(species), "_", sel_tissue), "Tau",
                  "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")
if (sel_sign == "neg") {
  cat("\n selecting negative loadings")
  expr1_selsign <- expr1[expr1$coef < 0, ]
} else {
  cat("\n selecting positive loadings")
  expr1_selsign <- expr1[expr1$coef >= 0, ]
}
expr1_selsign$score <- expr1_selsign$Tau * abs(expr1_selsign$coef)
expr1_selsign$rank <- rank(-1 * expr1_selsign$score)
expr1_selsign <- expr1_selsign[order(expr1_selsign$rank), ]
genes_fg <- row.names(expr1_selsign[expr1_selsign$rank <= top_genes, ])
if (species == "Homo_sapiens") {
  cat("\nUsing orthologs (human IDs) as background")
  genes_bg <- orthologs$hsapiens
}
if (species == "Mus_musculus") {
  cat("\nUsing orthologs (mouse IDs) as background")
  genes_bg <- orthologs$mmusculus
}
genes_bg <- genes_bg[!genes_bg %in% genes_fg]
genome <- "hg38"
total_genes <- unique(c(genes_fg, genes_bg))
up_genes <- as.integer(total_genes %in%  genes_fg)
names(up_genes) <- total_genes

## ---- echo = FALSE, warning = FALSE-------------------------------------------
# load("hg38 length data")
load("~/Dropbox/0.Analysis/0.paper/GOKeggFiles/hg38_length.EnsemblV94.RData")
lengthData.up <- lengthData[names(up_genes)]
# load("hg38 EnsembleID to GO data")
load("~/Dropbox/0.Analysis/0.paper/GOKeggFiles/hg38_geneID2GO.EnsemblID2GO.EnsembleV94.Robj")
pwf <- nullp(up_genes, bias.data = lengthData.up, plot.fit = F)
go <- goseq(pwf, "hg38", "ensGene", gene2cat = geneID2GO, test.cats = c("GO:BP"))
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
go_plot <- ggplot(go_out, aes(x = term, y = ratio, fill = padj)) + geom_col() +
  scale_y_continuous(expand = c(0, 0), breaks = breaks,
                     limits = c(0, max(go_out[["ratio"]] + 0.05))) +
  scale_x_discrete() + coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  ylab(paste0("Ratio of genes in GO category (",
              species_scientific(species), ")")) +
  xlab("") + customTheme() + theme(legend.position = "right",
                                   legend.direction = "vertical",
                                   plot.margin = unit(c(10, 5, 5, 5), "mm"))
go_plot

## -----------------------------------------------------------------------------
species <- "Mus_musculus"
expr <- combine_expr[[species]]
asbf_coef <- sbf$u[[species]]
expr[["coef"]] <- asbf_coef[, sel_dim, drop = T]
expr1 <- expr[, c(paste0(getSpeciesShortName(species), "_", sel_tissue), "Tau",
                  "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")
if (sel_sign == "neg") {
  cat("\n selecting negative loadings")
  expr1_selsign <- expr1[expr1$coef < 0, ]
} else {
  cat("\n selecting positive loadings")
  expr1_selsign <- expr1[expr1$coef >= 0, ]
}
expr1_selsign$score <- expr1_selsign$Tau * abs(expr1_selsign$coef)
expr1_selsign$rank <- rank(-1 * expr1_selsign$score)
expr1_selsign <- expr1_selsign[order(expr1_selsign$rank), ]
genes_fg <- row.names(expr1_selsign[expr1_selsign$rank <= top_genes, ])
if (species == "Homo_sapiens") {
  cat("\nUsing orthologs (human IDs) as background")
  genes_bg <- orthologs$hsapiens
}
if (species == "Mus_musculus") {
  cat("\nUsing orthologs (mouse IDs) as background")
  genes_bg <- orthologs$mmusculus
}
genes_bg <- genes_bg[!genes_bg %in% genes_fg]
genome <- "mm10"
total_genes <- unique(c(genes_fg, genes_bg))
up_genes <- as.integer(total_genes %in%  genes_fg)
names(up_genes) <- total_genes

## ---- echo = FALSE, warning = FALSE-------------------------------------------
# load("mm10 length data")
load("~/Dropbox/0.Analysis/0.paper/GOKeggFiles/mm10_length_EnsemblV94.RData")
lengthData.up <- lengthData[names(up_genes)]
# load("mm10 EnsembleID to GO data")
load("~/Dropbox/0.Analysis/0.paper/GOKeggFiles/mm10_geneID2GO.EnsemblID2GO.EnsembleV94.Robj")
pwf <- nullp(up_genes, bias.data = lengthData.up, plot.fit = F)
go <- goseq(pwf, "mm10", "ensGene", gene2cat = mouse_geneID2GO, test.cats = c("GO:BP"))
go.sub <- go[go$ontology == "BP", ]
go.sub$padj <- p.adjust(go.sub$over_represented_pvalue, method = "BH")
go.sub[["ratio"]] <- round(go.sub[["numDEInCat"]] / go.sub[["numInCat"]], 4)
go.sub <- go.sub[with(go.sub, order(padj, decreasing = c(FALSE))), ]
go.sub$over_represented_pvalue <- NULL
go.sub$under_represented_pvalue <- NULL

## -----------------------------------------------------------------------------
head(go.sub)

## -----------------------------------------------------------------------------
# GO enrichment plot for mouse
go_out <- head(go.sub, n = 8)
go_out$padj <- as.numeric(go_out$padj)
go_out$term <- factor(go_out$term, levels = go_out$term)
breaks <- round(c(0, 1 / 4, 2 / 4, 3 / 4, 1) * max(go_out[["ratio"]]), 2)
go_plot <- ggplot(go_out, aes(x = term, y = ratio, fill = padj)) + geom_col() +
  scale_y_continuous(expand = c(0, 0), breaks = breaks,
                     limits = c(0, max(go_out[["ratio"]] + 0.05))) +
  scale_x_discrete() + coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  ylab(paste0("Ratio of genes in GO category (",
              species_scientific(species), ")")) +
  xlab("") + customTheme() + theme(legend.position = "right",
                                   legend.direction = "vertical",
                                   plot.margin = unit(c(10, 5, 5, 5), "mm"))
go_plot
