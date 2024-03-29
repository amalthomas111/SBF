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
cat("\nOSBF with optimizing V = FALSE started\n")
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n")
t1 <- proc.time()
# decrease tol to minimize error
osbf_noVupdate <- SBF(avg_counts, transform_matrix = TRUE, orthogonal = TRUE,
                     optimizeV = FALSE, tol = 1e-3)#, tol = 1e-10)
t2 <- proc.time()
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n")
cat("OSBF with optimizing V = FALSE finished\n")
cat("Time taken:\n")
t2 - t1

## -----------------------------------------------------------------------------
cat("\nOSBF with optimizing V=TRUE started\n")
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n")
t1 <- proc.time()
# decrease tol to minimize error
osbf <- SBF(avg_counts, transform_matrix = TRUE, orthogonal = TRUE,
                     optimizeV = TRUE, tol = 1e-3)#, tol = 1e-10)
t2 <- proc.time()
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n")
cat("OSBF with optimizing V=TRUE finished\n")
cat("Time taken:\n")
t2 - t1

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
cat("\nPercentage for each delta [No V update]:")
percentInfo_noVupdate <- calcPercentInfo(osbf_noVupdate)
for (i in names(osbf_noVupdate$delta)) {
  cat("\n", sprintf("%-25s:", i), sprintf("%8.2f", percentInfo_noVupdate[[i]]))
}
percentInfo <- calcPercentInfo(osbf)
for (i in names(osbf$delta)) {
  cat("\n", sprintf("%-25s:", i), sprintf("%8.2f", percentInfo[[i]]))
}

## -----------------------------------------------------------------------------
# project profles using no V update estimates
# we can project both mean expression profiles as well as individual expression
# profiles
df_proj_avg_noVupdate <- projectCounts(avg_counts, osbf_noVupdate)
meta <- data.table::tstrsplit(row.names(df_proj_avg_noVupdate), "_")
df_proj_avg_noVupdate$tissue <- factor(meta[[2]])
df_proj_avg_noVupdate$species <- factor(meta[[1]])
df_proj_avg_noVupdate <- df_proj_avg_noVupdate %>%
  mutate(species = factor(species, levels = species_short))


df_proj_noVupdate <- projectCounts(counts_list_sub, osbf_noVupdate)
meta1 <- data.table::tstrsplit(row.names(df_proj_noVupdate), "_")
df_proj_noVupdate$tissue <- factor(meta1[[3]])
df_proj_noVupdate$species <- factor(meta1[[2]])
df_proj_noVupdate <- df_proj_noVupdate %>% mutate(species = factor(species,
                            levels = species_short))

## -----------------------------------------------------------------------------
# project using V update estimates
df_proj_avg <- projectCounts(avg_counts, osbf)
meta <- data.table::tstrsplit(row.names(df_proj_avg), "_")
df_proj_avg$tissue <- factor(meta[[2]])
df_proj_avg$species <- factor(meta[[1]])
df_proj_avg <- df_proj_avg %>% mutate(species = factor(species,
                            levels = species_short))


df_proj <- projectCounts(counts_list_sub, osbf)
meta1 <- data.table::tstrsplit(row.names(df_proj), "_")
df_proj$tissue <- factor(meta1[[3]])
df_proj$species <- factor(meta1[[2]])
df_proj <- df_proj %>% mutate(species = factor(species,
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
data_noVupdate <- df_proj_noVupdate
data_noVupdate$tissue <- NULL
data_noVupdate$species <- NULL
data_noVupdate <- as.matrix(data_noVupdate)
data_noVupdate_dist <- as.matrix(dist(data_noVupdate, method = "euclidean"))
meta <- data.table::tstrsplit(colnames(data_noVupdate_dist), "_")
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

myheatmap1 <- ComplexHeatmap::Heatmap(as.matrix(data_noVupdate_dist),
                                      cluster_rows = TRUE,
                                      clustering_method_rows = "centroid",
                                      cluster_columns = TRUE,
                                      clustering_method_columns = "centroid",
                                      top_annotation = ht, col = morecolors(50),
                                      show_row_names = FALSE,
                                      show_column_names = FALSE,
                                      name = "distance")
myheatmap1

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
# 2D plot for Dim1 and Dim2 [No V update]
sel_colors <- RColorBrewer::brewer.pal(8, "Dark2")[seq_len(length(unique(df_proj_noVupdate$tissue)))]
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
i <- 1
j <- 2
ggplot2::ggplot(df_proj, aes(x = df_proj[, i], y = df_proj[, j], col = tissue,
                             shape = species, fill = tissue)) +
  xlab(paste("OSBF Dim", i)) + ylab(paste("OSBF Dim", j)) +
  geom_point(size = 1.5) + scale_color_manual(values = sel_colors) +
  scale_shape_manual(values = c(21:25, 3:7)) +
  scale_fill_manual(values = sel_colors) +
  customTheme() +
  theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# 2D plot for Dim2 and Dim3 [No V update]
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
sel_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "mediumturquoise",
                "#E6AB02", "darkmagenta", "#666666", "black", "darkolivegreen2")
outputname <- "human_mouse_blood"
finished <- c()
for (i in 1:(ncol(df_proj) - 2)) {
  for (j in 1:(ncol(df_proj) - 2)) {
    if (i == j) next
    if (j %in% finished) next
    ggplot(df_proj, aes(x = df_proj[, i], y = df_proj[, j],
                            col = tissue, shape = species, fill = tissue)) +
      xlab(paste0("OSBF Dim ", i)) +
      ylab(paste0("OSBF Dim ", j)) +
      geom_point(size = 1.5) +
      scale_color_manual(values = sel_colors) +
      scale_shape_manual(values = c(21:25, 3:7)) +
      scale_fill_manual(values = sel_colors) +
      customTheme(base_size = 12) +
      theme(legend.title = element_blank())
    #ggsave(filename = paste0(outdir, "2Dplots/opt_2Dplot_Dim_", i, "-", j, "_",
    #                         outputname, ".pdf"), device = "pdf",
    #       width = 7, height = 7, useDingbats = FALSE)
  }
  finished <- c(finished, i)
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
sel_dim <- 2
sel_tissue <- "testis"
species <- "Homo_sapiens"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = TRUE]
expr1 <- expr[, c(paste0(getSpeciesShortName(species), "_", sel_tissue),
                   "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")
head(expr1)

## -----------------------------------------------------------------------------
# plot scatter
mid <- 0
p1 <- ggplot2::ggplot(expr1, aes(x = Tau, y = coef, col = tissue_zscore)) +
  theme_bw() +
  geom_point(size = 0.5) + xlab("Expression specificity") +
  ylab(paste0("Dim", sel_dim, " Coefficient")) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red", space = "Lab") +
  scale_y_continuous(limits = c(-1 * max(abs(expr1$coef)),
                                max(abs(expr1$coef))),
                     breaks = seq(-1 * round(max(abs(expr$coef)), 2),
                                  round(max(abs(expr$coef)), 2), by = 0.01)) +
  customTheme() +  theme(legend.position = "right",
                         legend.direction = "vertical") +
  labs(title = paste0(getScientificName(species)), color = "Z-score") +
  theme(legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(face = "italic"))
p1

## -----------------------------------------------------------------------------
sel_dim <- 2
sel_tissue <- "testis"
species <- "Sus_scrofa"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = TRUE]
expr1 <- expr[, c(paste0(getSpeciesShortName(species), "_", sel_tissue),
                   "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")

# plot scatter
mid <- 0
p2 <- ggplot2::ggplot(expr1, aes(x = Tau, y = coef, col = tissue_zscore)) +
  theme_bw() +
  geom_point(size = 0.5) + xlab("Expression specificity") +
  ylab(paste0("Dim", sel_dim, " Coefficient")) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red", space = "Lab") +
  scale_y_continuous(limits = c(-1 * max(abs(expr1$coef)),
                                max(abs(expr1$coef))),
                     breaks = seq(-1 * round(max(abs(expr$coef)), 2),
                                  round(max(abs(expr$coef)), 2), by = 0.01)) +
  customTheme() +  theme(legend.position = "right",
                         legend.direction = "vertical") +
  labs(title = paste0(getScientificName(species)), color = "Z-score") +
  theme(legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(face = "italic"))
p2

## -----------------------------------------------------------------------------
# load GO analysis files
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
file <- "allwayOrthologs_hsapiens-ptroglodytes-mmulatta-mmusculus-rnorvegicus-btaurus-sscrofa-ggallus_ens94.tsv"
orthologs <- read.table(file.path(path, "GOKeggFiles/", file), header = TRUE,
                        sep = "\t")
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
# axis positive (pos) or negative (neg)
sel_sign <- "neg"

## -----------------------------------------------------------------------------
species <- "Homo_sapiens"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = TRUE]
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
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
# load("hg38 length data")
load(paste0(path, "GOKeggFiles/hg38_length.EnsemblV94.RData"))
lengthData.up <- lengthData[names(up_genes)]
# load("hg38 EnsembleID to GO data")
load(paste0(path, "GOKeggFiles/hg38_geneID2GO.EnsemblID2GO.EnsembleV94.Robj"))
pwf <- goseq::nullp(up_genes, bias.data = lengthData.up, plot.fit = FALSE)
go <- goseq::goseq(pwf, "hg38", "ensGene", gene2cat = geneID2GO,
                   test.cats = c("GO:BP"))
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
go_plot <- ggplot2::ggplot(go_out, aes(x = term, y = ratio, fill = padj)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0), breaks = breaks,
                     limits = c(0, max(go_out[["ratio"]] + 0.05))) +
  scale_x_discrete() + coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  ylab(paste0("Ratio of genes in GO category (",
              getScientificName(species), ")")) +
  xlab("") + customTheme() + theme(legend.position = "right",
                                   legend.direction = "vertical",
                                   plot.margin = unit(c(10, 5, 5, 5), "mm"))
go_plot

## -----------------------------------------------------------------------------
species <- "Mus_musculus"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = TRUE]
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
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
# load("mm10 length data")
load(paste0(path, "/GOKeggFiles/mm10_length_EnsemblV94.RData"))
lengthData.up <- lengthData[names(up_genes)]
# load("mm10 EnsembleID to GO data")
load(paste0(path, "/GOKeggFiles/mm10_geneID2GO.EnsemblID2GO.EnsembleV94.Robj"))
pwf <- goseq::nullp(up_genes, bias.data = lengthData.up, plot.fit = FALSE)
go <- goseq::goseq(pwf, "mm10", "ensGene", gene2cat = mouse_geneID2GO,
                   test.cats = c("GO:BP"))
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
go_plot <- ggplot2::ggplot(go_out, aes(x = term, y = ratio, fill = padj)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0), breaks = breaks,
                     limits = c(0, max(go_out[["ratio"]] + 0.05))) +
  scale_x_discrete() + coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  ylab(paste0("Ratio of genes in GO category (",
              getScientificName(species), ")")) +
  xlab("") + customTheme() + theme(legend.position = "right",
                                   legend.direction = "vertical",
                                   plot.margin = unit(c(10, 5, 5, 5), "mm"))
go_plot

## -----------------------------------------------------------------------------
# set seed
s1 <- 135
s2 <- 13835
species <- c("Homo_sapiens", "Pan_troglodytes", "Macaca_mulatta",
             "Mus_musculus", "Rattus_norvegicus", "Bos_taurus",
             "Sus_scrofa", "Gallus_gallus")
species_short <- sapply(species, getSpeciesShortName)
tissues <-  c("brain", "heart", "kidney", "liver", "lung", "testis")
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
counts_list_shuff <- metadata_list_shuff <- avg_counts_shuff <- list()
for (sp in species) {
  # reading raw counts
  counts <- read.table(paste0(path, "counts/", sp, "_rawcounts.tsv"),
                       header = TRUE, sep = "\t", row.names = 1)
  info <- data.table::tstrsplit(colnames(counts), "_")
  metadata <- data.frame(project = info[[1]],
        species = info[[2]],
        tissue = info[[3]],
        gsm = info[[4]],
        name = colnames(counts),
        stringsAsFactors = FALSE)
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
  gene_length <- read.table(paste0(path, "ensembl94_annotation/", sp,
                                   "_genelength.tsv"), sep = "\t",
                            header = TRUE, row.names = 1,
                            stringsAsFactors = FALSE)
  if (!all(row.names(counts_avg) %in% row.names(gene_length))) stop("Error")
  gene_length$Length <- gene_length$Length / 1e3
  gene_length <- gene_length[row.names(counts_avg), , drop = TRUE]
  names(gene_length) <- row.names(counts_avg)
  counts_tpm <- normalizeTPM(rawCounts = counts_avg, gene_len = gene_length)
  min_tpm <- 1
  counts_tpm[counts_tpm < min_tpm] <- 1
  counts_tpm <- log2(counts_tpm)

  info <- data.table::tstrsplit(colnames(counts_tpm), "_")
  metadata <- data.frame(
        species = info[[1]],
        tissue = info[[2]],
        name = colnames(counts_tpm),
        stringsAsFactors = FALSE)
  metadata$key <- paste0(metadata$species, "_", metadata$tissue)
  avg_counts_shuff[[sp]] <- calcAvgCounts(counts_tpm, metadata)
  # write.table(avg_counts_shuff[[sp]],file=paste0(mainDir, "shuffledavgCounts/avg_logTPM_",
  #                               species,"_SHUFF_counts_seed",z,"_",z1,".tsv"),
  #           sep="\t",quote=F,col.names=NA)
  counts_list_shuff[[sp]] <- counts_tpm
  metadata_list_shuff[[sp]] <- metadata
}
sapply(counts_list_shuff, dim)
# remove zero counts
avg_counts_shuff <- lapply(avg_counts_shuff, removeZeros)
sapply(avg_counts_shuff, dim)

## -----------------------------------------------------------------------------
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n")
# decrease tol to minimize error
osbf_shuf <- SBF(avg_counts_shuff, transform_matrix = TRUE, orthogonal = TRUE,
                 tol = 1e-2)
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n")
osbf_shuf$error

## -----------------------------------------------------------------------------
Tau_null <- lapply(avg_counts_shuff, function(x) {calc_tissue_specificity(x)})
avg_counts_shuff_scaled <- lapply(avg_counts_shuff, function(x) {
  t(scale(t(x)))
  })
combine_expr_null <- list()
for (sp in names(avg_counts_shuff_scaled)) {
  x <- as.data.frame(avg_counts_shuff_scaled[[sp]])
  x[["Tau"]] <- Tau_null[[sp]]
  combine_expr_null[[sp]] <- x
}

## -----------------------------------------------------------------------------
sel_dim <- 2
sel_tissue <- "testis"
# axis positive (pos) or negative (neg)
sel_sign <- "neg"
species <- "Homo_sapiens"
species_short <- "hsapiens"
expr <- combine_expr[[species]]
osbf_coef <- osbf$u[[species]]
expr[["coef"]] <- osbf_coef[, sel_dim, drop = TRUE]
expr1 <- expr[, c(paste0(species_short, "_", sel_tissue),
                   "Tau", "coef")]
colnames(expr1) <- c("tissue_zscore", "Tau", "coef")

# null loadings for the same dimensions
expr_null <- combine_expr_null[[species]]
null_u <- osbf_shuf$u[[species]]
expr_null[["coef"]] <- null_u[, sel_dim, drop = TRUE]
expr1_null <- expr_null[, c(paste0(species_short, "_", sel_tissue), "Tau",
                            "coef")]
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
expr1$pvalue <- sapply(expr1$score, function(x) {
  sum(as.integer(expr1_null$score > x)) / length(expr1_null$score)
  })
head(expr1)

## -----------------------------------------------------------------------------
# cut off for the p-value
alpha <- 1e-3
summary(expr1$pvalue <= alpha)

## -----------------------------------------------------------------------------
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
gene_info <- read.table(paste0(path, "ensembl94_annotation/", species_short,
                               "_genes_completeinfo.tsv"),
                               sep = "\t", header = TRUE, quote = "\"")
gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id), ]
gene_info <- gene_info[gene_info$ensembl_gene_id %in% row.names(expr1), ]
row.names(gene_info) <- gene_info$ensembl_gene_id
gene_info <- gene_info[row.names(expr1), ]
expr1$gene_name <- gene_info$external_gene_name
expr1$biotype <- gene_info$gene_biotype
head(expr1, n = 10)

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
  ylab(paste0("Dim", sel_dim, " Coefficient")) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red", space = "Lab") +
  customTheme() +  theme(legend.position = "right",
                         legend.direction = "vertical") +
  labs(title = getScientificName(species), color = "Tau") +
  theme(legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(face = "italic"))
p1

## -----------------------------------------------------------------------------
sel_dim <- 1
combine_expr_Dim1 <- list()
for (sp in names(combine_expr)) {
  expr <- combine_expr[[sp]]
  osbf_coef <- osbf$u[[sp]]
  expr[["coef"]] <- osbf_coef[, sel_dim, drop = TRUE]
  expr_null <- combine_expr_null[[sp]]
  null_u <- osbf_shuf$u[[sp]]
  expr_null[["coef"]] <- null_u[, sel_dim, drop = TRUE]
  expr$score <- abs(expr$coef)
  expr$rank <- rank(-1 * expr$score)
  expr_null$score <- abs(expr_null$coef)
  expr[[paste0("Dim", sel_dim, "_pval")]] <- sapply(expr$score, function(x) {
    sum(as.integer(expr_null$score > x)) / length(expr_null$score)
    })
  combine_expr_Dim1[[sp]] <- expr
}

## -----------------------------------------------------------------------------
# cut off for the p-value
alpha <- 1e-3
sapply(combine_expr_Dim1, function(x) {summary(x$Dim1_pval <= alpha)})

## -----------------------------------------------------------------------------
# set the path to the working directory. Change this accordingly
path <- "~/Dropbox/0.Analysis/0.paper/"
file <- "allwayOrthologs_hsapiens-ptroglodytes-mmulatta-mmusculus-rnorvegicus-btaurus-sscrofa-ggallus_ens94.tsv"
orthologs <- read.table(file.path(path, "GOKeggFiles/", file), header = TRUE,
                        sep = "\t")

Dim1_genes <- list()
for (sp in names(combine_expr_Dim1)) {
  x <- combine_expr_Dim1[[sp]][combine_expr_Dim1[[sp]]$Dim1_pval <= alpha, ]
  x <- x[order(x$rank), c("rank", "score", "Tau", "Dim1_pval"), drop = FALSE]
  # get gene details
  gene_info <- read.table(paste0(path, "ensembl94_annotation/",
                                 getSpeciesShortName(sp),
                                 "_genes_completeinfo.tsv"),
                                 sep = "\t", header = TRUE, quote = "\"")
  gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id), ]
  gene_info <- gene_info[gene_info$ensembl_gene_id %in% row.names(x), ]
  row.names(gene_info) <- gene_info$ensembl_gene_id
  gene_info <- gene_info[row.names(x), ]
  gene_info <- gene_info[, c("ensembl_gene_id", "external_gene_name",
                             "gene_biotype", "chromosome_name")]
  colnames(gene_info) <- c("id", "gene_name", "biotype", "chr")
  species_allway_orthogs <- orthologs[, getSpeciesShortName(sp), drop = TRUE]
  ortholog_status <- rep(FALSE, nrow(x))
  ortholog_status[row.names(x) %in%  species_allway_orthogs] <- TRUE
  x$allwayOrtholog <- as.factor(ortholog_status)
  df <- cbind(gene_info, x)
  row.names(df) <- NULL
  Dim1_genes[[sp]] <- df
}

## -----------------------------------------------------------------------------
head(Dim1_genes[["Homo_sapiens"]])

## -----------------------------------------------------------------------------
head(Dim1_genes[["Mus_musculus"]])

## -----------------------------------------------------------------------------
cnames <- c("Dim", "species", "nTotal", "nCoding", "nMt", "nAllwayOrtholog")
summary_stats <- data.frame(matrix(NA, nrow = 0, ncol = length(cnames)))
colnames(summary_stats) <- cnames

for (sp in names(Dim1_genes)) {
  df <- Dim1_genes[[sp]]
  stats <- data.frame(matrix(0L, nrow = 1, ncol = length(cnames)))
  colnames(stats) <- cnames
  stats$Dim <- "Dim1"
  stats$species <- getSpeciesShortName(sp)
  stats$nTotal <- nrow(df)
  stats$nCoding <- round(nrow(df[df$biotype == "protein_coding",
                                 , drop = FALSE]) * 100 / nrow(df), 2)
  stats$nMt <- round(nrow(df[tolower(df$chr) == "mt",
                             , drop = FALSE]) * 100 / nrow(df), 2)
  stats$nAllwayOrtholog <- round(nrow(df[df$allwayOrtholog == TRUE,
                                , drop = FALSE]) * 100 / nrow(df), 2)
  summary_stats <- rbind(summary_stats, stats)
}
summary_stats

## -----------------------------------------------------------------------------
sessionInfo()

