#' Function to create short name for species
#'
#' @param species_list a vector of species name
#'
#' @return short version of species name
#' @export
#'
#' @examples
#' species <- c("Homo_sapiens", "Mus_musculus")
#' species_short <- sapply(species, getSpeciesShortName)
#' species <- "Homo_sapiens"
#' getSpeciesShortName(species)
getSpeciesShortName <- function(species_list) {
  shortname <- c()
  for (sp in species_list) {
    if (!grepl("_", sp)) {
      cat("\nspecies:", sp, "\n")
      stop("Invalid species full name: No '_' found!")
    }
    shortname <- c(shortname, paste0(tolower(substr(base::strsplit(sp, "_")[[1]][1],
                                                    1, 1)),
                                     base::strsplit(sp, "_")[[1]][2]))
  }
  return(shortname)
}
#' Function to return scientific name of species
#'
#' @param x species name
#'
#' @return species name in scientific format
#' @export
#'
#' @examples
#' species <- "Homo_sapiens"
#' getScientificName(species)
getScientificName <- function(x) {
  if (!grepl("_", x)) {
    cat("\nspecies:", x, "\n")
    stop("Invalid species full name: No '_' found!")
  }
  parts <- base::strsplit(x, "_")[[1]]
  return(paste0(substr(x, 1, 1), ".", parts[2]))
}
#' Calculate TPM normalization for a vector
#'
#' @param x a column of gene expression profiles
#' @param gene_len gene lengths (kb)
#'
#' @return
#' @export
#'
#' @examples
countToTpm <- function(x, gene_len) {
  rate <- log(x) - log(gene_len)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
#' TPM normalization for a gene expression matrix
#'
#' @param rawCounts Raw gene expression counts matrix
#' @param gene_len Matching gene lengths (kb)
#'
#' @return TPM normalized matrix
#' @export
#'
#' @examples
normalizeTPM <- function(rawCounts, gene_len = NULL) {
  common_genes <- intersect(rownames(rawCounts), names(gene_len))
  rawCounts <- rawCounts[common_genes, , drop = FALSE]
  gene_len <- gene_len[common_genes]
  glen <- gene_len[match(toupper(rownames(rawCounts)),
                            toupper(names(gene_len)))]
  if (nrow(rawCounts) == length(gene_len)) {
    cat("\nTPM counts returned")
    df <- as.data.frame(apply(rawCounts, 2, function(col) {
       countToTpm(x = col, gene_len = glen)}))
    df <- round(df, 4)
  } else {
    stop("\tUnequal number of genes for TPM calculation")
  }
}

#' Function to calculate mean expression values for each tissue/group
#'
#' @param counts gene expression matrix
#' @param metadata data frame with species tissue key column
#' @param ndecimal number of decimals to round
#'
#' @return
#' @export
#'
#' @examples
calcAvgCounts <- function(counts, metadata, ndecimal = 4) {
  counts_avg <- list()
  if ("key" %in% colnames(metadata)) {
    for (species_organ in sort(unique(metadata$key))) {
      nlibs <- length(colnames(counts)[grepl(species_organ,
                                                  colnames(counts))])
      if (nlibs > 1)
        counts_avg[[species_organ]] <- round(rowMeans(counts[, colnames(counts)[
          grepl(species_organ, colnames(counts))]],
          dim = 1), ndecimal)
      else if (nlibs == 1)
        counts_avg[[species_organ]] <- round(counts[, colnames(counts)[grepl(
          species_organ, colnames(counts))], drop = TRUE], ndecimal)
    }
    avg_counts <- as.data.frame(do.call(cbind, counts_avg))
    #row.names(avg_counts) <- row.names(counts)
    return(avg_counts)
  } else {
    stop("'key' column not found in the metadata!Exiting!")
  }
}

#' Compute information represented by each dimensions
#'
#' @param l SBF list object
#'
#' @return list with percentage of information represented by each dimension
#' according to their Delta_i values
#' @export
#'
#' @examples
#' #' @examples
#' # create test dataset
#' set.seed(1231)
#' mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)
#'
#' # SBF call. Estimate V using the sum of Di^TDi
#' sbf <- SBF(matrix_list = mymat)
#' calcPercentInfo(sbf)
calcPercentInfo <- function(l) {
  if (is.list(l) && ("delta" %in% names(l))) {
  lapply(l$delta, function(x) {
    round(x^2 / sum(x^2) * 100, 2)
  })
  } else {
    stop("Not an output list of SBF call")
  }
}
#' Invert delta_i entries and return a matrix
#'
#' @param l SBF list object
#'
#' @return a matrix with inverted delta_i values
#' @export
#'
#' @examples
invertDelta <- function(l) {
  if (!is.list(l) || !("delta" %in% names(l)))
    stop("Not an output list of SBF call")
  d_inv <- list()
  for (sp in names(l$delta)) {
    if (length(l$delta[[sp]]) == 1) {
      d_inv[[sp]] <- as.matrix(diag(as.matrix(1 / l$delta[[sp]])))
    } else {
      d_inv[[sp]] <- as.matrix(diag(1 / l$delta[[sp]]))
    }
  }
  return(d_inv)
}

#' Function to project counts to the common expression space
#' This function projects individual profiles or average counts to the
#' common space by computing D_i^T U_i Delta^{-1}
#'
#' @param counts mean expression list or counts list
#' @param obj SBF output list
#' @param combine If multiple libraries exist, the function
#' will row bind the projected counts if set TRUE
#' @return projected counts data.frame if combine is set TRUE else
#' list with projected counts
#' @export
#'
#' @examples
projectCounts <- function(counts, obj, combine = TRUE) {
  if (!is.list(counts) || !is.list(obj))
    stop("counts and sbf object should be lists")
  if (!all(c("delta", "u") %in% names(obj)))
    stop("Invalid sbf list")
  if (!all(names(counts) == names(obj$u)) ||
      !all(names(counts) == names(obj$delta)))
    stop("Names of counts not match sbf's u or delta")
  if (nrow(counts[[names(counts)[1]]]) != nrow(obj$u[[names(obj$u)[1]]]))
    stop("# of genes not matching counts and osbf$u")
  counts_projected <- list()
  counts_rnames <- c()
  d_inv <- invertDelta(obj)
  for (sp in names(counts)) {
    counts_projected[[sp]] <- as.matrix(t(counts[[sp]])) %*%
      as.matrix(obj$u[[sp]]) %*% d_inv[[sp]]
    counts_rnames <- c(counts_rnames, row.names(counts_projected[[sp]]))
  }
  if (combine) {
    df_proj <- as.data.frame(do.call(rbind, counts_projected))
    rownames(df_proj) <- counts_rnames
    colnames(df_proj) <- paste0("Dim", seq_len(ncol(df_proj)))
    return(df_proj)
  } else {
    return(counts_projected)
  }
}
