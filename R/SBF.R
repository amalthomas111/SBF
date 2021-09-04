#' Function to compute Shared Basis Factorization (SBF) and
#' Approximate Shared Basis Factorization (A-SBF)
#' @param avg_counts A list containing Di matrices for joint matrix
#' factorization. Column names of each Di matrix may or may not have information
#' about tissue or cell type.
#' @param check_col_matching if the column names have information about tissue/
#' cell type and one-to-one correspondence of tissue types across species has to
#' be checked, set this parameter to be TRUE. Otherwise, set it to FALSE.
#' @param colNameSep separator in column names to separate different fields.
#' E.g. for column names 'hsapiens_brain', 'hsapiens_heart' etc., the separator
#' is underscore. Set it to NULL if column matching across species has to be
#' performed and there is no separator in the column names.
#' Only checked if check_col_matching = TRUE. Default underscore.
#' @param colIndex If a separator separates information in column names,
#' the colIndex is the index in the column name corresponding to tissue or
#' cell type. E.g. for column name 'hsapiens_brain', colIndex is 2.
#' Only checked if check_col_matching = TRUE. Default NULL
#' @param approximate TRUE will compute A-SBF. Default TRUE
#' @param transformData if TRUE, then Di will be transformed to compute
#' correlation matrix, and V is computed based on this instead of
#' Di^TDi. An unbiased estimate of covariance (denominator n-1) is
#' used for the computing correlation. Default FALSE
#' @param verbose if TRUE print verbose lines. Default FALSE
#'
#' @return
#' @export
#'
#' @examples
#' # SBF call
#' avg_counts <- SBF::avg_counts
#' sbf <- SBF(avg_counts = avg_counts,colIndex = 2, approximate = FALSE,
#'          transformData = FALSE)
#'
#' # SBF call using correlation matrix
#' avg_counts <- SBF::avg_counts
#' sbf.cor <- SBF(avg_counts = avg_counts,colIndex = 2, approximate = FALSE,
#'               transformData = TRUE)
#'
#' # A-SBF call
#' avg_counts <- SBF::avg_counts
#' asbf <- SBF(avg_counts = avg_counts,colIndex = 2, approximate = TRUE,
#'            transformData = FALSE)
#'
#' # A-SBF call using correlation matrix
#' avg_counts <- SBF::avg_counts
#' asbf.cor <- SBF(avg_counts = avg_counts,colIndex = 2, approximate = TRUE,
#'                transformData = TRUE)
SBF <- function(avg_counts = NULL, check_col_matching = TRUE, colNameSep = "_",
                colIndex = NULL, approximate = TRUE, transformData = FALSE,
                verbose = FALSE) {
    if (length(avg_counts) >= 2 & !is.null(avg_counts)) {
        if (is.null(colIndex) | !is.numeric(colIndex))
            stop(paste0("\nInvalid Index to match Columns(Tissues/Celltypes).
                        Exiting!"))
        if (check_col_matching) {
            c_tissues <- as.data.frame(sapply(avg_counts, function(x)
              data.table::tstrsplit(colnames(x), colNameSep)[[colIndex]]))
            if (!all(apply(c_tissues, 1, function(x) all(x == x[1]))))
                stop(paste0("\nTissue names are not matching.Exiting!"))
        }
        if (!all(sapply(avg_counts, ncol) == ncol(avg_counts[[names(avg_counts)[1]]])))
            stop(paste0("\nAll data.frame should have same number of columns"))
        if (approximate)
            cat("\nA-SBF is computed")
        if (transformData)
            cat("\nV is computed using correlation matrix")
        species.names <- names(avg_counts)
        A <- list()
        S <- matrix(0L, nrow = ncol(avg_counts[[species.names[1]]]),
                    ncol = ncol(avg_counts[[species.names[1]]]))
        if (transformData) {
            for (species in species.names) {
                A[[species]] <- stats::cor(avg_counts[[species]], method = "pearson")
                S <- S + A[[species]]
            }
        } else {
            for (species in species.names) {
                A[[species]] <- t(as.matrix(avg_counts[[species]])) %*%
                  as.matrix(avg_counts[[species]])
            }
        }
        S <- S / length(species.names)
        # eigen decomposition of S
        ev <- eigen(S)
        lambda <- ev$values
        V <- ev$vectors
        B <- sigma <- sigma.ortho <- U <- U.ortho <- list()
        for (species in species.names) {
            B[[species]] <- as.matrix(avg_counts[[species]]) %*% V
            sigma[[species]] <- sqrt(colSums(B[[species]]^2))
            U[[species]] <- sweep(B[[species]], 2, sigma[[species]], "/")
            if (verbose) {
                cat("\nspecies:", species)
                cat("\nDim of avg counts:", dim(avg_counts[[species]]))
                cat("\nDim of V:", dim(V))
                cat("\nDim of B:", dim(B[[species]]))
                cat("\nSigma vec:", sigma[[species]])
                if (length(sigma[[species]]) == 1) {
                  cat("\nDim of sigma:",
                      dim(as.matrix(diag(as.matrix(sigma[[species]])))))
                } else {
                  cat("\nDim of sigma:", dim(as.matrix(diag(sigma[[species]]))))
                }
            }
            if (approximate) {
                if (length(sigma[[species]]) == 1) {
                  phi <- as.matrix(diag(as.matrix(sigma[[species]])))
                } else {
                  phi <- as.matrix(diag(sigma[[species]]))
                }
                dvsig.svd <- svd(as.matrix(avg_counts[[species]]) %*% V %*% phi)
                U.ortho[[species]] <- dvsig.svd$u %*% t(dvsig.svd$v)
                row.names(U.ortho[[species]]) <- row.names(U[[species]])
                sigma.ortho[[species]] <- diag(t(U.ortho[[species]]) %*%
                                                 as.matrix(avg_counts[[species]]) %*% solve(t(V)))
            }
        }
        if (approximate) {
            initial.error <- calcDecompError(avg_counts, sigma, U.ortho, V)
            initial.error.ortho <- calcDecompError(avg_counts, sigma.ortho,
                                                   U.ortho, V)
            if (verbose) {
                cat("\nInitial.Error =", round(initial.error$error, 2),
                    "after ortho sigma=", round(initial.error.ortho$error, 2))
            }
        }
        if (approximate) {
            output.list <- list(v = V, lambda = lambda, u = U,
                                u_ortho = U.ortho, sigma = sigma,
                                sigma.ortho = sigma.ortho, s = S,
                                error = initial.error$error,
                                error.ortho = initial.error.ortho$error)
        } else {
            output.list <- list(v = V, lambda = lambda, u = U,
                                sigma = sigma, s = S)
        }
        return(output.list)
    } else {
        stop("\nInvalid counts file. Counts file should be list
         with matching column names.Exiting!")
    }
}
