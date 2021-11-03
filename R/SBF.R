#' Compute Shared Basis Factorization (SBF) and Approximate Shared Basis
#' Factorization (A-SBF)
#'
#' Function to compute Shared Basis Factorization (SBF) and
#' Approximate Shared Basis Factorization (A-SBF)
#' @param matrix_list A list containing Di matrices for joint matrix
#' factorization. Column names of each Di matrix may or may not have information
#' about tissue or cell type.
#' @param check_col_matching if the column names have information about tissue/
#' cell type and one-to-one correspondence of tissue types across species has to
#' be checked, set this parameter to be TRUE. Default FALSE.
#' @param col_sep separator in column names to separate different fields.
#' Example for column names 'hsapiens_brain', 'hsapiens_heart' etc.,
#' the separator is underscore.
#' Set it to NULL if column matching across species has to be
#' performed and there is no separator in the column names.
#' Only checked if check_col_matching = TRUE. Default underscore.
#' @param col_index If a separator separates information in column names,
#' the col_index is the index in the column name corresponding to tissue or
#' cell type. E.g. for column name 'hsapiens_brain', col_index is 2.
#' Only checked if check_col_matching = TRUE. Default NULL.
#' @param weighted If TRUE each Di^TDi is scaled using inverse variance weights
#' Default FALSE.
#' @param approximate TRUE will compute A-SBF. Default FALSE.
#' @param transform_matrix if TRUE, then Di will be transformed to compute
#' correlation matrix, and V is computed based on this instead of
#' Di^TDi. An unbiased estimate of covariance (denominator n-1) is
#' used for the computing correlation. Default FALSE.
#' @param verbose if TRUE print verbose lines. Default FALSE.
#'
#' @return a list containing u, delta, v, m, lambda (eigenvalues of m), and
#' other outputs of SBF/A-SBF factorization.
#' @export
#'
#' @examples
#' # create test dataset
#' set.seed(1231)
#' mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)
#'
#' # SBF call. Estimate V using the sum of Di^TDi
#' sbf <- SBF(matrix_list = mymat)
#'
#' # SBF call. Estimate V using inverse-variance weighted Di^TDi
#' sbf <- SBF(matrix_list = mymat, weighted = TRUE)
#' # calculate decomposition error
#' decomperror <- calcDecompError(mymat, sbf$u, sbf$delta, sbf$v)
#'
#' # SBF call using correlation matrix
#' sbf_cor <- SBF(matrix_list = mymat, transform_matrix = TRUE)
#' decomperror <- calcDecompError(mymat, sbf_cor$u, sbf_cor$delta, sbf_cor$v)
#'
#' # SBF call for gene expression dataset using correlation matrix
#' avg_counts <- SBF::TissueExprSpecies
#' sbf_cor <- SBF(matrix_list = avg_counts, transform_matrix = TRUE)
#'
#' # A-SBF call for gene expression dataset using correlation matrix
#' avg_counts <- SBF::TissueExprSpecies
#' asbf_cor <- SBF(matrix_list = avg_counts, approximate = TRUE,
#'                 transform_matrix = TRUE)
#' # calculate decomposition error
#' decomperror <- calcDecompError(avg_counts, asbf_cor$u_ortho, asbf_cor$delta,
#'                                 asbf_cor$v)
SBF <- function(matrix_list = NULL, check_col_matching = FALSE, col_sep = "_",
                col_index = NULL, weighted = FALSE,
                approximate = FALSE, transform_matrix = FALSE,
                verbose = FALSE) {
    if (length(matrix_list) >= 2 && !is.null(matrix_list)) {
        if (check_col_matching) {
            if (is.null(col_index) || !is.numeric(col_index))
                stop(paste0("\nInvalid index to match columns. Exiting!"))
            col_selected <- as.data.frame(sapply(matrix_list, function(x) {
            sapply(strsplit(as.character(colnames(x)), col_sep), function(y) {
                    y[col_index] })
            }))
            if (!all(apply(col_selected, 1, function(x) all(x == x[1]))))
                stop(paste0("\nNames are not matching.Exiting!"))
        }
        if (!all(sapply(matrix_list, ncol) == ncol(matrix_list[[names(matrix_list)[1]]])))
            stop(paste0("\nAll matrices should have same number of columns"))

        # Compute Di^TDi or Ri^TRi----------------------------------------------
        matrix_names <- names(matrix_list)
        mat_list_trans <- list()
        mat_list_trans_sum <- matrix(0L, nrow =
                                     ncol(matrix_list[[matrix_names[1]]]),
                                    ncol = ncol(matrix_list[[matrix_names[1]]]))
        if (weighted == TRUE && transform_matrix == TRUE) {
            cat(paste("\nWith inter-sample correlation transfromation,",
                "no additional scaling is required.\n", sep = " "))
        }
        if (transform_matrix == TRUE) {
            if (verbose)
                cat("\nV is computed using inter-sample correlation\n")
            for (mat in matrix_names) {
              mat_list_trans[[mat]] <- stats::cor(matrix_list[[mat]],
                                           method = "pearson")
              mat_list_trans_sum <- mat_list_trans_sum + mat_list_trans[[mat]]
            }
            mat_list_trans_sum <- mat_list_trans_sum / length(matrix_names)
        } else if (weighted) {
            if (verbose)
                cat("\nInverse variance weighting applied\n")
            tot_var_sum <- 0
            for (mat in matrix_names) {
                w_i <- sum(diag(stats::cov(matrix_list[[mat]])))
                mat_list_trans[[mat]] <- (t(as.matrix(matrix_list[[mat]])) %*%
                                            as.matrix(matrix_list[[mat]])) / w_i
                mat_list_trans_sum <- mat_list_trans_sum + mat_list_trans[[mat]]
                tot_var_sum <- tot_var_sum + w_i^-1
            }
            mat_list_trans_sum <- mat_list_trans_sum / tot_var_sum
        } else {
            # no weights and no transformation applied
            for (mat in matrix_names) {
              mat_list_trans[[mat]] <- t(as.matrix(matrix_list[[mat]])) %*%
                  as.matrix(matrix_list[[mat]])
              mat_list_trans_sum <- mat_list_trans_sum + mat_list_trans[[mat]]
            }
            mat_list_trans_sum <- mat_list_trans_sum / length(matrix_names)
        }
        colnames(mat_list_trans_sum) <- NULL
        rownames(mat_list_trans_sum) <- NULL

        # Eigen decomposition of S to find V, compute U, delta----------------
        ev <- eigen(mat_list_trans_sum)
        lambda <- ev$values
        if (!all(duplicated(lambda) == FALSE)) {
            cat("\n[Warning] Non-defective M matrix\n")
        }
        V <- ev$vectors
        B <- delta <- U <- U_ortho <- list()
        for (mat in matrix_names) {
            B[[mat]] <- as.matrix(matrix_list[[mat]]) %*% V
            delta[[mat]] <- sqrt(colSums(B[[mat]]^2))
            U[[mat]] <- sweep(B[[mat]], 2, delta[[mat]], FUN = "/")
            if (verbose) {
                cat("\nMatrix:", mat)
                cat("\nDim of matrix:", dim(matrix_list[[mat]]))
                cat("\nDim of V:", dim(V))
                cat("\nDim of B:", dim(B[[mat]]))
                cat("\nDelta vec:", delta[[mat]])
                if (length(delta[[mat]]) == 1) {
                  cat("\nDim of delta:",
                      dim(as.matrix(diag(as.matrix(delta[[mat]])))))
                } else {
                  cat("\nDim of delta:", dim(as.matrix(diag(delta[[mat]]))))
                }
                cat("\n")
            }
            if (approximate) {
                if (length(delta[[mat]]) == 1) {
                  phi <- as.matrix(diag(as.matrix(delta[[mat]])))
                } else {
                  phi <- as.matrix(diag(delta[[mat]]))
                }
                dvsig_svd <- svd(as.matrix(matrix_list[[mat]]) %*% V %*% phi)
                U_ortho[[mat]] <- dvsig_svd$u %*% t(dvsig_svd$v)
                row.names(U_ortho[[mat]]) <- row.names(U[[mat]])
            }
        }
        if (approximate) {
            if (verbose)
                cat("\nA-SBF is computed\n")
            initial_error <- calcDecompError(matrix_list, U_ortho, delta, V)
            out <- list(v = V, lambda = lambda,
                        u = U, u_ortho = U_ortho,
                        delta = delta, m = mat_list_trans_sum,
                        error = initial_error)
        } else {
            out <- list(v = V, lambda = lambda, u = U, delta = delta,
                        m = mat_list_trans_sum)
        }
        return(out)
    } else {
        stop("\nInvalid list. List of matrices expected. Exiting!")
    }
}
