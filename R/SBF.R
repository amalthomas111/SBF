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
#' be checked, set this parameter to be TRUE. Otherwise, set it to FALSE.
#' @param col_sep separator in column names to separate different fields.
#' E.g. for column names 'hsapiens_brain', 'hsapiens_heart' etc., the separator
#' is underscore. Set it to NULL if column matching across species has to be
#' performed and there is no separator in the column names.
#' Only checked if check_col_matching = TRUE. Default underscore.
#' @param col_index If a separator separates information in column names,
#' the col_index is the index in the column name corresponding to tissue or
#' cell type. E.g. for column name 'hsapiens_brain', col_index is 2.
#' Only checked if check_col_matching = TRUE. Default NULL.
#' @param approximate TRUE will compute A-SBF. Default TRUE.
#' @param transform_matrix if TRUE, then Di will be transformed to compute
#' correlation matrix, and V is computed based on this instead of
#' Di^TDi. An unbiased estimate of covariance (denominator n-1) is
#' used for the computing correlation. Default FALSE.
#' @param verbose if TRUE print verbose lines. Default FALSE.
#'
#' @return a list containing u, v, lambda, m, and other outputs
#' of SBF/A-SBF factorization.
#' @export
#'
#' @examples
#' # SBF call
#' avg_counts <- SBF::TissueExprSpecies
#' sbf <- SBF(matrix_list = avg_counts, col_index = 2, approximate = FALSE,
#'            transform_matrix = FALSE)
#'
#' # SBF call using correlation matrix
#' avg_counts <- SBF::TissueExprSpecies
#' sbf_cor <- SBF(matrix_list = avg_counts, col_index = 2, approximate = FALSE,
#'                transform_matrix = TRUE)
#'
#' # A-SBF call
#' avg_counts <- SBF::TissueExprSpecies
#' asbf <- SBF(matrix_list = avg_counts, col_index = 2, approximate = TRUE,
#'             transform_matrix = FALSE)
#' # calculate decomposition error
#' decomperror <- calcDecompError(avg_counts, asbf$delta, asbf$u_ortho, asbf$v)
#'
#' # A-SBF call using correlation matrix
#' avg_counts <- SBF::TissueExprSpecies
#' asbf_cor <- SBF(matrix_list = avg_counts, col_index = 2, approximate = TRUE,
#'                 transform_matrix = TRUE)
#' # calculate decomposition error
#' decomperror <- calcDecompError(avg_counts, asbf_cor$delta, asbf_cor$u_ortho,
#'                                 asbf_cor$v)
SBF <- function(matrix_list = NULL, check_col_matching = TRUE, col_sep = "_",
                col_index = NULL, approximate = TRUE, transform_matrix = FALSE,
                verbose = FALSE) {
    if (length(matrix_list) >= 2 & !is.null(matrix_list)) {
        if (check_col_matching) {
            if (is.null(col_index) | !is.numeric(col_index))
                stop(paste0("\nInvalid index to match columns. Exiting!"))
            col_selected <- as.data.frame(sapply(matrix_list, function(x)
              data.table::tstrsplit(colnames(x), col_sep)[[col_index]]))
            if (!all(apply(col_selected, 1, function(x) all(x == x[1]))))
                stop(paste0("\nNames are not matching.Exiting!"))
        }
        if (!all(sapply(matrix_list, ncol) == ncol(matrix_list[[names(matrix_list)[1]]])))
            stop(paste0("\nAll matrices should have same number of columns"))
        if (approximate)
            cat("\nA-SBF is computed\n")
        if (transform_matrix)
            cat("\nV is computed using inter-sample correlation\n")
        # Compute Di^TDi or Ri^TRi----------------------------------------------
        matrix_names <- names(matrix_list)
        mat_list_trans <- list()
        mat_list_trans_sum <- matrix(0L, nrow =
                                     ncol(matrix_list[[matrix_names[1]]]),
                                    ncol = ncol(matrix_list[[matrix_names[1]]]))
        if (transform_matrix) {
            for (mat in matrix_names) {
              mat_list_trans[[mat]] <- stats::cor(matrix_list[[mat]],
                                           method = "pearson")
              mat_list_trans_sum <- mat_list_trans_sum + mat_list_trans[[mat]]
            }
        } else {
            for (mat in matrix_names) {
              mat_list_trans[[mat]] <- t(as.matrix(matrix_list[[mat]])) %*%
                  as.matrix(matrix_list[[mat]])
              mat_list_trans_sum <- mat_list_trans_sum + mat_list_trans[[mat]]
            }
        }
        mat_list_trans_sum <- mat_list_trans_sum / length(matrix_names)
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
            U[[mat]] <- sweep(B[[mat]], 2, delta[[mat]], "/")
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
            initial_error <- calcDecompError(matrix_list, delta, U_ortho, V)
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
