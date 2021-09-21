#' Calculate decomposition error
#'
#' Function to compute Forbenius norm between two sets of matrix.
#' Computes sum of Forbenius norm for
#' \code{matrix_initial[[i]] - u[[i]]  diag(delta[[i]])  t(v)} for
#' all i.
#'
#' @param matrix_initial list with initial Di matrices
#' @param delta delta_i values computed using SBF/A-SBF function
#' @param u  U_i values computed using SBF/A-SBF function
#' @param v V computed using SBF/A-SBF function
#'
#' @return a numeric value for the factorization error
#' @export
#'
#' @examples
#' # load dataset
#' avg_counts <- SBF::TissueExprSpecies
#' # call sbf
#' sbf <- SBF(matrix_list = avg_counts, check_col_matching = TRUE,
#'            col_index = 2, approximate = FALSE,
#'            transform_matrix = FALSE, verbose = FALSE)
#'
#' # calculate decomposition error
#' decomperror <- calcDecompError(avg_counts, sbf$delta, sbf$u, sbf$v)
#'
#' # e.g. 2
#' avg_counts <- SBF::TissueExprSpecies
#' asbf_cor <- SBF(matrix_list = avg_counts, check_col_matching = TRUE,
#'                 col_index = 2, approximate = TRUE,
#'                 transform_matrix = TRUE, verbose = FALSE)
#' decomperror <- calcDecompError(avg_counts, asbf_cor$delta, asbf_cor$u_ortho,
#'                                asbf_cor$v)
calcDecompError <- function(matrix_initial, delta, u, v) {
    decomp_error <- 0
    matrix_new <- list()
    if (! all(names(matrix_initial) == names(delta)))
        stop("Different names for matrix and delta list")
    for (mat in names(matrix_initial)) {
        if (length(delta[[mat]]) == 1)
            phi <- as.matrix(diag(as.matrix(delta[[mat]])))
        else phi <- as.matrix(diag(delta[[mat]]))
        matrix_new[[mat]] <- as.matrix(u[[mat]]) %*% phi %*%
                                                as.matrix(t(v))
        if (!identical(dim(matrix_new[[mat]]), dim(matrix_initial[[mat]])))
            stop("Dimension of the matrices not matching. Error")
        error <- as.matrix(matrix_initial[[mat]]) - matrix_new[[mat]]
        decomp_error <- decomp_error + norm(error, type = "F")^2
    }
    return(decomp_error)
}

#' Generate random integer matrices with full column rank
#'
#' Function to generate random matrices with user-specified columns
#' and rows
#'
#' @param n number of matrices to generate. Integer value. Default 3.
#' @param ncols number of columns for each matrix. All matrices will have
#' the same number of columns. Integer value. Default 3.
#' @param nrows possible values for the number of rows for the matrices.
#' The number of rows should be greater than the column value.
#' Default c(3,4,5,6).
#' @param max_iter maximum number of iteration. Default 1000
#'
#' @return a list containing random matrices
#' @export
#'
#' @examples
#' # generate 4 matrices with possible rows = c(4,8) with 3 columns
#' set.seed(123)
#' mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:8)
createRandomMatrices <- function(n = 3, ncols = 3, nrows = 3:6,
                                 max_iter = 1000) {
    if (!is.numeric(n) | !is.numeric(ncols) | !is.numeric(nrows)) {
        stop("Integer n, ncols, and nrows values expected")
    }
    if (!all(nrows == floor(nrows))) {
        stop("Integer values for rows required")
    }
    if (length(ncols) != 1) {
        stop("All matrices should have same number of columns")
    }
    if (ncols > max(nrows)) {
        stop("# of rows should be >= # of columns")
    }
    if (ncols < max(nrows) &  ncols > min(nrows)) {
        nrows <- nrows[nrows >= ncols]
    }
    i <- k <- 0
    matrix_list <- list()
    while (i <= n) {
        nc <- ncols
        if (length(nrows) > 1)
            nr <- sample(nrows, size = 1)
        else
            nr <- nrows
        if (nr * nc <= 100)
            max_size <- 100
        else
            max_size <- 10 * nr * nc
        mat1 <- matrix(sample(1:max_size, size = nr * nc), nrow = nr, ncol = nc)
        if (qr(mat1)$rank == ncols) {
            k <- k + 1
            matrix_list[[paste0("mat", k)]] <- mat1
        }
        if (k == n) break
        if (i == max_iter) {
            cat("\n# of matrices created:", length(matrix_list))
            cat("\nmax iteration reached. Try increasing max_iter")
            stop("Exiting")
        }
    }
    return(matrix_list)
}
