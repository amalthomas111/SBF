#' Compute GSVD factorization
#'
#' Function to GSVD factorization for two matrices D1 and D2
#' @param D1 A numeric matrix for GSVD factorization.
#' @param D2 A numeric matrix for GSVD factorization.
#'
#' @return a list containing u1, u2, v, d1, and d2
#' @export
#'
#' @examples
#' set.seed(1231)
#' mymat <- createRandomMatrices(n = 2, ncols = 3, nrows = 4:6)
#' gsvd <- GSVD(mymat$mat1, mymat$mat2)
#' print(gsvd$v)
GSVD <- function(D1, D2) {
  if (!is.numeric(D1) | !is.numeric(D2)) {
    stop("D1 and D2 must be numeric matrices")
  }
  if (ncol(D1) != ncol(D2))
    stop("Number of columns should be same")
  D1 <- as.matrix(D1)
  D2 <- as.matrix(D2)
  QR1 <- qr(D1)
  QR2 <- qr(D2)
  R1 <- qr.R(QR1)
  R2 <- qr.R(QR2)
  R <-  R1 %*% solve(R2)
  sv <- svd(R)
  V_notunit <- t(R1) %*% sv$u
  V <- sweep(V_notunit, 2, sqrt(colSums(V_notunit^2)), "/")

  # solving for  U1 and U2
  B1 <- D1 %*% solve(t(V))
  B2 <- D2 %*% solve(t(V))
  sigma1 <- sqrt(colSums(B1^2))
  sigma2 <- sqrt(colSums(B2^2))
  U1 <- sweep(B1, 2, sigma1, "/")
  U2 <- sweep(B2, 2, sigma2, "/")
  return(list(u1 = U1, u2 = U2, v = V, d1 = sigma1, d2 = sigma2))
}

#' Compute HOGSVD factorization
#'
#' Function to compute HOGSVD factorization
#' @param matrix_list A list containing Di matrices for joint matrix
#' factorization. Column names of each Di matrix may or may not have information
#' about tissue or cell type.
#' @param check_col_matching if the column names have information about tissue/
#' cell type and one-to-one correspondence of tissue types across species has to
#' be checked, set this parameter to be TRUE. Default FALSE.
#' @param col_sep separator in column names to separate different fields.
#' Example for column names 'hsapiens_brain', 'hsapiens_heart', etc., the
#' separator is underscore.
#' Set it to NULL if column matching across species has to be
#' performed and there is no separator in the column names.
#' Only checked if check_col_matching = TRUE. Default underscore.
#' @param col_index If a separator separates information in column names,
#' the col_index is the index in the column name corresponding to tissue or
#' cell type. E.g. for column name 'hsapiens_brain', col_index is 2.
#' Only checked if check_col_matching = TRUE. Default NULL.
#' @param verbose if TRUE print verbose lines. Default FALSE.
#'
#' @return a list containing u, v, lambda, d, and other outputs
#' of HOGSVD factorization.
#' @export
#'
#' @examples
#' set.seed(1231)
#' mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)
#' hogsvd <- HOGSVD(matrix_list = mymat)
#' print(hogsvd$v)
HOGSVD <- function(matrix_list = NULL, check_col_matching = FALSE,
                   col_index = NULL, col_sep = "_", verbose = FALSE) {
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
    A <- list()
    matrix_names <- names(matrix_list)
    for (mat in matrix_names) {
      A[[mat]] <- t(as.matrix(matrix_list[[mat]])) %*%
        as.matrix(matrix_list[[mat]])
    }
    S <- matrix(0L, nrow = nrow(A[[matrix_names[1]]]),
                ncol = ncol(A[[matrix_names[1]]]))
    for (i in 1:length(matrix_names)) {
      for (j in 1:length(matrix_names)) {
        if (j > i) {
          sum1 <- A[[matrix_names[i]]] %*% solve(A[[matrix_names[j]]])
          sum2 <- A[[matrix_names[j]]] %*% solve(A[[matrix_names[i]]])
          S <- S + sum1 + sum2
        }
      }
    }
    S <- S / (length(matrix_names) * (length(matrix_names) - 1))
    ev <- eigen(S)
    lambda <- ev$values
    V <- ev$vectors
    B <- sigma <- U <- list()
    for (mat in matrix_names) {
      B[[mat]] <- as.matrix(matrix_list[[mat]]) %*% solve(t(V))
      sigma[[mat]] <- sqrt(colSums(B[[mat]]^2))
      U[[mat]] <- sweep(B[[mat]], 2, sigma[[mat]], "/")
      if (verbose) {
        cat("\nmat:", mat)
        cat("\nDim of avg counts:", dim(matrix_list[[mat]]))
        cat("\nDim of V:", dim(V))
        cat("\nDim of B:", dim(B[[mat]]))
        cat("\nSigma vec:", sigma[[mat]])
        if (length(sigma[[mat]]) == 1)
          cat("\nDim of sigma:", dim(as.matrix(diag(as.matrix(sigma[[mat]])))))
        else
          cat("\nDim of sigma:", dim(as.matrix(diag(sigma[[mat]]))))
      }
    }
    return(list(v = V, lambda = lambda, u = U, d = sigma, s = S))
    } else {
    stop("\nInvalid matrix list. It should be a list. Exiting!")
  }
}
