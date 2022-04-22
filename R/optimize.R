#' Find the Nearest orthogonal matrix
#'
#' Function to compute nearest orthogonal matrix to a set of matrices.
#' @param mat_list A list of numeric matrices
#'
#' @return orthogonal matrix
#' @export
#'
#' @examples
#' set.seed(1231)
#' mymat <- createRandomMatrices(n = 3, ncols = 3, nrows = 3)
#' nearest_ortho <- calcNearestOrthoMatrix(mymat)
#' print(nearest_ortho %*% t(nearest_ortho))
calcNearestOrthoMatrix <- function(mat_list) {
  matrix_names <- names(mat_list)
  if (!all(sapply(mat_list, ncol) == ncol(mat_list[[matrix_names[1]]])))
    stop(paste0("\nAll matrices should have same number of columns"))
  if (!all(sapply(mat_list, nrow) == nrow(mat_list[[matrix_names[1]]])))
    stop(paste0("\nAll matrices should have same number of rows"))
  mat_sum <- matrix(0L, nrow =  nrow(mat_list[[matrix_names[1]]]),
                    ncol = nrow(mat_list[[matrix_names[1]]]))
  for (i in matrix_names) {
    mat_sum <- mat_sum + as.matrix(mat_list[[i]])
  }
  mysvd <- svd(mat_sum)
  return(mysvd$u %*% t(mysvd$v))
}

#' Compute U and Delta for a given set of matrices and a V
#'
#' Function to compute U and Delta for a given set of matrices and V.
#' @param mat_list A list containing numeric matrices
#' @param V A numeric matrix of dimension ncol(matrix1) * ncol(matrix1)
#' @param approximate Compute U with orthonormal columns. Default TRUE
#'
#' @return list of U and delta
#' @export
#'
#' @examples
#' set.seed(1231)
#' mymat <- createRandomMatrices(n = 3, ncols = 3, nrows = 3)
#' set.seed(325)
#' mymat_rand <- createRandomMatrices(n = 1, ncols = 3, nrows = 3)[[1]]
#' out <- computeUDelta(mymat, svd(mymat_rand)$v)
computeUDelta <- function(mat_list, V, approximate = TRUE) {
  matrix_names <- names(mat_list)
  B <- delta <- U <- U_ortho <- list()
  for (mat in matrix_names) {
    B[[mat]] <- as.matrix(mat_list[[mat]]) %*% V
    delta[[mat]] <- sqrt(colSums(B[[mat]]^2))
    U[[mat]] <- sweep(B[[mat]], 2, delta[[mat]], FUN = "/")
    if (approximate) {
      if (length(delta[[mat]]) == 1) {
        phi <- as.matrix(diag(as.matrix(delta[[mat]])))
      } else {
        phi <- as.matrix(diag(delta[[mat]]))
      }
      dvsig_svd <- svd(as.matrix(mat_list[[mat]]) %*% V %*% phi)
      U_ortho[[mat]] <- dvsig_svd$u %*% t(dvsig_svd$v)
      row.names(U_ortho[[mat]]) <- row.names(U[[mat]])
    }
  }
  if (approximate) {
    delta_new <- updateDelta(mat_list, U_ortho, V)
    error <- calcDecompError(mat_list, U_ortho, delta_new, V)
    return(list(u = U, u_ortho = U_ortho, d = delta, d_ortho = delta_new,
                error = error))
  } else {
    return(list(u = U, d = delta))
  }
}
#' Update U with orthonormal columns
#'
#' Function to update U using given delta and V. The new U is set as XY^T,
#' where SVD of D V Delta = X Sigma Y^T.
#' @param mat_list A list containing numeric matrices
#' @param d A list containing delta matrices
#' @param v V matrix
#'
#' @return a list containing new U matrices
#' @export
#'
#' @examples
#' set.seed(1231)
#' mymat <- createRandomMatrices(n = 3, ncols = 3, nrows = 3)
#' sbf <- SBF(matrix_list = mymat, approximate = FALSE,
#'            transform_matrix = FALSE)
#' newU <- updateU(mymat, sbf$delta, sbf$v)
updateU <- function(mat_list, d, v) {
  if (! all(names(mat_list) == names(d)))
    stop("Different names for matrix list and d list")
  u_ortho <- list()
  for (matrix_name in names(mat_list)) {
    if (length(d[[matrix_name]]) == 1) {
      phi <- as.matrix(diag(as.matrix(d[[matrix_name]])))
    } else {
      phi <- as.matrix(diag(d[[matrix_name]]))
    }
    mysvd <- svd(as.matrix(mat_list[[matrix_name]]) %*% v %*% phi)
    u_ortho[[matrix_name]] <- mysvd$u %*% t(mysvd$v)
  }
  return(u_ortho)
}

#' Update Delta
#'
#' Function to update Delta using given Uis and V. The new Delta is set as
#' diag(U^T D V).
#' @param mat_list A list containing numeric matrices
#' @param u A list containing U matrices
#' @param v V matrix
#'
#' @return a list containing new U matrices
#' @export
#'
#' @examples
#' set.seed(1231)
#' mymat <- createRandomMatrices(n = 3, ncols = 3, nrows = 3)
#' sbf <- SBF(matrix_list = mymat, approximate = FALSE,
#'            transform_matrix = FALSE)
#' newU <- updateU(mymat, sbf$delta, sbf$v)
#' newDelta <- updateDelta(mymat, newU, sbf$v)
updateDelta <- function(mat_list, u, v) {
  if (! all(names(mat_list) == names(u)))
    stop("Different names for matrix list and u list")
  delta_updated <- list()
  for (matrix_name in names(mat_list)) {
    delta_updated[[matrix_name]] <- diag(t(as.matrix(u[[matrix_name]])) %*%
                                    as.matrix(mat_list[[matrix_name]]) %*% v)
  }
  return(delta_updated)
}

#' Update V
#'
#' Function to update V using given Delta and U. The new V is set as
#' X Y^T, where SVD of D^T U Delta = X Sigma Y^T.
#' @param mat_list A list containing numeric matrices
#' @param d A list containing delta matrices
#' @param u A list containing U matrices
#'
#' @return a list containing new U matrices
#' @export
#'
#' @examples
#' set.seed(1231)
#' mymat <- createRandomMatrices(n = 3, ncols = 3, nrows = 3)
#' sbf <- SBF(matrix_list = mymat, approximate = FALSE,
#'            transform_matrix = FALSE)
#' newU <- updateU(mymat, sbf$delta, sbf$v)
#' newDelta <- updateDelta(mymat, newU, sbf$v)
#' newV <- updateV(mymat, newU, newDelta)
updateV <- function(mat_list, u, d) {
  matrix_names <- names(mat_list)
  if (! all(names(mat_list) == names(u)))
    stop("Different names for matrix list and u list")
  if (! all(names(mat_list) == names(d)))
    stop("Different names for matrix list and d list")
  mat_sum <- matrix(0L, nrow =  ncol(mat_list[[matrix_names[1]]]),
                    ncol = ncol(mat_list[[matrix_names[1]]]))
  for (name in matrix_names) {
    if (length(d[[name]]) == 1) {
      phi <- as.matrix(diag(as.matrix(d[[name]])))
    } else {
      phi <- as.matrix(diag(d[[name]]))
    }
    mat_sum <- mat_sum + (t(as.matrix(mat_list[[name]])) %*%
                            as.matrix(u[[name]]) %*% phi)
  }
  mysvd <- svd(mat_sum)
  v <-  mysvd$u %*% t(mysvd$v)
  return(v)
}
#' Function to minimize decomposition error
#'
#' Function iteratively updates u, delta, and V and finds the
#' minimum decomposition error.
#' @param mat_list A list containing numeric matrices
#' @param u A list containing U matrices
#' @param d A list containing delta matrices
#' @param v V matrix
#' @param optimizeV Whether initial V should be update or not. Default TRUE
#' @param initial_exact Whether the initial value of U, Delta,
#' and V gives exact factorization. Default FALSE
#' @param max_iter Maximum number of iterations. In each iteration u, d, and
#' v are updated. Default 1e4
#' @param tol Tolerance value. During the iterations, if the difference between
#' previous best and current best factorization error becomes less than tol,
#' no more iteration is performed. Default tol = 1e-10
#' @param verbose print error at each update. Default FALSE
#'
#' @return a list containing optimal U, delta, and V that
#' minimizes the factorization error
#' @export
#'
#' @examples
#' set.seed(1231)
#' mymat <- createRandomMatrices(n = 3, ncols = 3, nrows = 3)
#' sbf <- SBF(matrix_list = mymat, approximate = FALSE,
#'            transform_matrix = FALSE)
#' newU <- updateU(mymat, sbf$delta, sbf$v)
#' newDelta <- updateDelta(mymat, newU, sbf$v)
#' newV <- updateV(mymat, newU, newDelta)
#' opt <- optimizeFactorization(mymat, newU, newDelta, newV, max_iter = 1e4)
optimizeFactorization <- function(mat_list, u, d, v, optimizeV = TRUE,
                                  initial_exact = FALSE,
                                  max_iter = 1e4, tol = 1e-10,
                                  verbose = FALSE) {
  min_error <- calcDecompError(mat_list, u, d, v)
  if (initial_exact == TRUE) {
    u <- updateU(mat_list, d, v)
    min_error <- Inf
  }
  if (min_error < tol) {
    stop(paste("\nInitial factorization error is already < tol.",
                "Try setting initial_exact = TRUE"))
  }
  error_vec <- c()
  error_decreased <- FALSE
  k <- 0
  if (verbose) {
    cat("\t\n Update:", k, " Error:", round(min_error, 2))
    k <- k + 1
  }

  for (i in 1:max_iter) {
    d <- updateDelta(mat_list, u, v)
    error <- calcDecompError(mat_list, u, d, v)
    error_vec <- c(error_vec, error)
    if (verbose) {
      cat("\t\n Update:", k, " Error:", round(error, 2), "dt:",
          round(error - min_error, 2))
      k <- k + 1
    }
    if (error <= min_error) {
      u_opt <- u
      v_opt <- v
      d_opt <- d
      previous_min <- min_error
      min_error <- error
      error_decreased <- TRUE
      if (abs(min_error - previous_min) < tol) {
        break
      }
    }
    if (optimizeV == TRUE) {
      v <- updateV(mat_list, u, d)
      error <- calcDecompError(mat_list, u, d, v)
      error_vec <- c(error_vec, error)
      if (verbose) {
        cat("\t\n Update:", k, " Error:", round(error, 2), "dt:",
            round(error - min_error, 2))
        k <- k + 1
      }
      if (error <= min_error) {
        u_opt <- u
        v_opt <- v
        d_opt <- d
        previous_min <- min_error
        min_error <- error
        error_decreased <- TRUE
        if (abs(min_error - previous_min) < tol) {
          break
        }
      }
    }
    u <- updateU(mat_list, d, v)
    error <- calcDecompError(mat_list, u, d, v)
    error_vec <- c(error_vec, error)
    if (verbose) {
      cat("\t\n Update:", k, " Error:", round(error, 2), "dt:",
          round(error - min_error, 2))
      k <- k + 1
    }
    if (error <= min_error) {
      u_opt <- u
      v_opt <- v
      d_opt <- d
      previous_min <- min_error
      min_error <- error
      error_decreased <- TRUE
      if (abs(min_error - previous_min) < tol) {
        break
      }
    }
  }
  if (error_decreased == TRUE) {
    min_pos <- which.min(error_vec)
    if (min_pos == (3 * max_iter))
      cat("\nNot converged! Try increasing max_iter\n")
    return(list(u = u_opt,
                v = v_opt,
                d = d_opt,
                error = min_error,
                error_pos = min_pos,
                error_vec = error_vec))
  } else {
    cat("\nFatorization error could not be decreased further\n")
    return(NULL)
  }
}
