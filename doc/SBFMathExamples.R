## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# load SBF package
library(SBF)

## -----------------------------------------------------------------------------
set.seed(1231)
mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)
sapply(mymat, dim)

## -----------------------------------------------------------------------------
sapply(mymat, function(x) {
  qr(x)$rank
  })

## -----------------------------------------------------------------------------
sbf <- SBF(matrix_list = mymat)
sbf_inv <- SBF(matrix_list = mymat, weighted = TRUE)
sbf_cor <- SBF(matrix_list = mymat, transform_matrix = TRUE)

## -----------------------------------------------------------------------------
names(sbf)

## -----------------------------------------------------------------------------
sbf$v

## -----------------------------------------------------------------------------
printDelta <- function(l) {
  for (eachmat in names(l$delta)) {
  cat(eachmat, ":", l$delta[[eachmat]], "\n")
  }
}
cat("sbf\n");printDelta(sbf)
cat("sbf_inv\n");printDelta(sbf_inv)
cat("sbf_cor\n");printDelta(sbf_cor)

## -----------------------------------------------------------------------------
zapsmall(t(sbf$v) %*% sbf$v)

## -----------------------------------------------------------------------------
qr(sbf$v)$rank

## -----------------------------------------------------------------------------
sapply(sbf$u, dim)

## -----------------------------------------------------------------------------
t(sbf$u[[names(sbf$u)[1]]]) %*% sbf$u[[names(sbf$u)[1]]]

## -----------------------------------------------------------------------------
sbf$lambda

## -----------------------------------------------------------------------------
calcDecompError(mymat, sbf$u, sbf$delta, sbf$v)
calcDecompError(mymat, sbf_inv$u, sbf_inv$delta, sbf_inv$v)
calcDecompError(mymat, sbf_cor$u, sbf_cor$delta, sbf_cor$v)

## -----------------------------------------------------------------------------
sapply(mymat, function(x) sum(diag(cov(x))))

## -----------------------------------------------------------------------------
mat5 <- matrix(c(130, 183, 62, 97, 147, 94, 102, 192, 19), byrow = T,
                    nrow = 3, ncol = 3)
mat5_highvar <- matrix(c(406, 319, 388, 292, 473, 287, 390, 533, 452),
                       byrow = T, nrow = 3, ncol = 3)

mymat_new <- mymat
mymat_new[["mat5"]] <- mat5
sapply(mymat_new, function(x) sum(diag(cov(x))))
mymat_new_noisy <- mymat
mymat_new_noisy[["mat5"]] <- mat5_highvar
sapply(mymat_new_noisy, function(x) sum(diag(cov(x))))

## -----------------------------------------------------------------------------
sbf_new <- SBF(matrix_list = mymat_new)
sbf_inv_new <- SBF(matrix_list = mymat_new, weighted = TRUE)


sbf_new_noisy <- SBF(matrix_list = mymat_new_noisy)
sbf_inv_new_noisy <- SBF(matrix_list = mymat_new_noisy, weighted = TRUE)

## -----------------------------------------------------------------------------
e1 <- calcDecompError(mymat, sbf_new$u[1:4], sbf_new$delta[1:4], sbf_new$v)
e2 <- calcDecompError(mymat, sbf_new_noisy$u[1:4], sbf_new_noisy$delta[1:4],
                      sbf_new_noisy$v)
e2 / e1

## -----------------------------------------------------------------------------
e3 <- calcDecompError(mymat, sbf_inv_new$u[1:4], sbf_inv_new$delta[1:4],
                      sbf_inv_new$v)
e4 <- calcDecompError(mymat, sbf_inv_new_noisy$u[1:4],
                      sbf_inv_new_noisy$delta[1:4], sbf_inv_new_noisy$v)
e4 / e3

## -----------------------------------------------------------------------------
set.seed(1231)
mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)
sapply(mymat, dim)

## -----------------------------------------------------------------------------
asbf <- SBF(matrix_list = mymat, orthogonal = TRUE,
            minimizeError = FALSE)
asbf_inv <- SBF(matrix_list = mymat, weighted = TRUE, orthogonal = TRUE,
                minimizeError = FALSE)
asbf_cor <- SBF(matrix_list = mymat, orthogonal = TRUE,
                transform_matrix = TRUE, minimizeError = FALSE)

## -----------------------------------------------------------------------------
names(asbf)

## -----------------------------------------------------------------------------
# decomposition error
asbf$error
asbf_inv$error
asbf_cor$error

## -----------------------------------------------------------------------------
zapsmall(t(asbf$u_ortho[[names(asbf$u_ortho)[1]]]) %*%
           asbf$u_ortho[[names(asbf$u_ortho)[1]]])

## -----------------------------------------------------------------------------
zapsmall(t(asbf$v) %*% asbf$v)

## -----------------------------------------------------------------------------
set.seed(1231)
mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)
asbf <- SBF(matrix_list = mymat, orthogonal = TRUE)
asbf_inv <- SBF(matrix_list = mymat, weighted = TRUE, orthogonal = TRUE)
asbf_cor <- SBF(matrix_list = mymat, orthogonal = TRUE, transform_matrix = TRUE)

## -----------------------------------------------------------------------------
names(asbf)

## -----------------------------------------------------------------------------
# initial decomposition error
asbf$error_start
asbf_inv$error_start
asbf_cor$error_start

## -----------------------------------------------------------------------------
# final decomposition error
asbf$error
asbf_inv$error
asbf_cor$error

## ---- echo = FALSE------------------------------------------------------------
if ((round(asbf$error, 2) == round(asbf_inv$error, 2)) &&
    (round(asbf$error, 2) == round(asbf_cor$error, 2))) {
  cat("same (up to 2 decimals). The final error is", round(asbf$error, 2))
  } else {
  cat("not exactly the same")
    }

## -----------------------------------------------------------------------------
myopt <- optimizeFactorization(mymat, asbf$u_ortho_start, asbf$delta_start,
                               asbf$v_start)
names(myopt)

## -----------------------------------------------------------------------------
myopt$error

## -----------------------------------------------------------------------------
cat("For asbf, # iteration =", asbf$error_pos, "final error =", asbf$error)
cat("For asbf inv, # iteration =", asbf_inv$error_pos, "final error =",
    asbf_inv$error)
cat("For asbf cor, # iteration =", asbf_cor$error_pos, "final error =",
    asbf_cor$error)

## -----------------------------------------------------------------------------
set.seed(1231)
mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)

## -----------------------------------------------------------------------------
set.seed(111)
rand_mat <- createRandomMatrices(n = 1, ncols = 3, nrows = 3)
cat("\nRank is:", qr(rand_mat[[1]])$rank, "\n")
dim(rand_mat[[1]])

## -----------------------------------------------------------------------------
mysvd <- svd(rand_mat[[1]])
randV <- mysvd$v

## -----------------------------------------------------------------------------
# get Ui and Delta for this newV
out <- computeUDelta(mymat, randV)
names(out)

## -----------------------------------------------------------------------------
calcDecompError(mymat, out$u_ortho, out$d, randV)

## -----------------------------------------------------------------------------
newopt <- optimizeFactorization(mymat, out$u_ortho, out$d, randV)
# Number of updates taken
newopt$error_pos
# New error
newopt$error

## -----------------------------------------------------------------------------
mysvd <- svd(rand_mat[[1]])
randV <- mysvd$u
dim(randV)

## -----------------------------------------------------------------------------
# get Ui and Delta for this newV
out <- computeUDelta(mymat, randV)
calcDecompError(mymat, out$u_ortho, out$d, randV)

## -----------------------------------------------------------------------------
newopt <- optimizeFactorization(mymat, out$u_ortho, out$d, randV)
# Number of updates taken
newopt$error_pos
# New error
newopt$error

## -----------------------------------------------------------------------------
set.seed(111)
# new random v
newv <- createRandomMatrices(n = 1, ncols = 3, nrows = 3)[[1]]
# seed value
k <- 2392
newu <- newd <- list()
for (i in names(mymat)) {
  myrow <- nrow(mymat[[i]])
  mycol <- ncol(mymat[[i]])
  set.seed(k)
  # new random u_i
  newu[[i]] <- createRandomMatrices(n = 1, ncols = mycol, nrows = myrow)[[1]]
  set.seed(k * 2)
  # new random d_i
  newd[[i]] <- sample(1:1000, size = mycol)
  newmat <- newu[[i]] %*% diag(newd[[i]]) %*% t(newv)
  if (!qr(newmat)$rank == mycol)
    cat("\nNew matrix does not have full column rank")
  k <- k + 1
}
error <- calcDecompError(mymat, newu, newd, newv)
cat("\nInitial error = ", error, "\n")

## -----------------------------------------------------------------------------
newopt <- optimizeFactorization(mymat, newu, newd, newv)
newopt$error_pos
newopt$error

## -----------------------------------------------------------------------------
set.seed(171)
newmat <- createRandomMatrices(n = 1, ncols = 3, nrows = 3)
newmat

## -----------------------------------------------------------------------------
newu <- newd <- list()
newu[["mat1"]]  <- diag(3)
newd[["mat1"]] <- diag(newu[["mat1"]])
newu
newd

## -----------------------------------------------------------------------------
calcDecompError(newmat, newu, newd, diag(3))

## -----------------------------------------------------------------------------
opt_new <- optimizeFactorization(newmat, newu, newd, diag(3))
cat("\n # of updates:", opt_new$error_pos, "\n")
opt_new$error

## -----------------------------------------------------------------------------
newmat
opt_new$u[[1]] %*% diag(opt_new$d[[1]]) %*% t(opt_new$v)

## -----------------------------------------------------------------------------
opt_new1 <- optimizeFactorization(newmat, newu, newd, diag(3), tol = 1e-21)
cat("\n # of updates:", opt_new1$error_pos, "\n")
opt_new1$error

## -----------------------------------------------------------------------------
newmat
opt_new1$u[[1]] %*% diag(opt_new1$d[[1]]) %*% t(opt_new1$v)

## -----------------------------------------------------------------------------
newmat_svd <- svd(newmat[[1]])

## -----------------------------------------------------------------------------
newmat_svd$d
opt_new1$d

## -----------------------------------------------------------------------------
newmat_svd$u
opt_new1$u[[1]]

## -----------------------------------------------------------------------------
newmat_svd$v
opt_new1$v

## -----------------------------------------------------------------------------
set.seed(253)
randmat_new <- createRandomMatrices(n = 1, ncols = 3, nrows = 3)
randmat_new
newsvd <- svd(randmat_new[[1]])

## -----------------------------------------------------------------------------
newu <- newd <- list()
newu[[names(randmat_new)]] <- newsvd$u
newd[[names(randmat_new)]] <- newsvd$d

## -----------------------------------------------------------------------------
calcDecompError(newmat, newu, newd, newsvd$v)

## -----------------------------------------------------------------------------
opt_new <- optimizeFactorization(newmat, newu, newd, newsvd$v)
cat("\n # of updates:", opt_new$error_pos, "\n")
opt_new$error

## -----------------------------------------------------------------------------
newmat
opt_new$u[[1]] %*% diag(opt_new$d[[1]]) %*% t(opt_new$v)

## -----------------------------------------------------------------------------
opt_new1 <- optimizeFactorization(newmat, newu, newd, newsvd$v, tol = 1e-21)
cat("\n # of updates:", opt_new1$error_pos, "\n")
opt_new1$error
opt_new1$u[[1]] %*% diag(opt_new1$d[[1]]) %*% t(opt_new1$v)

## -----------------------------------------------------------------------------
newmat_svd <- svd(newmat[[1]])

## -----------------------------------------------------------------------------
newmat_svd$u
opt_new1$u[[1]]

## -----------------------------------------------------------------------------
newmat_svd$v
opt_new1$v

