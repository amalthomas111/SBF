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
sapply(mymat, function(x) {qr(x)$rank})

## -----------------------------------------------------------------------------
sbf <- SBF(matrix_list = mymat, check_col_matching = FALSE, weighted = FALSE,
            approximate = FALSE, transform_matrix = FALSE)
sbf_inv <- SBF(matrix_list = mymat, check_col_matching = FALSE, weighted = TRUE,
           approximate = FALSE, transform_matrix = FALSE)
sbf_cor <- SBF(matrix_list = mymat, check_col_matching = FALSE, weighted = FALSE,
               approximate = FALSE, transform_matrix = TRUE)


## -----------------------------------------------------------------------------
names(sbf)

## -----------------------------------------------------------------------------
sbf$v
sbf_inv$v
sbf_cor$v

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
sbf_inv$lambda
sbf_cor$lambda

## -----------------------------------------------------------------------------
names(mymat)
sbf$delta

## -----------------------------------------------------------------------------
calcDecompError(mymat, sbf$delta, sbf$u, sbf$v)
calcDecompError(mymat, sbf_inv$delta, sbf_inv$u, sbf_inv$v)
calcDecompError(mymat, sbf_cor$delta, sbf_cor$u, sbf_cor$v)

## -----------------------------------------------------------------------------
sapply(mymat, function(x) sum(diag(cov(x))))

## -----------------------------------------------------------------------------
mat5 <- matrix(c(130, 183, 62, 97, 147, 94, 102, 192, 19), byrow = T,
                    nrow = 3, ncol = 3)
mat5_highvar <- matrix(c(406, 319, 388, 292, 473, 287, 390, 533, 452), byrow = T,
                    nrow = 3, ncol = 3)

mymat_new <- mymat
mymat_new[["mat5"]] <- mat5
sapply(mymat_new, function(x) sum(diag(cov(x))))
mymat_new_noisy <- mymat
mymat_new_noisy[["mat5"]] <- mat5_highvar
sapply(mymat_new_noisy, function(x) sum(diag(cov(x))))

## -----------------------------------------------------------------------------
sbf_new <- SBF(matrix_list = mymat_new, check_col_matching = FALSE,
               weighted = FALSE, approximate = FALSE, transform_matrix = FALSE)
sbf_inv_new <- SBF(matrix_list = mymat_new, check_col_matching = FALSE,
                   weighted = TRUE, approximate = FALSE,
                   transform_matrix = FALSE)


sbf_new_noisy <- SBF(matrix_list = mymat_new_noisy, check_col_matching = FALSE,
                     weighted = FALSE, approximate = FALSE,
                     transform_matrix = FALSE)
sbf_inv_new_noisy <- SBF(matrix_list = mymat_new_noisy,
                         check_col_matching = FALSE, weighted = TRUE,
                         approximate = FALSE, transform_matrix = FALSE)

## -----------------------------------------------------------------------------
e1 <- calcDecompError(mymat, sbf_new$delta[1:4], sbf_new$u[1:4], sbf_new$v)
e2 <- calcDecompError(mymat, sbf_new_noisy$delta[1:4], sbf_new_noisy$u[1:4],
                      sbf_new_noisy$v)
e2 / e1

## -----------------------------------------------------------------------------
e3 <- calcDecompError(mymat, sbf_inv_new$delta[1:4],
                      sbf_inv_new$u[1:4], sbf_inv_new$v)
e4 <- calcDecompError(mymat, sbf_inv_new_noisy$delta[1:4],
                      sbf_inv_new_noisy$u[1:4], sbf_inv_new_noisy$v)
e4 / e3

## -----------------------------------------------------------------------------
asbf <- SBF(matrix_list = mymat, check_col_matching = FALSE, weighted = FALSE,
            approximate = TRUE, transform_matrix = FALSE)
asbf_inv <- SBF(matrix_list = mymat, check_col_matching = FALSE, weighted = TRUE,
           approximate = TRUE, transform_matrix = FALSE)
asbf_cor <- SBF(matrix_list = mymat, check_col_matching = FALSE, weighted = FALSE,
                approximate = TRUE, transform_matrix = TRUE)

## -----------------------------------------------------------------------------
names(asbf)

## -----------------------------------------------------------------------------
asbf$error
asbf_inv$error
asbf_cor$error

## -----------------------------------------------------------------------------
calcDecompError(mymat, asbf$delta, asbf$u_ortho, asbf$v)
calcDecompError(mymat, asbf_inv$delta, asbf_inv$u_ortho, asbf_inv$v)
calcDecompError(mymat, asbf_cor$delta, asbf_cor$u_ortho, asbf_cor$v)

## -----------------------------------------------------------------------------
zapsmall(t(asbf$u_ortho[[names(asbf$u_ortho)[1]]]) %*%
           asbf$u_ortho[[names(asbf$u_ortho)[1]]])

## -----------------------------------------------------------------------------
zapsmall(t(asbf$v) %*% asbf$v)

## -----------------------------------------------------------------------------
# load dataset
avg_counts <- SBF::TissueExprSpecies
# check the names of species
names(avg_counts)

## -----------------------------------------------------------------------------
# head for first species
avg_counts[[names(avg_counts)[1]]][1:5, 1:5]

## -----------------------------------------------------------------------------
sapply(avg_counts, dim)

## -----------------------------------------------------------------------------
# A-SBF call using correlation matrix
asbf_cor <- SBF(matrix_list = avg_counts, col_index = 2, weighted = FALSE,
                approximate = TRUE, transform_matrix = TRUE)
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, asbf_cor$delta, asbf_cor$u_ortho,
                                asbf_cor$v)
decomperror

