
<!-- README.md is generated from README.Rmd. Please edit that file -->

## SBF: A R package for Shared Basis Factorization

Shared Basis Factorization (SBF) is a joint matrix diagonalization
approach we developed for cross-species gene expression analysis.
Approximate Shared Basis Factorization (A-SBF) is an extension of the
SBF approach.

### Installation

-   Clone from Github

``` r
git clone https://github.com/amalthomas111/SBF.git
```

Inside an R console

    library(devtools)
    install("<path to SBF>/SBF")

-   \[OR\] install directly from Github via `devtools`

``` r
library(devtools)
devtools::install_github("amalthomas111/SBF", build_vignettes = TRUE)
```

Please contact us via e-mail or submit a [GitHub
issue](https://github.com/amalthomas111/SBF/issues) if there is any
trouble with installation.

### Quick demo

### Load SBF package

``` r
# load SBF package
library(SBF)
```

### Analysis for a test dataset

-   We will first create a test dataset. SBF package has the function
    `createRandomMatrices` to create matrices with full column rank and
    different number of rows. All the matrices will have the same number
    of columns.

``` r
set.seed(1231)
mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)
sapply(mymat, dim)
#>      mat1 mat2 mat3 mat4
#> [1,]    5    6    4    5
#> [2,]    3    3    3    3
```

``` r
sapply(mymat, function(x) {qr(x)$rank})
#> mat1 mat2 mat3 mat4 
#>    3    3    3    3
```

We will use the test dataset as our *D*<sub>*i*</sub> matrices.

### SBF computation

Estimating *V* using the sum of
*D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub>.

``` r
# SBF call. Estimate V using the sum of Di^TDi
sbf <- SBF(matrix_list = mymat)
```

`?SBF` help function shows all arguments for the SBF call.

``` r
names(sbf)
#> [1] "v"      "lambda" "u"      "delta"  "m"
```

Check whether the estimated *V* is orthogonal.

``` r
zapsmall(sbf$v %*% t(sbf$v))
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```

Let us check the factorization error.

``` r
# calculate decomposition error
decomperror <- calcDecompError(mymat, sbf$u, sbf$delta, sbf$v)
decomperror
#> [1] 1.851693e-26
```

Estimating *V* using sum of inverse-variance weighted
*D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub>.

``` r
# SBF call. Estimate V using inverse-variance weighted Di^TDi
sbf <- SBF(matrix_list = mymat, weighted = TRUE)
# calculate decomposition error
decomperror <- calcDecompError(mymat, sbf$u, sbf$delta, sbf$v)
decomperror
#> [1] 1.894292e-26
```

SBF computation based on inter-sample correlation. *V* is estimated
using sum of *R*<sub>*i*</sub><sup>*T*</sup>*R*<sub>*i*</sub>.

``` r
# SBF call using correlation matrix
sbf_cor <- SBF(matrix_list = mymat, transform_matrix = TRUE)
decomperror <- calcDecompError(mymat, sbf_cor$u, sbf_cor$delta, sbf_cor$v)
decomperror
#> [1] 2.835482e-26
```

### Approximate-SBF (A-SBF)

Estimating *V* using the sum of
*D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub> and estimating
*U*<sub>*i*</sub>’s such that columns are orthonormal.

``` r
# A-SBF call
asbf <- SBF(matrix_list = mymat, approximate = TRUE)
```

``` r
names(asbf)
#> [1] "v"       "lambda"  "u"       "u_ortho" "delta"   "m"       "error"
```

In A-SBF, *V* is orthogonal and columns of *U*<sub>*i*</sub>’s are
orthonormal (*U*<sub>*i*</sub><sup>*T*</sup>*U*<sub>*i*</sub> = *I*).

``` r
zapsmall(t(asbf$v) %*% asbf$v)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```

``` r
# check the columns of first matrix of U_ortho
zapsmall(t(asbf$u_ortho[[names(asbf$u_ortho)[[1]]]]) %*%
           asbf$u_ortho[[names(asbf$u_ortho)[[1]]]])
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```

A-SBF is not an exact factorization.

``` r
# calculate decomposition error
decomperror <- calcDecompError(mymat, asbf$u_ortho, asbf$delta, asbf$v)
decomperror
#> [1] 2329.73
```

This error is already computed and stored in `asbf$error`.

A-SBF call with inverse variance weighting. Estimating *V* using the sum
of inverse-variance weighted
*D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub> and estimating
*U*<sub>*i*</sub>’s such that the columns are orthonormal.

``` r
# A-SBF call with inverse variance weighting
asbf_inv <- SBF(matrix_list = mymat, weighted = TRUE, approximate = TRUE)
# calculate decomposition error
decomperror_inv <- asbf_inv$error
decomperror_inv
#> [1] 1651.901
```

A-SBF computation based on inter-sample correlation. *V* is estimated
using sum of *R*<sub>*i*</sub><sup>*T*</sup>*R*<sub>*i*</sub>.

``` r
# A-SBF call using correlation matrix
asbf_cor <- SBF(matrix_list = mymat, approximate = TRUE,
                transform_matrix = TRUE)
# calculate decomposition error
decomperror_cor <- asbf_cor$error
decomperror_cor
#> [1] 14045.99
```

### Reduce A-SBF factorization error

Let us optimize the factorization error using the
`optimizeFactorization` function, for the three cases of A-SBF
computation shown before. Depending upon the mymat and initial values,
optimization could take some time.

``` r
myopt <- optimizeFactorization(mymat, asbf$u_ortho, asbf$delta, asbf$v)
myopt_inv <- optimizeFactorization(mymat, asbf_inv$u_ortho, asbf_inv$delta,
                                   asbf_inv$v)
myopt_cor <- optimizeFactorization(mymat, asbf_cor$u_ortho, asbf_cor$delta,
                                   asbf_cor$v)
```

The number of iteration taken for optimizing and new factorization
error:

``` r
cat("For asbf, # iteration =", myopt$error_pos, "final error =", myopt$error)
#> For asbf, # iteration = 220 final error = 1411.555
cat("For asbf inv, # iteration =", myopt_inv$error_pos, "final error =",
    myopt_inv$error)
#> For asbf inv, # iteration = 202 final error = 1411.555
cat("For asbf cor, # iteration =", myopt_cor$error_pos, "final error =",
    myopt_cor$error)
#> For asbf cor, # iteration = 196 final error = 1411.555
```

After optimization, for all three A-SBF factorizations, the final error
is

    #> same (up to 2 decimals). The final error is 1411.56

### Analysis using a test gene expression dataset

-   SBF package has a sample gene expression data with the mean
    expression of nine tissues in five species. Now we will use this
    data as our *D*<sub>*i*</sub> matrices.

``` r
# load sample dataset from SBF package
avg_counts <- SBF::TissueExprSpecies
# dimension of matrices for different species
sapply(avg_counts, dim)
#>      Homo_sapiens Macaca_mulatta Mus_musculus
#> [1,]        58676          30807        54446
#> [2,]            5              5            5
# names of different tissue types in humans
names(avg_counts[["Homo_sapiens"]])
#> [1] "hsapiens_brain"  "hsapiens_heart"  "hsapiens_kidney" "hsapiens_liver" 
#> [5] "hsapiens_testis"
```

SBF computation based on inter-sample correlation. *V* is estimated
using the mean *R*<sub>*i*</sub><sup>*T*</sup>*R*<sub>*i*</sub>.

``` r
# SBF call using correlation matrix
sbf_gem <- SBF(matrix_list = avg_counts, transform_matrix = TRUE)
decomperror <- calcDecompError(avg_counts, sbf_gem$u, sbf_gem$delta, sbf_gem$v)
decomperror
#> [1] 9.938196e-25
```

``` r
# A-SBF call using correlation matrix
asbf_gem <- SBF(matrix_list = avg_counts, approximate = TRUE,
                transform_matrix = TRUE)
asbf_gem$error
#> [1] 65865.92
```

For gene-expression analysis, if we want the shared space to represent
inter-sample correlation relationship, we do not update/change *V* while
optimizing the factorization error. In such cases, while reducing the
factorization error we set `optimizeV = FALSE` in the
`optimizeFactorization` function.

``` r
# optimize by keeping the V capturing R_i relationship
asbf_gem_opt <- optimizeFactorization(avg_counts, asbf_gem$u_ortho,
                                      asbf_gem$delta, asbf_gem$v,
                                      optimizeV = FALSE)
asbf_gem_opt$error
#> [1] 63540.08
```

### Contacts

***Amal Thomas*** *<amalthom@usc.edu>*

***Andrew D Smith*** *<andrewds@usc.edu>*

Contact us via e-mail or submit a [GitHub
issue](https://github.com/amalthomas111/SBF/issues)

### Copyright

Copyright (C) 2018-2021 Amal Thomas and Andrew D. Smith

Authors: Amal Thomas and Andrew D. Smith

SBF is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

SBF is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.
