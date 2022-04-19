
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
#> 
#> A-SBF optimizing factorization error
```

``` r
names(asbf)
#>  [1] "v"             "u"             "d"             "error"        
#>  [5] "error_pos"     "error_vec"     "v_start"       "lambda_start" 
#>  [9] "u_start"       "u_ortho_start" "delta_start"   "m"            
#> [13] "error_start"
```

For A-SBF, the factorization error is optimized by default
(`minimizeError=TRUE`) and the `optimizeFactorization` function is
invoked. Depending upon the data matrix (`mymat`) and initial values,
optimization could take some time.

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
zapsmall(t(asbf$u[[names(asbf$u)[[1]]]]) %*%
           asbf$u[[names(asbf$u)[[1]]]])
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```

A-SBF is not an exact factorization. The decomposition error is
minimized and stored in `asbf$error`.

``` r
# initial decomposition error
asbf$error_start
#> [1] 2329.73
# final decomposition error
asbf$error
#> [1] 1411.555
```

The number of iteration taken for optimizing:

``` r
cat("For asbf, # iteration =", asbf$error_pos, "final error =", asbf$error)
#> For asbf, # iteration = 220 final error = 1411.555
```

Detailed explanation of the math and examples case of factorization can
be found in the vignettes/docs directory.

### Factorization for an example gene expression dataset

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

A-SBF computation based on inter-sample correlation. *V* is estimated
using the mean *R*<sub>*i*</sub><sup>*T*</sup>*R*<sub>*i*</sub>.

``` r
# A-SBF call using correlation matrix
asbf_gem <- SBF(matrix_list = avg_counts, approximate = TRUE,
                transform_matrix = TRUE, tol = 1e-4)
#> 
#> A-SBF optimizing factorization error
# initial decomposition error
asbf_gem$error_start
#> [1] 65865.92
# final decomposition error
asbf_gem$error
#> [1] 7088.058
```

Note: For high-dimensional datasets, vary the tolerance threshold
(`tol`) or maximum number of iteration parameter (`max_iter`) to reduce
the computing time.
<!-- For gene-expression analysis, if we want the shared space to represent -->
<!-- inter-sample correlation relationship, we do not update/change $V$ -->
<!-- while optimizing the factorization error. -->
<!-- In such cases, while reducing the factorization error -->
<!-- we set `optimizeV = FALSE` in the `optimizeFactorization` function. -->
<!-- We now optimize the factorization error to find the closest space. -->

Different cross-species analysis examples using A-SBF can be found in
the vignettes/docs directory.

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
