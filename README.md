
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
devtools::install_github("amalthomas111/SBF")
# load package
library(SBF)
```

Please contact us via e-mail or through a [GitHub
issue](https://github.com/amalthomas111/SBF/issues) if there is any
trouble with installation.

### Quick demo

### load SBF package

``` r
# load SBF package
library(SBF)
```

### load test dataset

-   SBF package has a sample gene expression data with the mean
    expression of nine tissues in five species

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

### SBF computation

Estimating V using the sum of
*D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub>

``` r
# SBF call. Estimate V using the sum of Di^TDi
avg_counts <- SBF::TissueExprSpecies
sbf <- SBF(matrix_list = avg_counts)
```

`?SBF` help function shows all arguments for the SBF call.

``` r
names(sbf)
#> [1] "v"      "lambda" "u"      "delta"  "m"
```

Check whether the estimated *V* is orthogonal

``` r
zapsmall(sbf$v %*% t(sbf$v))
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    0    0    1    0    0
#> [4,]    0    0    0    1    0
#> [5,]    0    0    0    0    1
```

Let us check the factorization error.

``` r
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, sbf$u, sbf$delta, sbf$v)
decomperror
#> [1] 3.251627e-25
```

Estimating V using inverse-variance weighted
*D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub>

``` r
# SBF call. Estimate V using inverse-variance weighted Di^TDi
avg_counts <- SBF::TissueExprSpecies
sbf <- SBF(matrix_list = avg_counts, weighted = TRUE)
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, sbf$u, sbf$delta, sbf$v)
decomperror
#> [1] 4.694125e-25
```

SBF computation based on inter-sample correlation. *V* is estimated
using *R*<sub>*i*</sub><sup>*T*</sup>*R*<sub>*i*</sub>.

``` r
# SBF call using correlation matrix
avg_counts <- SBF::TissueExprSpecies
sbf_cor <- SBF(matrix_list = avg_counts, transform_matrix = TRUE)
decomperror <- calcDecompError(avg_counts, sbf_cor$u, sbf_cor$delta, sbf_cor$v)
decomperror
#> [1] 9.938196e-25
```

### Approximate-SBF (A-SBF)

Estimating V using the sum of
*D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub> and estimating
orthogonal *U*<sub>*i*</sub>’s such that columns are orthonormal.

``` r
# A-SBF call
avg_counts <- SBF::TissueExprSpecies
asbf <- SBF(matrix_list = avg_counts, approximate = TRUE)
```

``` r
names(asbf)
#> [1] "v"       "lambda"  "u"       "u_ortho" "delta"   "m"       "error"
```

In A-SBF, *V* is orthogonal and columns of *U*<sub>*i*</sub>’s are
orthonormal (*U*<sub>*i*</sub><sup>*T*</sup>*U*<sub>*i*</sub> = *I*).

``` r
zapsmall(t(asbf$v) %*% asbf$v)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    0    0    1    0    0
#> [4,]    0    0    0    1    0
#> [5,]    0    0    0    0    1
```

``` r
zapsmall(t(asbf$u_ortho[[names(asbf$u_ortho)[[1]]]]) %*%
           asbf$u_ortho[[names(asbf$u_ortho)[[1]]]])
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    0    0    1    0    0
#> [4,]    0    0    0    1    0
#> [5,]    0    0    0    0    1
```

A-SBF is not an exact factorization.

``` r
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, asbf$u_ortho, asbf$delta, asbf$v)
decomperror
#> [1] 7131.625
```

This error is already computed and stored in `asbf$error`.

A-SBF call with inverse variance weighting. Estimating V using the sum
of inverse-variance weighted
*D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub> and estimating
*U*<sub>*i*</sub>’s such that columns are orthonormal.

``` r
# A-SBF call with inverse variance weighting
avg_counts <- SBF::TissueExprSpecies
asbf_inv <- SBF(matrix_list = avg_counts, weighted = TRUE, approximate = TRUE)
# calculate decomposition error
decomperror_inv <- asbf_inv$error
decomperror_inv
#> [1] 7270.009
```

A-SBF computation based on inter-sample correlation. *V* is estimated
using *R*<sub>*i*</sub><sup>*T*</sup>*R*<sub>*i*</sub>.

``` r
# A-SBF call using correlation matrix
avg_counts <- SBF::TissueExprSpecies
asbf_cor <- SBF(matrix_list = avg_counts, approximate = TRUE,
                transform_matrix = TRUE)
# calculate decomposition error
decomperror_cor <- asbf_cor$error
decomperror_cor
#> [1] 65865.92
```

### Reduce A-SBF factorization error

Let us optimize the factorization for the three cases of A-SBF using the
`optimizeFactorization` function. Depending upon the initial values,
optimization could take some time.

``` r
myopt <- optimizeFactorization(avg_counts, asbf$u_ortho, asbf$delta, asbf$v)
myopt_inv <- optimizeFactorization(avg_counts, asbf_inv$u_ortho, asbf_inv$delta,
                                   asbf_inv$v)
myopt_cor <- optimizeFactorization(avg_counts, asbf_cor$u_ortho, asbf_cor$delta,
                                   asbf_cor$v)
```

The number of iteration taken for optimizing and new factorization
error:

``` r
cat("\nFor asbf, # iteration =", myopt$error_pos, "final error =", myopt$error)
#> 
#> For asbf, # iteration = 952 final error = 7081.105
cat("\nFor asbf inv, # iteration =", myopt_inv$error_pos, "final error =",
    myopt_inv$error)
#> 
#> For asbf inv, # iteration = 883 final error = 7081.105
cat("\nFor asbf cor, # iteration =", myopt_cor$error_pos, "final error =",
    myopt_cor$error)
#> 
#> For asbf cor, # iteration = 1690 final error = 7081.105
```

    #> 
    #> For all three factorizations, after optimizing the final errors is the same(up to 2 decimals)
    #> The final error is 7081.11

### Contacts

***Amal Thomas*** *<amalthom@usc.edu>*

***Andrew D Smith*** *<andrewds@usc.edu>*

Contact us via e-mail or through a [GitHub
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
