
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

#### load test dataset

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

#### SBF computation

``` r
# SBF call. Estimate V using sum of Di^TDi
avg_counts <- SBF::TissueExprSpecies
sbf <- SBF(matrix_list = avg_counts, check_col_matching = TRUE, col_index = 2,
           weighted = FALSE, approximate = FALSE, transform_matrix = FALSE)
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, sbf$delta, sbf$u, sbf$v)
```

`?SBF` help function shows all arguments for the SBF call.

``` r
# SBF call. Estimate V using inverse-variance weighted Di^TDi
avg_counts <- SBF::TissueExprSpecies
sbf <- SBF(matrix_list = avg_counts, check_col_matching = TRUE, col_index = 2,
           weighted = TRUE, approximate = FALSE, transform_matrix = FALSE)
#> 
#> Inverse variance weighting applied
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, sbf$delta, sbf$u, sbf$v)
```

##### SBF computation based on inter-sample correlation

``` r
# SBF call using correlation matrix
avg_counts <- SBF::TissueExprSpecies
sbf_cor <- SBF(matrix_list = avg_counts, check_col_matching = TRUE,
               col_index = 2, weighted = FALSE,
               approximate = FALSE, transform_matrix = TRUE)
#> 
#> V is computed using inter-sample correlation
decomperror <- calcDecompError(avg_counts, sbf_cor$delta, sbf_cor$u, sbf_cor$v)
```

#### Approximate SBF (A-SBF) computation

``` r
# A-SBF call
avg_counts <- SBF::TissueExprSpecies
asbf <- SBF(matrix_list = avg_counts, check_col_matching = TRUE, col_index = 2,
            weighted = FALSE, approximate = TRUE, transform_matrix = FALSE)
#> 
#> A-SBF is computed
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, asbf$delta, asbf$u_ortho, asbf$v)
```

``` r
# A-SBF call with inverse variance weighting
avg_counts <- SBF::TissueExprSpecies
asbf <- SBF(matrix_list = avg_counts, check_col_matching = TRUE, col_index = 2,
            weighted = TRUE, approximate = TRUE, transform_matrix = FALSE)
#> 
#> Inverse variance weighting applied
#> 
#> A-SBF is computed
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, asbf$delta, asbf$u_ortho, asbf$v)
```

##### Approximate SBF (A-SBF) computation based on inter-sample correlation

``` r
# A-SBF call using correlation matrix
avg_counts <- SBF::TissueExprSpecies
asbf_cor <- SBF(matrix_list = avg_counts, check_col_matching = TRUE,
                col_index = 2, weighted = FALSE,
                approximate = TRUE, transform_matrix = TRUE)
#> 
#> V is computed using inter-sample correlation
#> 
#> A-SBF is computed
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, asbf_cor$delta, asbf_cor$u_ortho,
                                asbf_cor$v)
```

### Dependencies

-   `data.table`

### Contacts

***Amal Thomas*** *<amalthom@usc.edu>*

***Andrew D Smith*** *<andrewds@usc.edu>*

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
