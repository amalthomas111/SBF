## SBF: A R package for **S**hared **B**asis **F**actorization

**S**hared **B**asis **F**actorization (SBF) is a joint matrix diagonalization
approach we developed for cross-species gene expression analysis.
Approximate Shared Basis Factorization (A-SBF) is an extension of the SBF
approach.

### Installation
- Clone from Github
```
git clone https://github.com/amalthomas111/SBF.git
# Inside R
library(devtools)
install("<path to SBF>/SBF")
# load SBF package
library(SBF)
```
- install directly from Github via `devtools`
```
library(devtools)
devtools::install_github("amalthomas111/SBF")
# load package
library(SBF)
```

Please contact us via e-mail or through
a [GitHub issue](https://github.com/amalthomas111/SBF/issues)
if there is any trouble with installation.

### Examples

#### load test dataset

- SBF package has a sample gene expression data with the mean expression of
nine tissues in five species

```
# load sample dataset from SBF package
avg_counts <- SBF::avg_counts
# dimension of matrices for different species
sapply(avg_counts)
# names of different tissue types in humans
names(avg_counts[["Homo_sapiens"]])
```

#### SBF computation

```
# compute SBF factorization
avg_counts <- SBF::avg_counts
sbf <- SBF(matrix_list = avg_counts, col_index = 2, approximate = FALSE,
           transform_matrix = FALSE)
```
`?SBF` help function shows all arguments for the SBF call.

#### SBF computation based on inter-sample correlation

```
# SBF call using correlation matrix
avg_counts <- SBF::avg_counts
sbf_cor <- SBF(matrix_list = avg_counts, col_index = 2, approximate = FALSE,
               transform_matrix = TRUE)
```

#### Approximate SBF (A-SBF) computation

```
# A-SBF call
avg_counts <- SBF::avg_counts
asbf <- SBF(matrix_list = avg_counts, col_index = 2, approximate = TRUE,
            transform_matrix = FALSE)
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, asbf$delta, asbf$u, asbf$v)
````

#### Approximate SBF (A-SBF) computation based on inter-sample correlation

````
# A-SBF call using correlation matrix
avg_counts <- SBF::avg_counts
asbf_cor <- SBF(matrix_list = avg_counts, col_index = 2, approximate = TRUE,
                transform_matrix = TRUE)
# calculate decomposition error
decomperror <- calcDecompError(avg_counts, asbf_cor$delta, asbf_cor$u,
                                asbf_cor$v)
````

### Contacts ###

***Amal Thomas*** *amalthom@usc.edu*

***Andrew D Smith*** *andrewds@usc.edu*


### Copyright ###

Copyright (C) 2018-2021  Amal Thomas and Andrew D. Smith

Authors: Amal Thomas and Andrew D. Smith

SBF is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

SBF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.
