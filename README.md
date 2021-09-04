## SBF: A R package for **S**hared **B**asis **F**actorization

**S**hared **B**asis **F**actorization (SBF) is a joint matrix diagonalization
approach we developed for cross-species gene expression analysis.
Approximate Shared Basis Factorization (A-SBF) is an extention of the SBF
approach.

### Installation
- Clone from Github
```
git clone https://github.com/amalthomas111/SBF.git
# Inside R
library(devtools)
install("<path to SBF>/SBF")
# load package
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
data(avg_counts, package = "SBF")
# dimension of matrices for different species
sapply(avg_counts)
# names of different tissue types in humans
names(avg_counts[["Homo_sapiens"]])
```

#### SBF computation

```
# compute SBF factorization
data(avg_counts, package = "SBF")
sbf = SBF(avg_counts = avg_counts, colIndex = 2, approximate = FALSE,
              transformData = FALSE)
```
`?SBF` help function shows all arguments for the SBF call.

#### SBF computation based on inter-sample correlation

```
# SBF call using correlation matrix
sbf.cor = SBF(avg_counts = avg_counts,colIndex = 2, approximate = FALSE,
               transformData = TRUE)
```

#### Approximate SBF (A-SBF) computation

```
# A-SBF call
asbf = SBF(avg_counts = avg_counts,colIndex = 2, approximate = TRUE,
        transformData = FALSE)
````

#### Approximate SBF (A-SBF) computation based on inter-sample correlation

````
# A-SBF call using correlation matrix
asbf.cor = SBF(avg_counts = avg_counts,colIndex = 2, approximate = TRUE,
            transformData = TRUE)
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
