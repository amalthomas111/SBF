## SBF: A R package for **S**hared **B**asis **F**actorization 

### Installation
```
git clone https://github.com/amalthomas111/SBF.git
# Inside R
install.packages("<path to SBF>/SBF",repos=NULL)
library(SBF)
```

If you have trouble with the installation, please contact us via e-mail or
through a [GitHub issue](https://github.com/amalthomas111/issues).

### Examples

#### load test dataset

- Load a sample gene expession matrix with mean expression values for different
tissues from different species

```
# check sample dataset in from SBF package
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
`?SBF` shows all arguments for the SBF function

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

#### Approximate SBF (A-SBF) computation based on inter-sample

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

abismal is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.
