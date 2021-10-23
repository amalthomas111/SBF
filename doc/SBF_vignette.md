<style>
body {
text-align: justify}
</style>
Background
==========

Joint matrix factorization facilitates the comparison of expression
profiles from different species without using gene mapping. Transforming
gene expression profiles into reduced eigengene space using singular
value decomposition (SVD) has been shown to capture meaningful
biological information (Alter, Brown, and Botstein 2000). Tamayo et al.
(2007) used a non-negative matrix factorization approach to learn a
low-dimensional approximation of the microarray expression datasets and
used the reduced space for comparisons. Matrix factorization-based
methods are commonly used for gene expression analysis (Alter, Brown,
and Botstein 2000; Tamayo et al. 2007). An orthology independent matrix
factorization framework based on generalized singular value
decomposition (GSVD; Van Loan 1976) was used by Alter, Brown, and
Botstein (2003) to compare gene-expression profiles from two species.
This framework was later extended to develop higher-order generalized
singular value decomposition (HO GSVD) to analyze data from more than
two species (Ponnapalli et al. 2011).

This study developed a joint diagonalization approach called approximate
shared basis factorization (A-SBF) for cross-species expression
comparisons. This approach extends the exact factorization approach we
developed called shared basis factorization (SBF). We discuss the
details of the two methods in the following sections.

Shared basis factorization
==========================

Consider a set of real matrices
*D*<sub>*i*</sub> ∈ ℝ<sup>*m*<sub>*i*</sub> × *k*</sup>
(*i* = 1, …, *N*) with full column rank. We define shared basis
factorization (SBF) as

Here each *U*<sub>*i*</sub> ∈ ℝ<sup>*m*<sub>*i*</sub> × *k*</sup> is a
dataset-specific left basis matrix, each
*Δ*<sub>*i*</sub> ∈ ℝ<sup>*k* × *k*</sup> is a diagonal matrix with
positive values *δ*<sub>*i**k*</sub>, and *V* is a square invertible
matrix.

Estimating the shared right basis matrix
----------------------------------------

Let *M* be the scaled sum of the
*D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub>. We define *M* is
defined as
$$
 M = \\frac{\\sum\_{i=1}^{N} D\_i^T D\_i/w\_i}{\\alpha}.
$$

The scaling factor *w*<sub>*i*</sub> is the total variance explained by
the column vectors of *D*<sub>*i*</sub>, and *α* is the inverse sum of
the total variance of *D*<sub>*i*</sub>, for *i* = 1⋯*N*. The weights
*w*<sub>*i*</sub> and *α* are defined as

Here
$\\sum\_{j=1}^{k} \\sigma\_{jj}^{2\\mbox{ }(i)} = \\mathrm{tr}(D\_i^T D\_i) = \\mathrm{tr}(A\_i)$.
Using the *w*<sub>*i*</sub> and *α*, individual
*D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub> are standardized. If
all the variances are equal, *M* becomes the arithmetic mean of the sum
of *D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub>. The shared right
basis matrix *V* is then determined from the eigenvalue decomposition of
*M*, where *M* = *V**Θ**V*<sup>*T*</sup>. The shared right basis matrix
*V* is an orthogonal matrix as *M* is symmetric. Given *V*, we compute
*U*<sub>*i*</sub> and *Δ*<sub>*i*</sub> by solving the linear system
*D*<sub>*i*</sub>*V* = *U*<sub>*i*</sub>*Δ*<sub>*i*</sub> = *L*<sub>*i*</sub>.
By normalizing the columns of *L*<sub>*i*</sub>, we have
*δ*<sub>*i**k*</sub> = ∥*l*<sub>*i**k*</sub>∥ and
*Δ*<sub>*i*</sub> = diag(*δ*<sub>*i*1</sub>, …, *δ*<sub>*i**k*</sub>).

Approximate shared basis factorization
======================================

Consider a set of matrices
*D*<sub>*i*</sub> ∈ ℝ<sup>*m*<sub>*i*</sub> × *k*</sup>
(*i* = 1, …, *N*), each with full column rank. We define approximate
shared basis factorization (A-SBF) as

Each *U*<sub>*i*</sub> ∈ ℝ<sup>*m*<sub>*i*</sub> × *k*</sup> is a
species-specific left basis matrix with **orthonormal** columns
(eigengenes), *Δ*<sub>*i*</sub> ∈ ℝ<sup>*k* × *k*</sup> is a diagonal
matrix with positive values *Δ*<sub>*i**k*</sub> and *V* is a non
singular square matrix. The right basis matrix *V* is identical in all
the *N* matrix factorizations and defines the common space shared by all
species. We estimate the factorization such that the estimated *V* is
closest to that in the exact decomposition and by minimizing the total
decomposition error
$\\sum^{N}\_{i=1}\\epsilon\_i = {\\sum^{N}\_{i=1}\\|D\_i - U\_i\\Delta\_iV^T\\|^2}\_F$.

Usage cases
===========

SBF examples
------------

    # load SBF package
    library(SBF)

Let us create some random matrices using the `createRandomMatrices`
function from the SBF package. We will create four matrices, each with
three columns with rows varying from 4 to 6.

    set.seed(1231)
    mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)
    sapply(mymat, dim)
    #>      mat1 mat2 mat3 mat4
    #> [1,]    5    6    4    5
    #> [2,]    3    3    3    3

Rank of each of this matrices

    sapply(mymat, function(x) {
      qr(x)$rank
      })
    #> mat1 mat2 mat3 mat4 
    #>    3    3    3    3

Let us compute SBF using different approaches.

-   Estimate *V* using sum of
    *D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub>/*N*
-   Estimate *V* using sum of
    *D*<sub>*i*</sub><sup>*T*</sup>*D*<sub>*i*</sub>/*N* with inverse
    variance weighting
-   Estimate *V* using inter-sample correlation

<!-- -->

    sbf <- SBF(matrix_list = mymat)
    sbf_inv <- SBF(matrix_list = mymat, weighted = TRUE)
    sbf_cor <- SBF(matrix_list = mymat, transform_matrix = TRUE)

When *D*<sub>*i*</sub>’s are transformed to compute inter-sample
correlation, we do not need to scale it using inverse-variance weighting
anymore. We recommend using inverse variance weights, giving a more
robust estimate of *V* when noisy datasets are present. We estimate *V*
using inter-sample correlation when dealing with gene expression data
sets.

`?SBF` help function shows all arguments for the SBF function. Let us
inspect the output of the `SBF` call.

    names(sbf)
    #> [1] "v"      "lambda" "u"      "delta"  "m"

`sbf$u`, `sbf$v`, and `sbf$delta` correspond to the estimated left basis
matrix, shared right basis matrix, and diagonal matrices.

The estimated *V* ∈ *R*<sup>*k* × *k*</sup> has a dimension of
*k* × *k*, where *k* is the number of columns in *D*<sub>*i*</sub>.

    sbf$v
    #>           [,1]       [,2]       [,3]
    #> [1,] 0.4793022  0.8669998  0.1363110
    #> [2,] 0.7027353 -0.2860780 -0.6514004
    #> [3,] 0.5257684 -0.4080082  0.7463892

The delta values for each matrix for the three cases are shown below.

    printDelta <- function(l) {
      for (eachmat in names(l$delta)) {
      cat(eachmat, ":", l$delta[[eachmat]], "\n")
      }
    }
    cat("sbf\n");printDelta(sbf)
    #> sbf
    #> mat1 : 205.4915 29.6746 71.43295 
    #> mat2 : 206.5816 71.72548 55.682 
    #> mat3 : 189.9136 52.6758 42.36825 
    #> mat4 : 192.6911 80.22868 58.57913
    cat("sbf_inv\n");printDelta(sbf_inv)
    #> sbf_inv
    #> mat1 : 205.5109 22.4888 73.95623 
    #> mat2 : 206.5963 77.00352 48.05638 
    #> mat3 : 189.8719 58.60394 33.92988 
    #> mat4 : 192.6942 72.41955 67.98802
    cat("sbf_cor\n");printDelta(sbf_cor)
    #> sbf_cor
    #> mat1 : 200.4197 44.25134 77.99852 
    #> mat2 : 199.8621 80.34494 67.23723 
    #> mat3 : 176.1738 76.95449 60.64485 
    #> mat4 : 185.4579 67.51833 89.69184

The *V* ∈ *R*<sup>*k* × *k*</sup> estimated in SBF is orthogonal. So
*V*<sup>*T*</sup>*V* = *V**V*<sup>*T*</sup> = *I*.

    zapsmall(t(sbf$v) %*% sbf$v)
    #>      [,1] [,2] [,3]
    #> [1,]    1    0    0
    #> [2,]    0    1    0
    #> [3,]    0    0    1

The estimated *V* is an invertible matrix.

    qr(sbf$v)$rank
    #> [1] 3

The *U*<sub>*i*</sub> matrices estimated in the SBF do not have
orthonormal columns. Let us explore that.

    sapply(sbf$u, dim)
    #>      mat1 mat2 mat3 mat4
    #> [1,]    5    6    4    5
    #> [2,]    3    3    3    3

Let us take the first matrix
*U*<sub>*i*</sub> ∈ *R*<sup>*m*<sub>*i*</sub> × *k*</sup> to check this.
For this matrix, *U*<sub>*i*</sub><sup>*T*</sup>*U*<sub>*i*</sub> will
be *k* × *k* matrix where *k* = 3.

    t(sbf$u[[names(sbf$u)[1]]]) %*% sbf$u[[names(sbf$u)[1]]]
    #>             [,1]        [,2]       [,3]
    #> [1,]  1.00000000 -0.07071457  0.1405487
    #> [2,] -0.07071457  1.00000000 -0.6201468
    #> [3,]  0.14054867 -0.62014676  1.0000000

The estimated *M* matrix is stored `sbf$m` and `sbf$lambda` gives the
eigenvalues in the eigenvalue decomposition
(*M* = *V**Θ**V*<sup>*T*</sup>).

    sbf$lambda
    #> [1] 39524.940  3809.127  3357.434

SBF is an exact factorization. Let compute the factorization error for
the three cases using `calcDecompError` function.

    calcDecompError(mymat, sbf$u, sbf$delta, sbf$v)
    #> [1] 1.851693e-26
    calcDecompError(mymat, sbf_inv$u, sbf_inv$delta, sbf_inv$v)
    #> [1] 1.894292e-26
    calcDecompError(mymat, sbf_cor$u, sbf_cor$delta, sbf_cor$v)
    #> [1] 2.835482e-26

The errors are close to zero in all three cases.

### Adding new dataset

The total column variance of matrix 1-4 in `mymat` is approximately in
the same range.

    sapply(mymat, function(x) sum(diag(cov(x))))
    #>    mat1    mat2    mat3    mat4 
    #> 2076.80 2273.50 2375.25 2860.40

Now, let us create two new matrix lists containing the `mymat`. We will
add a dataset with a similar variance to the first list and a high
variance to the second.

    mat5 <- matrix(c(130, 183, 62, 97, 147, 94, 102, 192, 19), byrow = T,
                        nrow = 3, ncol = 3)
    mat5_highvar <- matrix(c(406, 319, 388, 292, 473, 287, 390, 533, 452),
                           byrow = T, nrow = 3, ncol = 3)

    mymat_new <- mymat
    mymat_new[["mat5"]] <- mat5
    sapply(mymat_new, function(x) sum(diag(cov(x))))
    #>     mat1     mat2     mat3     mat4     mat5 
    #> 2076.800 2273.500 2375.250 2860.400 2299.667
    mymat_new_noisy <- mymat
    mymat_new_noisy[["mat5"]] <- mat5_highvar
    sapply(mymat_new_noisy, function(x) sum(diag(cov(x))))
    #>     mat1     mat2     mat3     mat4     mat5 
    #>  2076.80  2273.50  2375.25  2860.40 22915.00

Let us compute SBF with the new datasets.

    sbf_new <- SBF(matrix_list = mymat_new)
    sbf_inv_new <- SBF(matrix_list = mymat_new, weighted = TRUE)


    sbf_new_noisy <- SBF(matrix_list = mymat_new_noisy)
    sbf_inv_new_noisy <- SBF(matrix_list = mymat_new_noisy, weighted = TRUE)

Let us take the newly estimated values *U*<sub>*i*</sub>,
*Δ*<sub>*i*</sub>, and *V* for the four initial matrices in `mymat`. We
will then compare the decomposition error for the two cases with and
without inverse variance weighting.

    e1 <- calcDecompError(mymat, sbf_new$u[1:4], sbf_new$delta[1:4], sbf_new$v)
    e2 <- calcDecompError(mymat, sbf_new_noisy$u[1:4], sbf_new_noisy$delta[1:4],
                          sbf_new_noisy$v)
    e2 / e1
    #> [1] 3.887268

    e3 <- calcDecompError(mymat, sbf_inv_new$u[1:4], sbf_inv_new$delta[1:4],
                          sbf_inv_new$v)
    e4 <- calcDecompError(mymat, sbf_inv_new_noisy$u[1:4],
                          sbf_inv_new_noisy$delta[1:4], sbf_inv_new_noisy$v)
    e4 / e3
    #> [1] 1.657059

With inverse variance weighting, the deviation is smaller.

A-SBF examples
--------------

    set.seed(1231)
    mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)
    sapply(mymat, dim)
    #>      mat1 mat2 mat3 mat4
    #> [1,]    5    6    4    5
    #> [2,]    3    3    3    3

Now let us compute Approximate-SBF for the same datasets.

-   A-SBF
-   A-SBF with inverse variance weighting
-   A-SBF with inter-sample correlation

<!-- -->

    asbf <- SBF(matrix_list = mymat, approximate = TRUE)
    asbf_inv <- SBF(matrix_list = mymat, weighted = TRUE, approximate = TRUE)
    asbf_cor <- SBF(matrix_list = mymat, approximate = TRUE, transform_matrix = TRUE)

    names(asbf)
    #> [1] "v"       "lambda"  "u"       "u_ortho" "delta"   "m"       "error"

A-SBF is not an exact factorization. A-SBF output has two additional
values. `asbf$u_ortho` is the left basis matrix with orthonormal columns
and `asbf$error` gives the decomposition error.

    asbf$error
    #> [1] 2329.73
    asbf_inv$error
    #> [1] 1651.901
    asbf_cor$error
    #> [1] 14045.99

The same error can also be computed using the `calcDecompError`
function.

    calcDecompError(mymat, asbf$u_ortho, asbf$delta, asbf$v)
    #> [1] 2329.73
    calcDecompError(mymat, asbf_inv$u_ortho, asbf_inv$delta, asbf_inv$v)
    #> [1] 1651.901
    calcDecompError(mymat, asbf_cor$u_ortho, asbf_cor$delta, asbf_cor$v)
    #> [1] 14045.99

In A-SBF factorization, *U*<sub>*i*</sub> has orthonormal columns, and
*V* is orthogonal.

    zapsmall(t(asbf$u_ortho[[names(asbf$u_ortho)[1]]]) %*%
               asbf$u_ortho[[names(asbf$u_ortho)[1]]])
    #>      [,1] [,2] [,3]
    #> [1,]    1    0    0
    #> [2,]    0    1    0
    #> [3,]    0    0    1

    zapsmall(t(asbf$v) %*% asbf$v)
    #>      [,1] [,2] [,3]
    #> [1,]    1    0    0
    #> [2,]    0    1    0
    #> [3,]    0    0    1

Minimizing A-SBF factorization error
====================================

Here we provide details on minimizing the decomposition error of the
A-SBF method. In A-SBF, we want to solve the following optimization
problem:
$$
\\begin{array}{rl}
  \\mbox{minimize:} & \\sum^{n}\_{i=1} \\|D\_i - U\_i\\Delta\_iV^T\\|^2\_F,\\\\\[0.5em\]
  \\mbox{subject to:} & U\_i^TU\_i = I, \\mbox{ for } 1\\leq i \\leq n,\\\\
  & V^TV = VV^T = I,\\\\
  & \\Delta\_i = \\mathrm{diag}(\\delta\_{i,1},\\ldots,\\delta\_{i,k}),\\\\
  & \\delta\_{i,j} \\geq 0, \\mbox{ for } 1\\leq i\\leq n.
  \\end{array}
$$

The objective function can be reformulated as
$$
\\textstyle
\\sum\_{i=1}^n
\\mathrm{tr}(D\_i^TD\_i)
- 2\\mathrm{tr}\\big(D\_i^TU\_i\\Delta\_iV^T\\big)
+ \\mathrm{tr}(\\Delta\_i^2)
$$
Let *Φ*, *Ψ*<sub>*i*</sub>, and *Θ*<sub>*i*</sub> be matrices with
Lagrange multipliers for constraints *V*<sup>*T*</sup>*V* − *I* = 0,
*U*<sub>*i*</sub><sup>*T*</sup>*U*<sub>*i*</sub> − *I* = 0 and
*δ*<sub>*i*, *j*</sub> ≥ 0. The Lagrange ℒ is

$$
\\mathcal{L}(U, \\Delta, V) = \\sum\_{i=1}^n \\mathrm{tr}(D\_i^T D\_i - 2 D\_i^T U\_i \\Delta\_i V^T + \\Delta\_i^2) +
\\mathrm{tr}(\\Phi(V^T V - I)) + \\mathrm{tr}(\\Psi\_i(U\_i^T U\_i - I)) + \\mathrm{tr}(\\Theta\_i \\Delta\_i).
$$

Solving for an individual *U*<sub>*i*</sub> using the partial derivative
of ℒ with respect to *U*<sub>*i*</sub> gives

$$
\\begin{aligned}
  \\frac{\\partial}{\\partial U\_i} \\mathcal{L}(U, \\Delta, V) &= -2D\_iV\\Delta\_i + U\_i(\\Psi\_i^T + \\Psi\_i)\\\\
  U\_i &= D\_iV\\Delta\_i((D\_iV\\Delta\_i)^T(D\_iV\\Delta\_i))^{-1/2}. \\label{eqU}
\\end{aligned}
$$

Solving *V* using the partial derivative of ℒ with respect to *V* gives

$$
\\begin{aligned}
  \\frac{\\partial}{\\partial V} \\mathcal{L}(U, \\Delta, V)
  &= \\sum\_{i=1}^n -2D\_i^TU\_i\\Delta\_i + V(\\Phi^T + \\Phi) \\\\
  V &= \\sum\_{i=1}^n D\_i^TU\_i\\Delta\_i ( (\\sum\_{i=1}^n D\_i^TU\_i\\Delta\_i)^T (\\sum\_{i=1}^n D\_i^TU\_i\\Delta\_i))^{-1/2}. \\label{eqV}
\\end{aligned}
$$

Solving for an individual *Δ*<sub>*i*</sub> using the partial derivative
of ℒ with respect to *Δ*<sub>*i*</sub> gives

$$
\\begin{aligned}
\\frac{\\partial}{\\partial \\Delta\_i} \\mathcal{L}(U, \\Delta, V) &= -2 U\_i^T D\_iV + 2\\Delta\_i  + \\Theta\_i \\\\
  \\Delta\_i &= U\_i^T D\_i V. \\label{eqDelta}
\\end{aligned}
$$

In practice, we do not compute the square root inverse to update
*U*<sub>*i*</sub> and *V*. Now we show that we make the optimal update
of our matrix in each iteration based on the singular value
decomposition (SVD) of other matrices. The first two results below show
that, given one matrix factor, we can obtain the other as required in
our algorithm.

Proposition I
-------------

Consider matrices *A* ∈ ℝ<sup>*m* × *n*</sup> with rank *n* and
*Q* ∈ ℝ<sup>*n* × *n*</sup>. If *X**Σ**Y*<sup>*T*</sup> is the SVD of
*A**Q*<sup>*T*</sup>, then among matrices *B* ∈ ℝ<sup>*m* × *n*</sup>
with orthonormal columns, ∥*A* − *B**Q*∥<sub>*F*</sub> is minimized when
*B* = *X**Y*<sup>*T*</sup>.

The norm and its square have the same minimum value, so we consider the
minimum of
∥*A* − *B**Q*∥<sub>*F*</sub><sup>2</sup> = *t**r*(*A*<sup>*T*</sup>*A*) + *t**r*(*Q*<sup>*T*</sup>*Q*) − 2*t**r*(*A**Q*<sup>*T*</sup>*B*<sup>*T*</sup>),

which coincides with the maximum value of
*t**r*(*A**Q*<sup>*T*</sup>*B*<sup>*T*</sup>). With
*X**Σ**Y*<sup>*T*</sup> the SVD of *A**Q*<sup>*T*</sup>, define
*Z* = *Y*<sup>*T*</sup>*B*<sup>*T*</sup>*X*. Since
*B*<sup>*T*</sup>*B* = *I*, we have *Z**Z*<sup>*T*</sup> = *I*. We can
bound the maximum value of *t**r*(*A**Q*<sup>*T*</sup>*B*<sup>*T*</sup>)
as follows:

$$
\\mathrm{tr}(AQ^TB^T) = \\mathrm{tr}(X\\Sigma Y^TB^T)
    = \\mathrm{tr}(Z\\Sigma) =
    \\sum\_{i=1}^{n} Z\_{ii}\\sigma\_i \\leq \\sum\_{i=1}^{n} \\sigma\_i,
$$
and this bound is attained when *Z* = *I*, and thus
*B* = *X**Y*<sup>*T*</sup>.

In proposition we substitute
*A* = *D*<sub>*i*</sub>, *B* = *U*<sub>*i*</sub>, and
*Q* = *Δ*<sub>*i*</sub>*V*<sup>*T*</sup>. Now we have
*A**Q*<sup>*T*</sup> = *D*<sub>*i*</sub>*V**Δ*<sub>*i*</sub> = *X**Σ**Y*<sup>*T*</sup>.
The *U*<sub>*i*</sub> with orthonormal columns is given by
*X**Y*<sup>*T*</sup>. The same can also be obtained by substituting
*D*<sub>*i*</sub>*V**Δ*<sub>*i*</sub> = *X**Σ**Y*<sup>*T*</sup> in
equation 1.

$$
\\begin{aligned}
    U\_i &= D\_iV\\Delta\_i((D\_iV\\Delta\_i)^T(D\_iV\\Delta\_i))^{-1/2} \\\\
    U\_i &= X \\Sigma Y^T((X \\Sigma Y^T)^T(X \\Sigma Y^T))^{-1/2} = X Y^T.
\\end{aligned}
$$

Proposition II
--------------

Consider matrices *A* ∈ ℝ<sup>*m* × *n*</sup> with rank *n* and
*B* ∈ ℝ<sup>*m* × *n*</sup>. If *X**Σ**Y*<sup>*T*</sup> is the SVD of
*A*<sup>*T*</sup>*B*, then among orthogonal matrices
*Q* ∈ ℝ<sup>*n* × *n*</sup>, ∥*A* − *B**Q*<sup>*T*</sup>∥<sub>*F*</sub>
is minimized when *Q* = *X**Y*<sup>*T*</sup>.

The norm and its square have the same minimum value, so we consider the
minimum of
∥*A* − *B**Q*<sup>*T*</sup>∥<sub>*F*</sub><sup>2</sup> = *t**r*(*A*<sup>*T*</sup>*A*) + *t**r*(*B*<sup>*T*</sup>*B*) − 2*t**r*(*A*<sup>*T*</sup>*B**Q*<sup>*T*</sup>),
which coincides with the maximum value of
*t**r*(*A*<sup>*T*</sup>*B**Q*<sup>*T*</sup>). With
*X**Σ**Y*<sup>*T*</sup> the SVD of *A*<sup>*T*</sup>*B*, define
orthogonal matrix *Z* = *Y*<sup>*T*</sup>*Q*<sup>*T*</sup>*X*. We bound
the maximum of *t**r*(*A*<sup>*T*</sup>*B**Q*<sup>*T*</sup>) as

$$
\\mathrm{tr}(A^TBQ^T) = \\mathrm{tr}(X\\Sigma Y^TQ^T)
    = \\mathrm{tr}(Z\\Sigma)
    = \\sum\_{i=1}^{n} Z\_{ii}\\sigma\_i \\leq \\sum\_{i=1}^{n} \\sigma\_i,
$$
and this bound is attained when *Z* = *I* and thus
*Q* = *X**Y*<sup>*T*</sup>.

In proposition II, we substitute
*A* = *D*<sub>*i*</sub>, *B* = *U*<sub>*i*</sub>*Δ*<sub>*i*</sub>, and
*Q* = *V*<sub>*i*</sub>. Now we have
*A*<sup>*T*</sup>*B* = *D*<sub>*i*</sub><sup>*T*</sup>*U*<sub>*i*</sub>*Δ*<sub>*i*</sub> = *X**Σ**Y*<sup>*T*</sup>
and orthogonal *V*<sub>*i*</sub> is given by *X**Y*<sup>*T*</sup>. In
A-SBF, *V* is shared across the *n* factorizations. To find the shared
*V*, we can use the following proposition.

Proposition III
---------------

For a set of *n* matrices
*A*<sub>1</sub>, *A*<sub>2</sub>, …, *A*<sub>*n*</sub>, where
*A*<sub>*i*</sub> ∈ ℝ<sup>*k* × *k*</sup>, the orthogonal matrix
*Q* ∈ ℝ<sup>*k* × *k*</sup> minimizing
$\\sum\_{i=1}^n\\|A\_i- Q\\|^2\_F$ is given by
*Q* = *X**Y*<sup>*T*</sup>, where SVD of
$\\sum\_{i=1}^nA\_i=X\\Sigma Y^T$.

We seek to minimize
$$
 \\sum\_{i=1}^n\\|A\_i - Q\\|^2\_F ~~\\mbox{subject to}~~ Q^TQ=I.
$$

Defining $A = \\sum\_{i=1}^n A\_i$, this objective function can be
rewritten

$$
\\mathrm{tr}(\\sum\_{i=1}^nA\_i^TA\_i) - 2\\mathrm{tr}(A^TQ) + n\\mathrm{tr}(Q^TQ).
$$

Let *Φ* be a matrix of Lagrange multipliers for
*Q*<sup>*T*</sup>*Q* − *I* = 0. The Lagrangian ℒ is then
$$
  \\mathcal{L}(Q,\\Phi) = \\mathrm{tr}(\\sum\_{i=1}^nA\_i^TA\_i) - 2\\mathrm{tr}(A^TQ) + \\mathrm{tr}(Q^TQ) + \\mathrm{tr}(\\Phi(Q^TQ-I)).
$$

Since *Q* is orthogonal, the partial derivatives of ℒ(*Q*, *Φ*) with
respect to *Q* are as follows:
$$
\\frac{\\partial \\mathcal{L}(Q,\\Phi)}{\\partial Q} = -2A + Q(\\Phi^T + \\Phi).
$$

Setting this equation to 0 yields:

$$
  Q = A\\left(\\frac{\\Phi^T + \\Phi}{2}\\right)^{-1} ~~\\mbox{and}~~ A = Q\\left(\\frac{\\Phi^T + \\Phi}{2}\\right).
$$

Again because *Q* is orthogonal,

$$
A^TA = \\left(\\frac{\\Phi^T + \\Phi}{2} \\right)^2 ~~\\mbox{and}~~
    \\left(\\frac{\\Phi^T + \\Phi}{2}\\right) = (A^TA)^{1/2},
$$
which implies *Q* = *A*(*A*<sup>*T*</sup>*A*)<sup> − 1/2</sup>. Since
*X**Σ**Y*<sup>*T*</sup> is the SVD of *A* we conclude

*Q* = *X**Σ**Y*<sup>*T*</sup>(*Y**Σ*<sup>2</sup>*Y*<sup>*T*</sup>)<sup> − 1/2</sup> = *X**Y*<sup>*T*</sup>.
In proposition , we substitute
*A*<sub>*i*</sub> = *V*<sub>*i*</sub> = *D*<sub>*i*</sub><sup>*T*</sup>*U*<sub>*i*</sub>*Δ*<sub>*i*</sub>
and *Q* = *V*. Now we have
$\\sum\_{i=1}^n A\_i = \\sum\_{i=1}^n D\_i^T U\_i \\Delta\_i = X \\Sigma Y^T$
and the shared basis in A-SBF is given by *V* = *X**Y*<sup>*T*</sup>.

Iterative update
----------------

Using these three propositions, we iteratively update *U*<sub>*i*</sub>,
*Δ*<sub>*i*</sub>, and *V*. The steps of the iterative algorithm are
shown below.

1.  Input: *D*<sub>*i*</sub>, *i* = 1…*N*
2.  Output: *U*<sub>*i*</sub>, *Δ*<sub>*i*</sub>, and *V*, where
    *U*<sub>*i*</sub><sup>*T*</sup>*U*<sub>*i*</sub> = *I* and
    *V*<sup>*T*</sup>*V* = *V**V*<sup>*T*</sup> = *I*
3.  Initialize *U*<sub>*i*</sub><sup>*k*</sup>,
    *Δ*<sub>*i*</sub><sup>*k*</sup>, and *V*<sup>*k*</sup>. Compute
    $\\epsilon^{k,k,k} = {\\sum^{N}\_{i=1}\\|D\_i - U\_i^{k} \\Delta\_i^{k} {V^{k}}^{\\!T}\\|^{2}}\_F$
4.  Update **U**<sub>**i**</sub><sup>**k** **+** **1**</sup>:  
       *U*<sub>*i*</sub><sup>*k* + 1</sup> = *Z**Y*<sup>*T*</sup>, where
    SVD of
    *D*<sub>*i*</sub>*V*<sup>*k*</sup>*Δ*<sub>*i*</sub><sup>*k*</sup> = *Z**Σ**Y*<sup>*T*</sup>.
    Compute *ϵ*<sup>*k* + 1, *k*, *k*</sup>  
       If
    *ϵ*<sup>*k* + 1, *k*, *k*</sup> &lt; *ϵ*<sup>*k*, *k*, *k*</sup>:  
          *U*<sub>*i*</sub> ← *U*<sub>*i*</sub><sup>*k* + 1</sup>.
5.  Update **Δ**<sub>**i**</sub><sup>**k** **+** **1**</sup>:  
       **Δ**<sub>**i**</sub><sup>**k** **+** **1**</sup> = diag((*U*<sub>*i*</sub><sup>*k* + 1</sup>)<sup>*T*</sup>*D*<sub>*i*</sub>*V*<sup>*k*</sup>).
    Compute *ϵ*<sup>*k* + 1, *k* + 1, *k*</sup>  
       If
    *ϵ*<sup>*k* + 1, *k* + 1, *k*</sup> &lt; *ϵ*<sup>*k* + 1, *k*, *k*</sup>:  
          *Δ*<sub>*i*</sub> ← *Δ*<sub>*i*</sub><sup>*k* + 1</sup>
6.  Update **V**<sup>**k** **+** **1**</sup>:  
       *V*<sup>*k* + 1</sup> = *M**Q*<sup>*T*</sup>, where SVD of
    $\\sum\_{i}^N D\_i^T U\_i^{k+1} \\Delta\_i^{k+1} = M \\Phi Q^T$.
    Compute *ϵ*<sup>*k* + 1, *k* + 1, *k* + 1</sup>  
       If
    *ϵ*<sup>*k* + 1, *k* + 1, *k* + 1</sup> &lt; *ϵ*<sup>*k* + 1, *k* + 1, *k*</sup>:  
          *V* ← *V*<sup>*k* + 1</sup>
7.  Repeat steps 4-6 until convergence.

Our iterative approach is a block-coordinate descent algorithm. In
gradient descent, we have the general update:
*θ*<sub>*t* + 1</sub> = *θ*<sub>*t*</sub> − *η*∇ℱ(*θ*<sub>*t*</sub>),
where *t* is the iteration counter, *η* is the learning rate and ∇ℱ(*θ*)
is the gradient of the cost function. In A-SBF, the gradient of the cost
function with respect to *U* is  − 2*D**V**Δ*. So we have

*U*<sub>*t* + 1</sub> = *U*<sub>*t*</sub> + *η**D**V**Δ*.
By setting *η* = *I* − *U*<sub>*t*</sub>(*D**V**Δ*)<sup> − 1</sup>, we
have *U*<sub>*t* + 1</sub> = *D**V**Δ*. Since we require orthonormal
columns, we find the closest orthogonal matrix to *D**V**Δ*, and set the
final value of *U*<sub>*t* + 1</sub> = *Z**Y*<sup>*T*</sup>, where SVD
of *D**V**Δ* = *Z**Σ**Y*<sup>*T*</sup>. This approach employs orthogonal
Procrustes solution along with the gradient descent algorithm. We can
also achieve this by directly setting
$\\eta = \\frac{Z Y^T - U\_t}{D V \\Delta}$, so that
*U*<sub>*t* + 1</sub> has all the properties we need. The learning rate
determines the size of the steps. If the value of *U*<sub>*t*</sub> is
very different from *Z**Y*<sup>*T*</sup>, *η* will be high, and we take
larger steps. As it becomes closer, the learning rate also decreases. We
have a similar case for updating *V* and *Δ*.

Examples
--------

### Optimizing A-SBF error

Let us optimize the factorization error using the
`optimizeFactorization` function for the three cases of A-SBF
computation.

    set.seed(1231)
    mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)
    asbf <- SBF(matrix_list = mymat, approximate = TRUE)
    asbf_inv <- SBF(matrix_list = mymat, weighted = TRUE, approximate = TRUE)
    asbf_cor <- SBF(matrix_list = mymat, approximate = TRUE, transform_matrix = TRUE)

Depending upon the data matrix and initial values of *U*<sub>*i*</sub>,
*Δ*<sub>*i*</sub>, and *V*, optimization could take some time.

    myopt <- optimizeFactorization(mymat, asbf$u_ortho, asbf$delta, asbf$v)
    myopt_inv <- optimizeFactorization(mymat, asbf_inv$u_ortho, asbf_inv$delta,
                                       asbf_inv$v)
    myopt_cor <- optimizeFactorization(mymat, asbf_cor$u_ortho, asbf_cor$delta,
                                       asbf_cor$v)

The number of iteration taken for optimizing and new factorization
error:

    cat("For asbf, # iteration =", myopt$error_pos, "final error =", myopt$error)
    #> For asbf, # iteration = 220 final error = 1411.555
    cat("\nFor asbf inv, # iteration =", myopt_inv$error_pos, "final error =",
        myopt_inv$error)
    #> 
    #> For asbf inv, # iteration = 202 final error = 1411.555
    cat("\nFor asbf cor, # iteration =", myopt_cor$error_pos, "final error =",
        myopt_cor$error)
    #> 
    #> For asbf cor, # iteration = 196 final error = 1411.555

    #> After optimization, for all three A-SBF factorizations, the final error
    #> is the same (up to 2 decimals).
    #> The final error is 1411.56

### Using different initial values

    set.seed(1231)
    mymat <- createRandomMatrices(n = 4, ncols = 3, nrows = 4:6)

1.  Let us initialize the `optimizeFactorization` function with a random
    orthogonal matrix and check the final optimization error. The *V*
    matrix estimated from the `mymat` matrix has a dimension of 3 × 3.
    First we will create a random 3 × 3 matrix and obtain an orthogonal
    matrix based on this.

<!-- -->

    set.seed(111)
    rand_mat <- createRandomMatrices(n = 1, ncols = 3, nrows = 3)
    cat("\nRank is:", qr(rand_mat[[1]])$rank,"\n")
    #> 
    #> Rank is: 3
    dim(rand_mat[[1]])
    #> [1] 3 3

Get an orthogonal *V* matrix using SVD. We will set *V* as the right
basis matrix from the SVD.

    mysvd <- svd(rand_mat[[1]])
    randV <- mysvd$v

Now for this *V*, we will first compute *U*<sub>*i*</sub>’s and
*Δ*<sub>*i*</sub> for different *D*<sub>*i*</sub> matrices in the
`mymat`. We achieve this by solving the linear equations:
*D*<sub>*i*</sub> = *U*<sub>*i*</sub>*Δ*<sub>*i*</sub>*V*<sup>*T*</sup>
for *i* = 1, …, 4. We then orthonormalize the columns of
*U*<sub>*i*</sub> using Proposition I.

    # get Ui and Delta for this newV
    out <- computeUDelta(mymat, randV)
    names(out)
    #> [1] "u"       "u_ortho" "d"       "d_ortho" "error"

The initial decomposition error is :

    calcDecompError(mymat, out$u_ortho, out$d, randV)
    #> [1] 22879.08

Now we will try to optimize using the new random *V* and corresponding
*U*<sub>*i*</sub>’s and *Δ*<sub>*i*</sub>’s.

    newopt <- optimizeFactorization(mymat, out$u_ortho, out$d, randV)
    # Number of updates taken
    newopt$error_pos
    #> [1] 220
    # New error
    newopt$error
    #> [1] 1411.555

We achieve the same factorization error (`1411.5550218`) after the
`optimizeFactorization` function call.

1.  Now instead of the right basis matrix from the SVD, we will set *V*
    as the left basis matrix.

<!-- -->

    mysvd <- svd(rand_mat[[1]])
    randV <- mysvd$u
    dim(randV)
    #> [1] 3 3

    # get Ui and Delta for this newV
    out <- computeUDelta(mymat, randV)
    calcDecompError(mymat, out$u_ortho, out$d, randV)
    #> [1] 13903.45

Now we will try to optimize with these matrices as our initial values.

    newopt <- optimizeFactorization(mymat, out$u_ortho, out$d, randV)
    # Number of updates taken
    newopt$error_pos
    #> [1] 283
    # New error
    newopt$error
    #> [1] 1411.555

Again we get the same decomposition error after optimizing.

1.  Instead of initial value being an orthogonal matrix, we will
    initialize *U*<sub>*i*</sub>’s, *Δ*<sub>*i*</sub>, and *V* with
    random matrices such that it does not guarantee

-   orthogonal property for *V* and
-   orthonormal columns for *U*<sub>*i*</sub>’s.

<!-- -->

    set.seed(111)
    # new random v
    newv <- createRandomMatrices(n = 1, ncols = 3, nrows = 3)[[1]]
    # seed value
    k <- 2392
    newu <- newd <- list()
    for(i in names(mymat)){
      myrow <- nrow(mymat[[i]])
      mycol <- ncol(mymat[[i]])
      set.seed(k)
      # new random u_i
      newu[[i]] <- createRandomMatrices(n = 1, ncols = mycol, nrows = myrow)[[1]]
      set.seed(k*2)
      # new random d_i
      newd[[i]] <- sample(1:1000, size = mycol)
      newmat <- newu[[i]] %*% diag(newd[[i]]) %*% t(newv)
      if (!qr(newmat)$rank == mycol)
        cat("\nNew matrix does not have full column rank")
      k = k+1
    }
    error <- calcDecompError(mymat, newu, newd, newv)
    cat("\nInitial error = ", error,"\n")
    #> 
    #> Initial error =  2.062531e+15

We see a very high factorization error because of the random
initialization.

    newopt <- optimizeFactorization(mymat, newu, newd, newv)
    newopt$error_pos
    #> [1] 142
    newopt$error
    #> [1] 1411.555

Again, we get the same factorization error after optimizing. Try
changing the seed value and compare the results.

This shows that the iterative update procedure converges and achieves
the same decomposition error regardless of the initial values.

### Estimating SVD

We will further demonstrate the case for *N* = 1 when we have just one
matrix. The `optimizeFactorization` function gives *U*<sub>*i*</sub>’s
with orthonormal column, *Δ*<sub>*i*</sub> a diagonal matrix, and an
orthogonal *V*. If the function converges, the results should be
identical to a standard SVD, except for the sign changes corresponding
to *U* and *V* columns. So we will compare the results from the
`optimizeFactorization` function with the standard SVD output. Let us
generate one example matrix say `newmat`.

    set.seed(171)
    newmat <- createRandomMatrices(n = 1, ncols = 3, nrows = 3)
    newmat
    #> $mat1
    #>      [,1] [,2] [,3]
    #> [1,]   41   10    6
    #> [2,]   64   85    8
    #> [3,]   82   87   57

1.  We will estimate the SVD of `newmat` using our iterative update
    function by setting the initial values to be identity matrix

<!-- -->

    newu <- newd <- list()
    newu[["mat1"]]  <- diag(3)
    newd[["mat1"]] <- diag(newu[["mat1"]])
    newu
    #> $mat1
    #>      [,1] [,2] [,3]
    #> [1,]    1    0    0
    #> [2,]    0    1    0
    #> [3,]    0    0    1
    newd
    #> $mat1
    #> [1] 1 1 1

The factorization error when initializing using identity matrix:

    calcDecompError(newmat, newu, newd, diag(3))
    #> [1] 30381

Let us optimize.

    opt_new <- optimizeFactorization(newmat, newu, newd, diag(3))
    cat("\n # of updates:", opt_new$error_pos,"\n")
    #> 
    #>  # of updates: 163
    opt_new$error
    #> [1] 0.0005937236

Error is close to zero. Let us compare the original matrix with the
reconstructed matrix based on the estimated *u*, *d* and *v* using the
`optimizeFactorization` function.

    newmat
    #> $mat1
    #>      [,1] [,2] [,3]
    #> [1,]   41   10    6
    #> [2,]   64   85    8
    #> [3,]   82   87   57
    opt_new$u[[1]] %*% diag(opt_new$d[[1]]) %*% t(opt_new$v)
    #>          [,1]     [,2]      [,3]
    #> [1,] 40.99568 10.00850  5.989020
    #> [2,] 64.01161 84.99197  7.993729
    #> [3,] 81.99198 87.00407 57.007920

    opt_new1 <- optimizeFactorization(newmat, newu, newd, diag(3), tol = 1e-21)
    cat("\n # of updates:", opt_new1$error_pos,"\n")
    #> 
    #>  # of updates: 559
    opt_new1$error
    #> [1] 7.402498e-13

    newmat
    #> $mat1
    #>      [,1] [,2] [,3]
    #> [1,]   41   10    6
    #> [2,]   64   85    8
    #> [3,]   82   87   57
    opt_new1$u[[1]] %*% diag(opt_new1$d[[1]]) %*% t(opt_new1$v)
    #>      [,1] [,2] [,3]
    #> [1,]   41   10    6
    #> [2,]   64   85    8
    #> [3,]   82   87   57

The reconstructed matrix is the same as the original matrix. Let us
compare the *U* and *V* with that from the standard SVD.

    newmat_svd <- svd(newmat[[1]])

    newmat_svd$d
    #> [1] 170.70126  31.96746  24.14876
    opt_new1$d
    #> $mat1
    #> [1]  24.14876 170.70126  31.96746

    newmat_svd$u
    #>            [,1]       [,2]        [,3]
    #> [1,] -0.2067092  0.1766241 -0.96232802
    #> [2,] -0.6071027 -0.7944740 -0.01541011
    #> [3,] -0.7672664  0.5810465  0.27145407
    opt_new1$u[[1]]
    #>             [,1]      [,2]       [,3]
    #> [1,]  0.96232803 0.2067092  0.1766241
    #> [2,]  0.01541005 0.6071027 -0.7944740
    #> [3,] -0.27145402 0.7672664  0.5810465

    newmat_svd$v
    #>            [,1]       [,2]       [,3]
    #> [1,] -0.6458388  0.1264120 -0.7529358
    #> [2,] -0.7054605 -0.4758901  0.5252181
    #> [3,] -0.2919209  0.8703727  0.3965270
    opt_new1$v
    #>            [,1]      [,2]       [,3]
    #> [1,]  0.7529358 0.6458388  0.1264119
    #> [2,] -0.5252182 0.7054605 -0.4758901
    #> [3,] -0.3965269 0.2919209  0.8703727

The results agree except for the sign and ordering of columns.

1.  Now we will estimate the SVD of `newmat` using our iterative update
    function from another random matrix with the same dimension.

<!-- -->

    set.seed(253)
    randmat_new <- createRandomMatrices(n = 1, ncols = 3, nrows = 3)
    randmat_new
    #> $mat1
    #>      [,1] [,2] [,3]
    #> [1,]   94   30   77
    #> [2,]   60   35  100
    #> [3,]   67   84   58
    newsvd <- svd(randmat_new[[1]])

Let us create a list for the *u* and *δ* matrices we just obtained from
the SVD of the random matrix. This allows us to use these matrices as
the initial values for the `optimizeFactorization` function.

    newu <- newd <- list()
    newu[[names(randmat_new)]] <- newsvd$u
    newd[[names(randmat_new)]] <- newsvd$d

The factorization error

    calcDecompError(newmat, newu, newd, newsvd$v)
    #> [1] 19465

Let us optimize.

    opt_new <- optimizeFactorization(newmat, newu, newd, newsvd$v)
    cat("\n # of updates:", opt_new$error_pos,"\n")
    #> 
    #>  # of updates: 235
    opt_new$error
    #> [1] 0.0006461125

Error is close to zero. Let us compare the original matrix with the
reconstructed matrix based on the estimated *u*, *d* and *v* using the
`optimizeFactorization` function.

    newmat
    #> $mat1
    #>      [,1] [,2] [,3]
    #> [1,]   41   10    6
    #> [2,]   64   85    8
    #> [3,]   82   87   57
    opt_new$u[[1]] %*% diag(opt_new$d[[1]]) %*% t(opt_new$v)
    #>          [,1]     [,2]      [,3]
    #> [1,] 40.99550 10.00886  5.988544
    #> [2,] 64.01211 84.99162  7.993461
    #> [3,] 81.99163 87.00424 57.008260

The estimated value is very close.

We can further improve our estimate by decreasing the tolerance
parameter (`tol`) in the optimization function.

    opt_new1 <- optimizeFactorization(newmat, newu, newd, newsvd$v, tol = 1e-21)
    cat("\n # of updates:", opt_new1$error_pos,"\n")
    #> 
    #>  # of updates: 673
    opt_new1$error
    #> [1] 9.156226e-14
    opt_new1$u[[1]] %*% diag(opt_new1$d[[1]]) %*% t(opt_new1$v)
    #>      [,1] [,2] [,3]
    #> [1,]   41   10    6
    #> [2,]   64   85    8
    #> [3,]   82   87   57

The reconstructed matrix is the same as the original matrix. Let us
compare the *U* and *V* with that from the standard SVD.

    newmat_svd <- svd(newmat[[1]])

    newmat_svd$u
    #>            [,1]       [,2]        [,3]
    #> [1,] -0.2067092  0.1766241 -0.96232802
    #> [2,] -0.6071027 -0.7944740 -0.01541011
    #> [3,] -0.7672664  0.5810465  0.27145407
    opt_new1$u[[1]]
    #>            [,1]       [,2]        [,3]
    #> [1,] -0.2067092 -0.1766241  0.96232803
    #> [2,] -0.6071027  0.7944740  0.01541009
    #> [3,] -0.7672664 -0.5810465 -0.27145405

    newmat_svd$v
    #>            [,1]       [,2]       [,3]
    #> [1,] -0.6458388  0.1264120 -0.7529358
    #> [2,] -0.7054605 -0.4758901  0.5252181
    #> [3,] -0.2919209  0.8703727  0.3965270
    opt_new1$v
    #>            [,1]       [,2]       [,3]
    #> [1,] -0.6458388 -0.1264120  0.7529358
    #> [2,] -0.7054605  0.4758901 -0.5252181
    #> [3,] -0.2919209 -0.8703727 -0.3965269

The results agree!

Cross-species gene expression analysis using A-SBF
==================================================

For cross-species gene expression datasets, we estimate the common space
*V* based on correlation (*R*<sub>*i*</sub>) between column phenotypes
(such as tissues, cell types, etc.) within a species. In our study, we
have shown that the inter-tissue gene expression correlation is similar
across species. Let
*X*<sub>*i*</sub> ∈ ℝ<sup>*m*<sub>*i*</sub> × *k*</sup> be a
standardized gene expression matrix where
*X*<sub>*i*</sub> = *C*<sub>*i*</sub>*D*<sub>*i*</sub>*S*<sub>*i*</sub><sup> − 1</sup>.
Here
*C*<sub>*i*</sub> = *I*<sub>*m*<sub>*i*</sub></sub> − *m*<sub>*i*</sub><sup> − 1</sup>1<sub>*m*<sub>*i*</sub></sub>1<sub>*m*<sub>*i*</sub></sub><sup>*T*</sup>
is a centering matrix and
*S*<sub>*i*</sub> = diag(*s*<sub>1</sub>, …, *s*<sub>*k*</sub>) is a
diagonal scaling matrix, where *s*<sub>*p*</sub> is the standard
deviation of *p*-th column of *D*<sub>*i*</sub>. The matrix
*X*<sub>*i*</sub> is a matrix with columns of *D*<sub>*i*</sub>
mean-centered and scaled by the standard deviation. The correlation
between expression profiles of *k* tissue types in species *i* is given
by
*R*<sub>*i*</sub> = *X*<sub>*i*</sub><sup>*T*</sup>*X*<sub>*i*</sub>/*m*<sub>*i*</sub>.
We then define an expected correlation matrix (𝔼(*R*<sub>*i*</sub>))
across *N* species as *M*, where *M* is defined as

$$
  M = \\frac{\\sum\_{i=1}^{N} R\_i}{N}.
$$

The shared right basis matrix *V* capturing the inter-tissue gene
expression correlation is determined from the eigenvalue decomposition
of *M*, where *M* = *V**Θ**V*<sup>*T*</sup>. Once the *V* and
*Δ*<sub>*i*</sub> are estimated using the SBF factorization, we compute
*U*<sub>*i*</sub> with orthonormal columns using proposition I. The
estimated *V* space captures inter-tissue gene expression correlation
relationship. For gene expression analysis, if we want the shared space
to represent inter-sample correlation relationship, we do not
update/change *V* while optimizing the factorization error. In such
cases, while reducing the factorization error we set `optimizeV = FALSE`
in the `optimizeFactorization` function.

Usage examples
--------------

Let us load the SBF package’s in-built gene expression dataset. The
dataset contains the average gene expression profile of five similar
tissues in three species.

    # load dataset
    avg_counts <- SBF::TissueExprSpecies
    # check the names of species
    names(avg_counts)
    #> [1] "Homo_sapiens"   "Macaca_mulatta" "Mus_musculus"

    # head for first species
    avg_counts[[names(avg_counts)[1]]][1:3, 1:3]
    #>                 hsapiens_brain hsapiens_heart hsapiens_kidney
    #> ENSG00000000003         2.3109         1.9414          5.2321
    #> ENSG00000000005         0.0254         0.2227          0.5317
    #> ENSG00000000419         5.2374         5.3901          5.5659

The number of genes annotated in different species is different. As a
result, the number of rows (genes) in the expression data will be
different for different species.

    sapply(avg_counts, dim)
    #>      Homo_sapiens Macaca_mulatta Mus_musculus
    #> [1,]        58676          30807        54446
    #> [2,]            5              5            5

Let us compute A-SBF with inter-tissue correlation.

    # A-SBF call using correlation matrix
    asbf_cor <- SBF(matrix_list = avg_counts, check_col_matching = TRUE,
                    col_index = 2, approximate = TRUE, transform_matrix = TRUE)
    # decomposition error
    asbf_cor$error
    #> [1] 65865.92

Optimize factorization to reduce decomposition error but by not updating
*V*.

    myopt_gef <- optimizeFactorization(avg_counts, asbf_cor$u_ortho, asbf_cor$delta,
                                       asbf_cor$v, optimizeV = FALSE)
    names(myopt_gef)
    #> [1] "u"         "v"         "d"         "error"     "error_pos" "error_vec"

    # new error
    myopt_gef$error
    #> [1] 63540.08
    # number of iterations
    myopt_gef$error_pos
    #> [1] 10

The number of iterations taken to optimize = 10.

    identical(asbf_cor$v, myopt_gef$v)
    #> [1] TRUE

Check whether estimated *U*<sub>*i*</sub>’s have orthonormal columns.

    zapsmall(t(as.matrix(myopt_gef$u[[names(myopt_gef$u)[1]]])) %*%
               as.matrix(myopt_gef$u[[names(myopt_gef$u)[1]]]))
    #>      [,1] [,2] [,3] [,4] [,5]
    #> [1,]    1    0    0    0    0
    #> [2,]    0    1    0    0    0
    #> [3,]    0    0    1    0    0
    #> [4,]    0    0    0    1    0
    #> [5,]    0    0    0    0    1

References
==========

Alter, Orly, Patrick O Brown, and David Botstein. 2000. “Singular Value
Decomposition for Genome-Wide Expression Data Processing and Modeling.”
*Proceedings of the National Academy of Sciences* 97 (18): 10101–6.

———. 2003. “<span class="nocase">Generalized singular value
decomposition for comparative analysis of genome-scale expression data
sets of two different organisms</span>.” *Proceedings of the National
Academy of Sciences* 100 (6): 3351–56.

Ponnapalli, Sri Priya, Michael A Saunders, Charles F Van Loan, and Orly
Alter. 2011. “<span class="nocase">A higher-order generalized singular
value decomposition for comparison of global mRNA expression from
multiple organisms</span>.” *PloS One* 6 (12): e28072.

Tamayo, Pablo, Daniel Scanfeld, Benjamin L Ebert, Michael A Gillette,
Charles WM Roberts, and Jill P Mesirov. 2007. “<span
class="nocase">Metagene projection for cross-platform, cross-species
characterization of global transcriptional states</span>.” *Proceedings
of the National Academy of Sciences* 104 (14): 5959–64.

Van Loan, Charles F. 1976. “Generalizing the Singular Value
Decomposition.” *SIAM Journal on Numerical Analysis* 13 (1): 76–83.
