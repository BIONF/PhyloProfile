## ------------------------------------------------------------------------
set.seed(123)
n = 100  # matrix size
k = 5    # number of eigenvalues to calculate
# Some random data
M = matrix(rnorm(n^2), n)
# Make it symmetric
A = crossprod(M)
# Show its largest 5 eigenvalues and associated eigenvectors
head(eigen(A)$values, 5)
head(eigen(A)$vectors[, 1:5])

## ------------------------------------------------------------------------
library(RSpectra)
res = eigs_sym(A, k, which = "LM")  # "LM" is the default
res$values
head(res$vectors)

## ------------------------------------------------------------------------
eigs_sym(A, k, opts = list(retvec = FALSE))

## ------------------------------------------------------------------------
library(Matrix)
Msp = as(M, "dgCMatrix")
Asp = as(A, "dgRMatrix")

eigs(Msp, k, which = "LR", opts = list(retvec = FALSE))$values  # largest real part
eigs_sym(Asp, k, opts = list(retvec = FALSE))$values

## ------------------------------------------------------------------------
# Implicitly define the matrix by a function that calculates A %*% x
# Below represents a diagonal matrix whose elements are stored in the `args` parameter
f = function(x, args)
{
    # diag(args) %*% x == x * args
    return(x * args)
}
eigs_sym(f, k = 3, n = 10, args = 1:10)

## ------------------------------------------------------------------------
eigs_sym(A, k, which = "LM", sigma = 0)$values  # recommended way
eigs_sym(A, k, which = "SM")$values             # not recommended

## ------------------------------------------------------------------------
set.seed(123)
m = 100
n = 20
k = 5
A = matrix(rnorm(m * n), m)

str(svds(A, k, nu = k, nv = k))

## ------------------------------------------------------------------------
Asp = as(A, "dgCMatrix")
svds(Asp, k, nu = 0, nv = 0)

## ------------------------------------------------------------------------
Af      = function(x, args)  as.numeric(args %*% x)
Atransf = function(x, args)  as.numeric(crossprod(args, x))
str(svds(Af, k, Atrans = Atransf, dim = c(m, n), args = Asp))

