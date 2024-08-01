
# DeconmiR
##Cell Type Deconvolution using miRNA Expression Profiling"

================

<!-- badges: start -->
<!-- badges: end -->

The **DeconmiR** package provides tools to infer the cell fractions a sample representing a mixture of such cell-types by using miRNA expression profiling. Inference proceeds via one of 5 methods (Robust Partial Correlations-**RPC**, Support Vector Regression-**SVR**, Constrained Projection-**CP**, Non-Negative Least Squares-**NNLS**, Solution of Linear Equations-**SLE**), as determined by the user. For now, the package contains 4 references, including three whole blood subtypes reference and one reference for breast tissue. This package  provides a unified interface to rapidly run and benchmark multiple deconvolution methods for miRNA expression based deconvolution.

## Installation

You can install the development version of DeconmiR like so:

``` r
devtools::install_github("zyangx/DeconmiR")
```

## Example

This is a basic example which shows you how to estimate cell proportions from miRNA profiling:

``` r
library(DeconmiR)
## basic example code

data(miRExpCount_LAML.m)
data(refMatrix_Blood1.m)
miRexpRPM_LAML.m <- getRPM(miRExpCount_LAML.m)
miRexpNorm_LAML.m <- logNormalize(miRexpRPM_LAML.m)
est.o <- DeconmiR(miRexpNorm_LAML.m, refMatrix_Blood1.m, method = 'RPC')
BloodFrac.m <- est.o$estF
BloodFrac.m 
```

Detailed example for package usage can be seen in the vignettes.

