
# DeconmiR
##Cell Type Deconvolution using miRNA Expression Profliling

================

<!-- badges: start -->
<!-- badges: end -->

The **DeconmiR** package provides tools to infer the cell fractions a sample representing a mixture of such cell-types by using miRNA expression profiling. Inference proceeds via one of 5 methods (Robust Partial Correlations-**RPC**, Support Vector Regression-**SVR**, Constrained Projection-**CP**, Non-Negative Least Squares-**NNLS**, Solution of Linear Equations-**SLE**), as determined by the user. For now, the package contains three references for whole blood tissue and nine references for solid tissues, including breast, lung, liver, kidney, prostate, esophagus, respiratory tract, vessel and skin. In addition, it also provides a generic reference which was used for cell type deconvolution of any solid tissue. This package provides a unified interface to rapidly run and benchmark multiple deconvolution methods for miRNA expression based deconvolution.

## Installation

You can install the development version of DeconmiR like so:

``` r
devtools::install_github("zyangx/DeconmiR")
```

## Example
We show an example of using this package to estimate immune cell-type fractions in adult whole blood.

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

DeconmiR also provides functions for cell type deconvolutions of solid tissues. The detailed example for package usage can be seen in the vignettes.
