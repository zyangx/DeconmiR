#' Title Log-transformation and normalization of miRNA expression
#'
#' @description Log-transform and further normalize the miRNA expression, which is
#' suitable normalization for deconvolution.
#'
#' @param miRExp.m miRNA expression table: a miRNAs (rows) by samples (columns) matrix.
#' Either RPM value or microarray based expression value.
#'
#' @param quantile Choice of either 'TRUE' or 'FALSE' for quantile normalization.
#' We suggest use normalization for microarray based expression.
#'
#' @return miRNA expression table with log2 transformed and normalized: a miRNAs
#' (rows) by samples (columns) matrix.
#'
#' @export
#'
#' @examples
#' data(miRExpCount_LAML.m)
#' miRexpRPM_LAML.m <- getRPM(miRExpCount_LAML.m)
#' miRexpNorm_LAML.m <- logNormalize(miRexpRPM_LAML.m)
#'
#' @importFrom limma normalizeBetweenArrays
logNormalize <- function(miRExp.m, quantile=FALSE){
  if (quantile == "TRUE") {
    logExp.m <- log2(miRExp.m + 1)
    logExpNorm.m <- normalizeBetweenArrays(logExp.m, method="quantile")
    return(logExpNorm.m)
  }else if(quantile == "FALSE") {
    logExpNorm.m <- log2(miRExp.m + 1)
    return(logExpNorm.m)
  }
}
