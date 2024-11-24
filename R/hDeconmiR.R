#' @title
#' Hierarchical mode of DeconmiR (hDeconmiR)
#'
#' @aliases hDeconmiR
#'
#' @description
#' hDeconmiR is the hierarchical mode of DeconmiR for cell type proportion estimation.
#' It uses two generic miRNA references, a primary reference for the estimation of
#' solid tissues cell-types fractions for epithelial, fibroblast and the immune component,
#' and a separate secondary non-overlapping reference for the estimation of underlying
#' subtype fractions of one of the cell-type in the primary reference (mainly the immune
#' cell types).
#'
#' @param data.m miRNA Expression profile: A data matrix with rows labeling the
#' miRNAs (should use miRBase ID as in ref.m) and columns labeling samples
#' (e.g. primary tumor specimens). Missing value is not allowed and all values
#' should be log-transformed.
#'
#' @param ref1.m Reference profile: The signature matrix should be a miRNAs (rows)
#' by cell types (columns) matrix containing log-normalized expression
#' values of signature miRNAs.
#'
#' @param ref2.m
#' Similar to ref1.m, but now a matrix of secondary reference. For example, ref1.m
#' contains reference for epithelial cells, fibroblasts and total immune cells. Then
#' the ref2.m can be subtypes of immune cells, such as B-cells, NK cells, monocytes and etc.
#'
#' @param h.CT.idx
#' An index tells which cell-type in ref1.m is the higher order cell-types
#' in ref2.m. For example, ref1.m contains reference for the epithelial cells,
#' fibroblasts and total immune cells. Whereas ref2.m contains subtypes of immune
#' cells, then the h.CT.idx should be 3, corresponding to immune cells in ref1.m.
#'
#' @param method Deconvolution method: Choice of the deconvolution methods to
#' be used: ("RPC", "SVR", "CP", "NNLS", "SLE"). By default is the "RPC".
#'
#' @param maxit The limit of the number of IWLS iterations, only used for "RPC" mode.
#'
#' @param nu.v A vector of several candidate nu values, which is needed for
#' nu-classification, nu-regression, and one-classification in svm, only used
#' for "SVR" mode. The best estimation results among all candidate nu will be
#' automatically returned.
#'
#' @param constraint Choice of either 'inequality' or 'equality' normalization
#' constraint, only used in "CP" mode, By default is 'inequality' (i.e sum of
#' weights adds to a number less or equal than 1), which was implemented in
#' Houseman et al (2012).
#'
#' @return A matrix of the estimated fractions
#'
#' @examples
#' data(miRExpCount_BRCA.m)
#' miRexpRPM_BRCA.m <- getRPM(miRExpCount_BRCA.m)
#' miRexpNorm_BRCA.m <- logNormalize(miRexpRPM_BRCA.m)
#'
#' data(GenericMatrix)
#' ref1.m <- GenericMatrix$ref1
#' ref2.m <- GenericMatrix$ref2
#'
#' frac.m <- hDeconmiR(data.m = miRexpNorm_BRCA.m, ref1.m = ref1.m,
#'         ref2.m = ref2.m, h.CT.idx = 3, method = 'RPC')
#'
#'
#' @export
#'

hDeconmiR <- function(data.m, ref1.m, ref2.m, h.CT.idx, method =  c("RPC", "SVR", "CP", "NNLS", "SLE"),
    maxit = 50, nu.v = c(0.25, 0.5, 0.75), constraint = c("inequality", "equality")) {
    method <- match.arg(method)
    constraint <- match.arg(constraint)
    if (!method %in% c("RPC", "SVR", "CP", "NNLS", "SLE"))
        stop("Input a valid method!")
    if (method == "RPC") {
        frac1.m <- DoRPC(data.m, ref1.m, maxit)$estF
        frac2.m <- DoRPC(data.m, ref2.m, maxit)$estF
        frac.m <- cbind(frac1.m[, -h.CT.idx, drop = FALSE], frac1.m[, h.CT.idx] * frac2.m)
    } else if (method == "CBS") {
        frac1.m <- DoSVR(data.m, ref1.m, nu.v)$estF
        frac2.m <- DoSVR(data.m, ref2.m, nu.v)$estF
        frac.m <- cbind(frac1.m[, -h.CT.idx, drop = FALSE], frac1.m[, h.CT.idx] * frac2.m)
    } else if (method == "CP") {
        if (!constraint %in% c("inequality", "equality")) {
            # make sure constraint must be inequality or equality
            stop("constraint must be inequality or equality when using CP!")
        } else {
            frac1.m <- DoCP(data.m, ref1.m, constraint)$estF
            frac2.m <- DoCP(data.m, ref2.m, constraint)$estF
            frac.m <- cbind(frac1.m[, -h.CT.idx, drop = FALSE], frac1.m[, h.CT.idx] * frac2.m)
        }
    }else if (method == "NNLS") {
        frac1.m <- DoNNLS(data.m, ref1.m)$estF
        frac2.m <- DoNNLS(data.m, ref2.m)$estF
        frac.m <- cbind(frac1.m[, -h.CT.idx, drop = FALSE], frac1.m[, h.CT.idx] * frac2.m)
    }else if (method == "SLE") {
        frac1.m <- DoSLE(data.m, ref1.m)$estF
        frac2.m <- DoSLE(data.m, ref2.m)$estF
        frac.m <- cbind(frac1.m[, -h.CT.idx, drop = FALSE], frac1.m[, h.CT.idx] * frac2.m)
  }

    return(frac.m)
}





