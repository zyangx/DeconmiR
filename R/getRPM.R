#' Title Convert read counts to RPM
#'
#' @description Convert raw read count to RPM value
#'
#' @param count.m  miRNA-seq reads count table: a miRNAs (rows) by samples (columns) matrix.
#'
#' @return  a RPM (reads-per-million-miRNA-mapped) -normalized matrix.
#'
#' @export
#'
#' @examples
#' data(miRExpCount_LAML.m)
#' miRexpRPM_LAML.m <- getRPM(miRExpCount_LAML.m)
#'
getRPM <- function(count.m){
  RPM.m <- apply(count.m, 2, function(x) x/sum(x) * 1e6)
  return(RPM.m)
}
