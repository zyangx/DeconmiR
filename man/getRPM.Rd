% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRPM.R
\name{getRPM}
\alias{getRPM}
\title{Title Convert read counts to RPM}
\usage{
getRPM(count.m)
}
\arguments{
\item{count.m}{miRNA-seq reads count table: a miRNAs (rows) by samples (columns) matrix.}
}
\value{
a RPM (reads-per-million-miRNA-mapped) -normalized matrix.
}
\description{
Convert raw read count to RPM value
}
\examples{
data(miRExpCount_LAML.m)
miRexpRPM_LAML.m <- getRPM(miRExpCount_LAML.m)

}
