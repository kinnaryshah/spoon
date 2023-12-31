% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighted_nnSVG.R
\name{weighted_nnSVG}
\alias{weighted_nnSVG}
\title{Weighted nnSVG}
\usage{
weighted_nnSVG(input, assay_name = "logcounts", w)
}
\arguments{
\item{input}{either a SpatialExperiment object which contains a logcounts matrix, or a matrix}

\item{assay_name}{name of the assay in which the logcounts matrix is stored}

\item{w}{weights matrix}
}
\value{
either spe with weighted nnSVG statistics, or matrix with weighted nnSVG statistics
}
\description{
Run nnSVG for SVG detection using the weights
}
\details{
This function incorporates the weights for each observation to run nnSVG
}
\examples{
weighted_nnSVG(spe, "logcounts", w)

}
