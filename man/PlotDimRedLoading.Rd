% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TockyPCA.R
\name{PlotDimRedLoading}
\alias{PlotDimRedLoading}
\title{Set a gate to define negative (and positive) for each marker expression}
\usage{
PlotDimRedLoading(x, reduction = "PCA")
}
\arguments{
\item{x}{A TockyPrepData.}

\item{reduction}{Choose whether to use PCA or CCA as a reduction method.}
}
\value{
A TockyPrepData (unchanged) for safety.
}
\description{
Set a gate to define negative (and positive) for each marker expression
}
\examples{
\dontrun{
PlotDimRedLoading(x)
}
}
