% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{quantile_cut}
\alias{quantile_cut}
\title{Stratify a vector into a factor with a desired number of levels by quantiles}
\usage{
quantile_cut(x, n)
}
\arguments{
\item{x}{A numeric vector}

\item{n}{The number of breaks}
}
\value{
A factor
}
\description{
Stratify a vector into a factor with a desired number of levels by quantiles
}
\examples{
x <- runif(1000)
g <- quantile_cut(x, n = 4)
}
\keyword{internal}
