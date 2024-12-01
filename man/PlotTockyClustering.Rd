% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TockyClustering.R
\name{PlotTockyClustering}
\alias{PlotTockyClustering}
\title{Plot Tocky Clusters}
\usage{
PlotTockyClustering(
  x,
  jpeg = FALSE,
  max_cells_displayed = 30000,
  filename = NULL
)
}
\arguments{
\item{x}{A TockyPrepData after running the function TockyClustering
This function will generate PCA plots with cluster and Angle data}

\item{jpeg}{Whether to out a jpeg file. The default is pdf = FALSE, by which a jpeg file is produced.}

\item{max_cells_displayed}{The number of cells displayed in plots.}

\item{filename}{A character string for file name}
}
\value{
The slot Reduction will contain the new slot Tocky_clusters, which includes kmeans clustering result
}
\description{
Plot Tocky Clusters
}
\examples{
\dontrun{
PlotTockyClustering(x)
}
}
