% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotPercentageClusters.R
\name{GetStatsClusterPercentage}
\alias{GetStatsClusterPercentage}
\title{Retrieve Cluster Percentage Data and Stats}
\usage{
GetStatsClusterPercentage(
  x,
  writeResults = FALSE,
  filename = "cluster_percentage"
)
}
\arguments{
\item{x}{A TockyPrepData that has been processed with the PlotClusterPercentage function.}

\item{writeResults}{A logical value. If TRUE, two files will be generated containing group statistics and p-values, respectively. If FALSE, these statistical results are displayed in the Terminal.}

\item{filename}{(optional) Base name for the output files when writeResults is TRUE.}
}
\value{
The same TockyPrepData passed as input, for consistency in function design, though the function primarily focuses on data retrieval and display or file writing.
}
\description{
This function retrieves the cluster percentage data from a TockyPrepData that has already been processed with the PlotClusterPercentage function. It can display the statistics in the Terminal or write them to CSV files.
}
\examples{
\dontrun{
  GetStatsClusterPercentage(x)
}
}
