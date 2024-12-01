% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{get_path}
\alias{get_path}
\title{Retrieve a Path from Predecessors Array}
\usage{
get_path(predecessors, origin_node, destination)
}
\arguments{
\item{predecessors}{A numeric vector where each element is the predecessor node
for a given node index, typically derived from a shortest path algorithm
like Dijkstra's.}

\item{origin_node}{The starting node index of the path.}

\item{destination}{The destination node index of the path.}
}
\value{
Returns a vector of node indices representing the path from the origin node
        to the destination node. If a path is not found (i.e., disconnected components),
        returns NULL.
}
\description{
This function reconstructs a path in a network from a vector of predecessors,
starting from a specified destination node and tracing back to the origin node.
The function returns the path as a sequence of node indices.
}
\examples{
\dontrun{
# Create a predecessors vector from a hypothetical Dijkstra's algorithm
predecessors <- c(6, 5, 5, 3, 6, -1, -1)
# Retrieve the path from node 1 to node 5
get_path(predecessors, origin_node = 1, destination = 5)
}
}
\keyword{internal}
