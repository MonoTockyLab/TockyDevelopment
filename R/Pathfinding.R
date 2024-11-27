# ==============================================================================
# Title:       TockyDevelopment
# Description: Multidimensional Analysis and Pathfinding Analysis of High-Dimensional Flow Cytometric Fluorescent Timer Data
# Version:     0.1.0
# Author:      Masahiro Ono
# Created:     27 November 2024
# Modified:    27 November 2024
#
# Copyright (C) 2024 Masahiro Ono
#
# NOTICE:  All rights are reserved, including all intellectual property and patent rights.
# A patent application has been filed related to the methodologies employed within this code.
#
# The code is available on GitHub without a standard licensing option, intended for
# public viewing and verification related to the associated academic publication. No
# rights are granted for the use, modification, or distribution of the code for any
# purposes without explicit permission from Masahiro Ono, Imperial College London.
#
# For permissions or inquiries, please contact: m.ono@imperial.ac.uk
#
# This software is distributed WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# ==============================================================================

#' Dijkstra-Tocky Pathfinding
#'
#' This function identifies paths within a Tocky Network using Dijkstra-Tocky Pathfinding
#'
#' @param x A TockyPrepData containing the network and relevant data.
#' @param origin_node The starting node (cluster) for the path.
#' @param destination_node The ending node (cluster) for the path.
#' @param verbose Logical indicating whether to print progress messages and outputs.

#'
#' @return The original TockyPrepData is returned with an additional list attached
#'         containing the ordered paths and the closest path.
#'
#' @details The function computes mean TimerAngles for each cluster, identifies
#'          all possible increasing paths from the origin to the destination,
#'          and then calculates the total weight for each path to find the closest one.
#'          Paths are determined by increasing TimerAngle values, ensuring that each
#'          step along a path moves to a node with a higher TimerAngle.
#'
#' @examples
#' \dontrun{
#' x <- DijkstraTockyPath(x, origin_node = '3', destination_node = '21')
#' }
#' @export
#' @importFrom methods show
#' @importClassesFrom TockyPrep TockyPrepData
DijkstraTockyPath <- function(x, origin_node, destination_node, verbose = TRUE) {
    
    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
        
    }
    
    if(is.null(x@Tocky$network[['CCA']])) {
        stop("No network data found. Run TockyNetwork using CCA. \n")
    }
    set.seed(123)
    
    network <- x@Tocky$network[['CCA']]
    data <- x@Data
    lg_na <- is.na(data$Angle)
    data <- data[!lg_na,]
    clusters <- x@Tocky$PCAclusters[!lg_na]
    
    timerAngle <- tapply(data$Angle, clusters, mean)
    node_names <- as.numeric(names(timerAngle))
    sorted_indices <- order(node_names)
    timerAngle <- timerAngle[sorted_indices]

    adjMatrix <- x@Tocky$adjacency_matrix[['CCA']]
    node_names <- as.numeric(colnames(adjMatrix))
    sorted_indices <- order(node_names)
    adjMatrix <- adjMatrix[sorted_indices, sorted_indices]
    
    adjMatrix <- adjMatrix
    diag(adjMatrix) <- 0

    dijkstra_with_timerangle <- function(adjMatrix, timerAngle, source) {
      numVertices <- nrow(adjMatrix)
      distances <- rep(Inf, numVertices)
      distances[source] <- 0
      visited <- rep(FALSE, numVertices)
      predecessors <- rep(-1, numVertices)

      for (count in 1:(numVertices-1)) {
        u <- which.min(ifelse(visited, Inf, distances))
        visited[u] <- TRUE
        
        for (v in 1:numVertices) {
          if (!visited[v] && adjMatrix[u, v] > 0 &&  distances[u] + adjMatrix[u, v] < distances[v] && timerAngle[v] > timerAngle[u]) {#
            distances[v] <- distances[u] + adjMatrix[u, v]
            predecessors[v] <- u
          }
        }
      }
      
      list(distances = distances, predecessors = predecessors)
    }


    result <- dijkstra_with_timerangle(adjMatrix, timerAngle, source = as.numeric(origin_node))
    closest_path <- get_path(result$predecessors, origin_node = as.numeric(origin_node), destination = as.numeric(destination_node))
    x@Tocky$network$DijkstraTockyPath <- result
    x@Tocky$network$DijkstraTockyPath$closest_path <- closest_path
    
    if (is.null(closest_path)) {
        cat("No Dijkstra-Tocky Path was identified. \n")
        return(invisible(x))
    } else {

        if(verbose){
            cat(paste('The optimised Tocky Angle gadient path is: \n'))
            show(closest_path)
        }

        
        path_name <- paste(origin_node, 'to', destination_node)
        
        plot_network_path(network, path = closest_path, title = path_name)
        
    }
    
    return(invisible(x))
}


#' Display Closest Tocky Angle Gradient Path
#'
#' Retrieves and displays the closest Tocky Angle Gradient Path from a TockyPrepData.
#' This function is intended for use after running DijkstraTockyPath to quickly
#' access and review the computed closest path based on increasing TimerAngle values.
#'
#' @param x A TockyPrepData that has already been processed with DijkstraTockyPath
#'        to compute the Tocky Angle Gradient Path.
#' @param origin_node The starting node (cluster) for the path.
#' @param destination_node The ending node (cluster) for the path.

#'
#' @return Prints the closest path and also returns it for further analysis or usage.
#'
#' @examples
#' \dontrun{
#' showDijkstraTockyPath(x)
#' }
#' @export
#' @importClassesFrom TockyPrep TockyPrepData
showDijkstraTockyPath <- function(x, origin_node, destination_node) {
  if (!inherits(x, "TockyPrepData")) {
    stop("Use a TockyPrepData \n.")
  }
  if (is.null(x@Tocky$network$DijkstraTockyPath)) {
    stop("Run DijkstraTockyPath \n.")
  }
  result <- x@Tocky$network$DijkstraTockyPath
  closest_path <- get_path(result$predecessors, origin_node = as.numeric(origin_node), destination = as.numeric(destination_node))
  cat("Dijkstra-Tocky Path is:\n")
  print(closest_path)
  
  network <- x@Tocky$network[['CCA']]
  plot_network_path(network, path = closest_path)
  
  return(closest_path)
}




