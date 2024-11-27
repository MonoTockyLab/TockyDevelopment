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

#' Plot Network Path
#'
#' Visualizes a network using the Fruchterman-Reingold layout, with the option to highlight a specific path and customize vertex colors.
#' Nodes and edges along the path are highlighted to distinguish them from the rest of the network.
#'
#' @param network An igraph network object to be visualized. The function expects this object to represent a graph with nodes and edges properly defined.
#' @param path Optional; a vector of node names indicating the sequence of nodes along the path to be highlighted. If provided, this path will be visually distinguished in the plot.
#' @param vertex_col A character string indicating the default color for vertices not in the path. This allows for customization of the network's appearance.
#' @param title A character string providing a title for the network plot. This allows for basic customization of the plot's appearance.
#'
#' @return Invisibly returns NULL. The function is used for its side effect of plotting the network.
#'
#' @examples
#' \dontrun{
#' plot_network_path(net, path = c("node1", "node2", "node3"),
#' vertex_col = "lightblue", title = "Sample Network Path")
#' }
#' @importFrom igraph layout_with_fr V E vcount ecount get_edge_ids
#' @keywords internal
plot_network_path <- function(network, path = NULL, vertex_col = 'grey', title = 'Angle Path') {
    layout <- layout_with_fr(network)

    vertex.color <- rep(vertex_col, vcount(network))
    names(vertex.color) <- V(network)$name
    edge.color <- rep("grey", ecount(network))
    vertex.size <- rep(10, vcount(network))
    names(vertex.size) <- V(network)$name###

    edge.width <- rep(1, ecount(network))
    label.dist <- 0
    label.degree <- 0
    
    if (!is.null(path) && length(path) > 1) {
        path <- as.character(path)
        path_edges <- c()
        for (i in 1:(length(path)-1)) {
            path_edges <- c(path_edges, get_edge_ids(network, c(path[i], path[i+1])))
        }
        
        edge.color[path_edges] <- "red"
        edge.width[path_edges] <- 2
        
        vertex.color[path] <- "pink"
        vertex.size[path] <- 10
    }

    plot(network, vertex.color = vertex.color, vertex.frame.color = "black", vertex.shape = "circle",
         vertex.size = vertex.size, vertex.label.color = "black", vertex.label.family = "Helvetica",
         vertex.label.font = 1, vertex.label.cex = 0.7, vertex.label.dist = 0, vertex.label.degree = 0,
         edge.color = edge.color, edge.width = edge.width, edge.arrow.size = 1, edge.arrow.width = 1,
         edge.lty = "solid", edge.curved = 0, main = title, layout = layout)
}



#' Retrieve a Path from Predecessors Array
#'
#' This function reconstructs a path in a network from a vector of predecessors,
#' starting from a specified destination node and tracing back to the origin node.
#' The function returns the path as a sequence of node indices.
#'
#' @param predecessors A numeric vector where each element is the predecessor node
#'        for a given node index, typically derived from a shortest path algorithm
#'        like Dijkstra's.
#' @param origin_node The starting node index of the path.
#' @param destination The destination node index of the path.
#'
#' @return Returns a vector of node indices representing the path from the origin node
#'         to the destination node. If a path is not found (i.e., disconnected components),
#'         returns NULL.
#'
#' @examples
#' \dontrun{
#' # Create a predecessors vector from a hypothetical Dijkstra's algorithm
#' predecessors <- c(6, 5, 5, 3, 6, -1, -1)
#' # Retrieve the path from node 1 to node 5
#' get_path(predecessors, origin_node = 1, destination = 5)
#' }
#' @keywords internal
get_path <- function(predecessors, origin_node, destination) {
    path <- numeric()
    current_vertex <- destination

    while (current_vertex != origin_node) {
      if (current_vertex == -1) {
        return(NULL)
      }
      
      path <- c(current_vertex, path)
      current_vertex <- predecessors[current_vertex]
    }
    
    path <- c(origin_node, path)
    
    return(path)
}



#' Calculate the number of cells in each cluster in each sample using a Tocky Object.
#' @param x A TockyPrepData after running the function TockyClustering.
#' @return The slot Reduction will contain the new slot cluster_cellnumbers, which includes the cell number of each cluster in each sample
#' @keywords internal
#' @examples
#' \dontrun{
#' x <- cluster_cell_num(x)
#'}
#'


cluster_cell_num <- function(x){
    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
    }
    if(is.null(x@Data$Locus)){
        stop("Perform TockyLocus before using this function. \n")
    }
    
    if(is.null(x@Tocky$PCAclusters)){
        stop("Perform TockyPCA and TockyClustering before using this function. \n")
    }
    
    data <- x@Data
    clusters <- x@Tocky$PCAclusters
    
    X <- split(x@Data, clusters)
    total_file <- x@Data$file
    total_length <- tapply(x@Data$Angle, total_file, length)
    
    tdf <- c()
    for(i in 1:length(X)){
        t_file <- X[[i]]$file
        t_len <- tapply(X[[i]]$Angle, t_file, length)
        tdf <- rbind(tdf, t_len)
    }
    
    cluster_cellnum <- tdf
    rownames(cluster_cellnum) <- names(X)
    colnames(cluster_cellnum) <- names(t_len)
    x@Tocky[['cluster_cellnumbers']] <- cluster_cellnum
    
    return(x)
    
}


#' To produce distinct colour codes
#' @param x A numeric vector
#' @return a character vector for colour code
#' @keywords internal
#' @examples
#' \dontrun{
#' col  <-  color_palette(x)
#'}
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom grDevices rainbow

color_palette <- function(x){
    if(!is.numeric(x)){
        stop("Use a numeric vector. \n")
    }
    l <- length(x)
    n <- round(l/11, digits = 0)
    pallett <- rainbow(11)
    for(j in 1:(n+3)){
            j <- sample(1:nrow(brewer.pal.info), 1)
            pallett <- c(pallett, brewer.pal(rownames(brewer.pal.info)[j], n=brewer.pal.info[j,'maxcolors']))
        }
    pallett <- unique(pallett)
    pallett <- sample(pallett)
    return(pallett)
    
}



#' Stratify a vector into a factor with a desired number of levels by quantiles
#' @param x A numeric vector
#' @param n The number of breaks
#' @return A factor
#' @export
#' @examples
#' x <- runif(1000)
#' g <- quantile_cut(x, n = 4)
#' @keywords internal

quantile_cut <- function(x, n) {
    breaks <- quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE)
    
    out <- cut(x, breaks = breaks, include.lowest = TRUE)
    
    return(out)
    
}



#' Convert Angle into colour code
#' @param x Angle numeric vector
#' @param alpha A numeric vector between 0 and 1 for transparency
#' @return a character vector for colour code
#' @export
#' @examples
#' \dontrun{
#' col  <-  angle_to_colour(x)
#'}
#'
#' @keywords internal
angle_to_colour <- function(x, alpha= 0.3){
    nalogic  <-  is.na(x)
    y <- x
    na_lg <- is.na(x)
    angle <- x[!na_lg]
    angle <- angle/90
    r <- sin(angle)
    b <- cos(angle)
    colour <- rgb(r, 0, b, alpha = alpha)
    y[!na_lg] <- colour
    y[na_lg] <- rgb(0, 0, 0, alpha = 0)
    return(y)
}


#' To produce colour key as a gradient rectangle
#' @param expression A numeric vector defining the expression to be plotted.
#' @param x A number. The x position of the lower left corner of the rectangle
#' @param y A number. The y position of the lower left corner of the rectangle
#' @param width A number. The width of the rectangle
#' @param height A number. The height of the rectangle
#' @param colour Either 'Spectral' or 'BlueRed' for Angle colour key.
#' @param method Either 'Angle' or 'Spectral' for Angle and any expression, respectively.
#' @export
#' @examples
#' \dontrun{
#' plot(1:5)
#' plot_color_code(x = 3, y = 4.5, width = 1, height = 0.3, colour = 'Spectral', method = 'Angle')
#' }
#' @importFrom RColorBrewer brewer.pal
plot_color_code <- function(expression = NULL, x, y, width, height, colour = 'Spectral', method = 'Expression') {
    
    min_expr <- min(expression, na.rm = TRUE)
    max_expr <- max(expression, na.rm = TRUE)
    
    if (colour == 'BlueRed') {
        gradient_colors <- colorRampPalette(c("blue", "white", "red"))(50)
    } else if (colour == 'Spectral') {
        gradient_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = colour)))(50)
    }
    
    rect(x, y, x + width, y + height, col = gradient_colors, border = NA)
    rect_width <- width / length(gradient_colors)
    
    for (i in 1:length(gradient_colors)) {
        rect_x <- x + (i - 1) * rect_width
        rect_y <- y
        rect(rect_x, rect_y, rect_x + rect_width, rect_y + height, col = gradient_colors[i], border = NA)
    }
    
    text(x + 0.5 * width, y + height * 1.2, labels = 'Colour Key', cex = 0.7)
    
    if (method == 'Expression') {
        min_lab <- round(min_expr, 1)
        mid_lab <- round((min_expr + max_expr) / 2, 1)
        max_lab <- round(max_expr, 1)
        text(x, y, labels = min_lab, cex = 0.5)
        text(x + 0.5 * width, y, labels = mid_lab, cex = 0.5)
        text(x + width, y, labels = max_lab, cex = 0.5)
    }
}
