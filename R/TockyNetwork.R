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


#' Produce A Network of Tocky Clusters
#'
#' @param x A TockyPrepData after running the function TockyClustering.
#' @param reduction Choose whether to use PCA or CCA as a reduction method. When using CCA, TockyCCA must have been applied to the TockyPrepData.
#' @param cut_off Threshold value as a quantile percentage for edge connection. For example, the default 0.2 will set the threshold at 20\% quantile of distance between clusters, connecting neighbor clusters.
#' @return A revised TockyPrepData containing an igraph network object. Network is constructed based on the neighbor proximity.
#' @export
#' @examples
#' \dontrun{
#' x <- TockyNetwork(x)
#' }
#' @importFrom igraph E V plot.igraph components
#' @importFrom stats quantile dist
#' @importClassesFrom TockyPrep TockyPrepData

TockyNetwork <- function(x, reduction = 'CCA', cut_off = 0.2){
    
    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
        
    }
    
    if(reduction == 'CCA'){
        if(is.null(x@Tocky$CCA)) {
            stop("Please run TockyCCA first. \n")
        }
    }
    
    if(is.null(x@Tocky$PCAclusters)) {
        stop("Please run TockyClustering first.")
    }
    
    total_indices <- seq_len(nrow(x@Data))
    clusters <- x@Tocky$PCAclusters
    clusters <- as.vector(clusters)
    clusters <- as.factor(clusters)
    names(clusters) <- total_indices
    
    if (reduction == 'PCA') {
        data_pca <- x@Tocky$PCA$x
        lg <- !is.na(x@Data$Angle)
        angle_data <- x@Data$Angle
    } else if (reduction == 'CCA') {
        used_indices <- x@Tocky$CCA$used_indices
        data_pca <- x@Tocky$CCA$CCA$wa
        data_pca <- as.matrix(data_pca)
        lg <- !is.na(x@Data[used_indices, 'Angle'])
        angle_data <- x@Data[used_indices, 'Angle']
        clusters <- clusters[used_indices]
        clusters <- as.factor(as.vector(clusters))
    }
    
    barycentre_x <- tapply(data_pca[, 1], clusters, mean)
    barycentre_y <- tapply(data_pca[, 2], clusters, mean)
    barycentre <- data.frame(CCA1 = barycentre_x, CCA2 = barycentre_y)
    tmp <- as.matrix(dist(barycentre))
    threshold <- quantile(tmp, cut_off)
    adjacency_matrix <- tmp
    
    for (i in 1:nrow(tmp)) {
        for (j in 1:nrow(tmp)) {
            if (tmp[i, j] <= threshold) {
                adjacency_matrix[i, j] <- 1 / adjacency_matrix[i, j]
            } else {
                adjacency_matrix[i, j] <- 0
            }
        }
    }

    network <- graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE, mode = "undirected", diag = FALSE)
    
    col <- tapply(angle_data[lg], clusters[lg], mean)
    col <- angle_to_colour(col, alpha = 0.4)
    col[!lg] <- 0
    
    components_info <- components(network)
    num_modules <- sum(components_info$csize > 1)
    num_isolated_clusters <- sum(components_info$csize == 1)
    
    print(paste("Number of modules:", num_modules))
    print(paste("Number of isolated clusters:", num_isolated_clusters))

    x@Tocky$network[[reduction]] <- network
    x@Tocky$adjacency_matrix[[reduction]] <- adjacency_matrix
    return(x)
}

 


#' Plot A Tocky Network
#' @param x A TockyPrepData after running the function TockyNetwork.
#' @param reduction Choose whether to use PCA or CCA as a reduction method.
#' @param edge_scale A scale factor for edge thickness.
#' @param select_variable Whether to choose a variable (marker expression). The default is FALSE and produces Tocky Angle and Intensity.
#' @param log2fc Logical. If TRUE, the colours of nodes are determined by log2 fold change between the two groups.
#' @param p_adjust If p_adjust is "none", p-value adjustment is not used.
#' Note that p-values from two-group comparisons are stored within the TockyPrepData object and have been calculated by 'PlotClusterPercentage'.
#' @return Network plot is produced, showing clusters as vertices and distance between clusters as edges.
#' @param mds If Multidimensional Scaling is used for constructing network graph.
#' @param reflect If Multidimensional Scaling is used, graph can be refelected by either 'horizontal' or 'vertical'.
#' @export
#' @examples
#' \dontrun{
#' PlotTockyNetwork(x)
#'}
#' @importFrom igraph E V E<- V<- graph_from_adjacency_matrix plot.igraph
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices col2rgb
#' @importFrom igraph layout_with_mds
#' @importClassesFrom TockyPrep TockyPrepData

PlotTockyNetwork <- function(x, reduction = 'CCA', select_variable = FALSE, edge_scale = NULL, mds = FALSE, reflect = NULL, log2fc = FALSE, p_adjust= NULL){

    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
        
    }
    if(length(x@Tocky$network[[reduction]])==0){
        stop("Apply TockyNetwork to x before using this function. \n")
        
    }
    
    intensity_to_colour <- function(intensity, low_col = "blue", high_col = "red", alpha = 0.3) {
        if (any(is.na(intensity))) {
            intensity[is.na(intensity)] <- 0  # Handle NA values, if necessary
        }
        
        colour <- colorRampPalette(c(low_col, high_col))(length(intensity))
        rgba_col <- apply(col2rgb(colour), 2, function(col) {
            rgb(col[1], col[2], col[3], maxColorValue = 255, alpha = alpha)
        })
        
        return(rgba_col)
    }
    
    prepare_layout <- function(network, use_mds = FALSE, reflect = FALSE) {
        if (!use_mds) {
            return(NULL)
        } else {
            layout <- layout_with_mds(network, dim = 2)
            if (!is.null(reflect)) {
                if (reflect == 'horizontal') {
                    layout[, 1] <- -layout[, 1]
                } else if (reflect == 'vertical') {
                    layout[, 2] <- -layout[, 2]
                } else {
                    stop("Choose either 'horizontal' or 'vertical' for the argument reflect.")
                }
            }
            return(layout)
        }
    }

    calculate_log2fc_color <- function(group_stats, low_col = "blue", mid_col = "white", high_col = "red", group_red = NULL) {
        groups <- unique(group_stats$group)
        if(!is.null(group_red)){
            tmpgroup <- groups
            groups <- c(setdiff(tmpgroup, group_red), group_red)

        }
        if (length(groups) != 2) {
            stop("The function is designed to work with exactly two groups.")
        }
        
        group1_stats <- group_stats[group_stats$group == groups[1], ]
        group2_stats <- group_stats[group_stats$group == groups[2], ]
        
        group1_stats$cluster <- as.numeric(as.character(group1_stats$cluster))
        group2_stats$cluster <- as.numeric(as.character(group2_stats$cluster))
        group1_stats <- group1_stats[order(group1_stats$cluster), ]
        group2_stats <- group2_stats[order(group2_stats$cluster), ]
        
        log2fc <- log2(group2_stats$average_percentage + 1) - log2(group1_stats$average_percentage + 1)
        
        max_fc <- max(abs(log2fc))
        color_scale <- colorRampPalette(c(low_col, mid_col, high_col))(100)
        scaled_fc <- round((log2fc + max_fc) / (2 * max_fc) * 99) + 1
        colors <- color_scale[scaled_fc]

        
        out <- list(log2fc = log2fc, colours = colors)
        return(out)
    }
    
    cat(paste(reduction, 'is used to plot a TockyNetwork... \n'))
    network <- x@Tocky$network[[reduction]]
    
    clusters <- x@Tocky$PCAclusters
    clusters <- as.vector(clusters)
    clusters <- as.factor(clusters)
    data <- x@Data
    
    if (reduction == 'CCA') {
        used_indices <- which(rownames(x@Data) %in% rownames(x@Tocky$CCA$CCA$wa))
        clusters <- clusters[used_indices]
        data <- data[used_indices,]
    }
    
    if (reduction == 'PCA') {
        lg <- !is.na(x@Data$Angle)
        angle_data <- x@Data$Angle[lg]
        clusters_for_color <- clusters[lg]
    } else if (reduction == 'CCA') {
        angle_data <- x@Data[used_indices, 'Angle']
        clusters_for_color <- as.factor(clusters)
    }
    
    col <- angle_values <- tapply(angle_data, clusters_for_color, mean)
    cat("Plotting Tocky Angle Network... \n")
    
    col2 <- cut(col, breaks = 49)
    col2 <- as.factor(col2)
    levels(col2) <- colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)
    col2 <- as.character(col2)
    col2[is.na(col2)] <- rgb(0, 0, 0, alpha = 0.01)
    names(col2) <- names(col)
    col2 <- col2[names(V(network))]
    angle_col <- col2
    
    par(mfrow=c(1,2))
    E(network)$weight[is.infinite(E(network)$weight)] <- 20
    
    tmp <- exp(E(network)$weight)
    tmp[tmp> 5] <- 5
    z_scores <- tmp / max(tmp)
    if(is.null(edge_scale)){
        edge_scale <- 1
    }
    rescaled_weights <- edge_scale * z_scores
    
    plot_network <- function(network, vertex_col, title, layout) {
        set.seed(123)
        plot(network, vertex.color = vertex_col, vertex.frame.color = "black", vertex.shape = "circle", vertex.size = 10,
        vertex.label.color = "black", vertex.label.family = "Helvetica", vertex.label.font = 1, vertex.label.cex = 0.7,
        vertex.label.dist = 0, vertex.label.degree = 0, edge.color = "black", edge.width = rescaled_weights,
        edge.arrow.size = 1, edge.arrow.width = 1, edge.lty = "solid", edge.curved = 0, main = title, layout = layout)
    }

    layout <- prepare_layout(network, mds, reflect)
    key_rect_x <- par()$usr[2] - 0.41
    key_rect_y <- par()$usr[4] - 0.13
    key_rect_width <- 0.4
    key_rect_height <- 0.125

    if (!select_variable && !log2fc) {
        plot_network(network, angle_col, 'Tocky Angle', layout)
        par(xpd = TRUE)
        plot_color_code(expression = angle_values, method = 'Expression', x = key_rect_x, y = key_rect_y, width = key_rect_width, height = key_rect_height, colour = 'Spectral')

    } else  if (log2fc) {
        group_stats <- x@Tocky$PlotClusterPercentage$group_statistics
        tmpgroup <- unique(group_stats$group)
        
        group_red <- select.list(tmpgroup, graphics = FALSE, title = "Select a group for 'Red' on blue-red scale to indicate increased cluster values.", multiple = FALSE)
        
        groups <- c(setdiff(tmpgroup, group_red), group_red)

        p_values <- x@Tocky$PlotClusterPercentage$p_values

        p_values$cluster <- as.numeric(as.character(p_values$cluster))

        
        p_values <- p_values[order(p_values$cluster), ]

        tmp <- calculate_log2fc_color(group_stats, group_red = group_red)
        log2fc_colors <- tmp[['colours']]
        
        
        
        log2fc_cluster_ids <- group_stats$cluster[group_stats$group == groups[1]]
        names(log2fc_colors) <- log2fc_cluster_ids

        if(!is.null(p_adjust)){
            if(p_adjust == 'none'){
                show(summary(p_values))
                significant_clusters <- p_values$cluster[p_values$p_value < 0.05]
            }else{
                significant_clusters <- p_values$cluster[p_values$significance]
            }
            
        }else{
            significant_clusters <- p_values$cluster[p_values$significance]
        }
        
        
        log2fc_colors[!(names(log2fc_colors) %in% significant_clusters)] <- "grey"

        expression <- tmp[['log2fc']]
        V(network)$color <- log2fc_colors[as.character(V(network)$name)]

        plot_network(network, V(network)$color, 'Log2FC Based Visualization', layout)
        par(xpd = TRUE)
        key_rect_x_log2fc <- par()$usr[2] - 0.41
        key_rect_y_log2fc <- par()$usr[4] - 0.13
        key_rect_width_log2fc <- 0.4
        key_rect_height_log2fc <- 0.125
        plot_color_code(expression = expression, x = key_rect_x_log2fc, y = key_rect_y_log2fc, width = key_rect_width_log2fc, height = key_rect_height_log2fc, colour = 'BlueRed', method = 'Expression')
        
        cat(paste("Red colours for clusters enriched in", groups[2]), '\n')
        cat(paste("Blue colours for clusters enriched in", groups[1]), '\n')


    }


    if (select_variable) {
        var <- colnames(data)[grepl(pattern = 'logdata', colnames(data))]
        var <- sub(var, pattern = '.logdata', replacement = '')
        variable <- select.list(var, graphics = TRUE, title = "Data to be plotted on network", multiple = FALSE)
        variable <- paste(variable, '.logdata', sep= '')
    } else {
        variable <- 'Intensity'
    }
    
    df <- x@QCdata$negative_gate_def
    negative_gate <- df[df$variable == variable,]$negative.gate
    expression <- data[,variable]
    expression[expression < negative_gate] = 0
    expression <- tapply(expression, clusters, mean)
    cat(variable)
    cat("\n")
    
    col2 <- cut(expression, breaks = 49)
    lg <- is.na(col2)
    col2 <- as.factor(col2)
    levels(col2) <- colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)
    col2 <- as.character(col2)
    col2[is.na(col2)] <- rgb(0,0,0,alpha = 0.01)
    names(col2) <- names(expression)
    col2 <- col2[names(V(network))]
    
    title <- sub(pattern = '\\.logdata', replacement = '', variable)
    plot_network(network, col2, title, layout)
    
    par(xpd = TRUE)
    key_rect_x <- par()$usr[2] - 0.41
    key_rect_y <- par()$usr[4] - 0.13
    key_rect_width <- 0.4
    key_rect_height <- 0.125
    
    expression <- expression[!is.na(expression)]
    plot_color_code(expression, method = 'Expression', x=key_rect_x, y=key_rect_y, width=key_rect_width, height = key_rect_height, colour = 'Spectral')
    
    return(x)

}

