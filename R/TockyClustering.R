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


#' Clustering Cells Using Dimensional Reduction by TockyPCA.
#' @param x A TockyPrepData after running the function TockyPCA.
#' @param choose_dimension Whether the number of dimension is specified by PCA plot. Default is FALSE.
#' @param num_centre Whether the number of clusters is specified. Default is FALSE and the number of variables will be used.
#' @param max_cells_displayed The number of cells displayed in plots for interatcive PCA plots.
#' @return The slot Reduction will contain the new slot PCAclusters, which includes kmeans clustering result
#' @export
#' @examples
#' \dontrun{
#' x <- TockyClustering(x)
#'}
#' @importFrom bigmemory as.big.matrix
#' @importFrom biganalytics bigkmeans

TockyClustering <- function(x, choose_dimension = FALSE, max_cells_displayed = 30000, num_centre = FALSE){

    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
    }
    if(is.null(x@Data$Locus)){
        stop("Perform TockyLocus before using this function. \n")
    }
    
    
    if(is.null(x@Tocky$PCA)){
        stop("Perform TockyPCA before using this function. \n")
    }
    data <- x@Data
    pca_out <- x@Tocky$PCA$x

    t_file <- data$file
    
    if(choose_dimension|num_centre){
            if(nrow(pca_out)> max_cells_displayed){
                tlg <- sample(1:nrow(pca_out), max_cells_displayed)
                t_pca_out <- pca_out[tlg, ]
                t_file <- t_file[tlg]
            }else{
                t_pca_out <- pca_out
                
            }
            plot(t_pca_out[,1:2], pch= '.', col = rgb(0,0,0,alpha=0.1), main ='PCA')
            
    }
    
    if(num_centre){
        n <- readline(prompt ="How many centres are used for k-means clustering? \n")
    }else{
        
        n <- nrow(x@Tocky$PCA$rotation)

        }
        
    if(choose_dimension){

        dim <- readline(prompt ="How many dimensions are used for k-means clustering? \n")
    }else{
        dim <- 3
    }
    n <- as.numeric(n); dim <- as.numeric(dim)
    
    big_matrix <- as.big.matrix(pca_out[,1:dim])

    clustering_result <-  bigkmeans(big_matrix, centers = n)

    cat(paste(n, "centres has been used for k-means. \n"))
    clusters <- clustering_result$cluster; names(clusters) = as.character(clusters)
     
     
    x@Tocky[['PCAclusters']] <- as.factor(clusters)
    return(x)
    
}


#' Plot Tocky Clusters
#' @param x A TockyPrepData after running the function TockyClustering
#' This function will generate PCA plots with cluster and Angle data
#' @param jpeg Whether to out a jpeg file. The default is pdf = FALSE, by which a jpeg file is produced.
#' @param filename A character string for file name
#' @param max_cells_displayed The number of cells displayed in plots.

#' @return The slot Reduction will contain the new slot Tocky_clusters, which includes kmeans clustering result
#' @export
#' @examples
#' \dontrun{
#' PlotTockyClustering(x)
#'}
#' @importFrom TockyLocus Locus_to_colour TockyLocusLegend
#' @importFrom grDevices pdf cm.colors
#' @importClassesFrom TockyPrep TockyPrepData
PlotTockyClustering <- function(x,  jpeg = FALSE, max_cells_displayed = 30000, filename = NULL){
    
    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
    }
    
    if(is.null(x@Tocky[['PCAclusters']])){
        stop("Perform TockyClustering before using this function. \n")
    }
    
    if(jpeg){
        if(is.null(filename)){
            filename <- 'PlotTockyClustering'
        }
        filename <- paste(filename, 'jpeg', sep='.')
        jpeg(filename = filename, width = 480*3, height = 480*3)
        par(mfrow = c(3,3))
        par(mar = c(5,5,5,5))
    }else{
        par(mfrow= c(3,3))
        par(mar = c(4,4,3,3))
    }
    locus_colour <- Locus_to_colour(x@Data$Angle)
    clusters <- x@Tocky[['PCAclusters']]
    pca_out <- x@Tocky$PCA$x
    
    if(nrow(pca_out)> max_cells_displayed){
        select_logic <- sample(1:nrow(pca_out), max_cells_displayed)
        pca_out  <-  pca_out[select_logic,]
        clusters <- clusters[select_logic]
        locus_colour <- locus_colour[select_logic]
    }
    
    tempcol <- as.factor(clusters); levels(tempcol) <- 1:length(levels(tempcol))
    tempcol <- as.vector(tempcol); tempcol <- as.numeric(tempcol)
    
    col <- color_palette(unique(tempcol))
    tempcol <- as.factor(tempcol)
    levels(tempcol) <- col
    tempcol <- as.character(tempcol)
    anglena_lg <- is.na(x@Data$Angle)
    tmp <- tapply(x@Data$Angle[!anglena_lg], clusters[!anglena_lg], mean)
    tmp <- Locus_to_colour(tmp)
    palette <- as.factor(clusters); levels(palette) <- tmp
    palette <- as.character(palette)
    
    if(!jpeg){
        for(i in 1:2){
            plot(pca_out[,c(i, i+1)], col=locus_colour, pch='.', main = 'Tocky Locus of Each Cell')
            plot(pca_out[,c(i, i+1)], col=tempcol, pch='.', main = 'Clusters (Colour-Coded)')

            mean_x <- tapply(pca_out[,i], clusters, mean)
            mean_y <- tapply(pca_out[,i+1], clusters, mean)

            plot(pca_out[,c(i, i+1)], col = rgb(0,0,0,alpha=0.01), pch='.', main = 'Clusters (ID in Barycentre)')
            text(mean_x, mean_y, labels = names(mean_x), col = 2, cex= 1.5)
        }
        TockyLocusLegend(cex = 2)
        plot(x = c(-1,1), y = c(-1,1),  ann = FALSE, xaxt = 'n', yaxt= 'n',  bty = "n",col = 0)
        text(-0.75, 1, labels = "Clusters", cex = 1)
        clusters <- as.factor(clusters)
        collevel <- color_palette(as.numeric(levels(clusters)))
        df <- data.frame(Cluster = levels(clusters))
        legend(-0.9,0.9, legend = levels(clusters), col = collevel, pch = 19, cex = 0.5)

    }
    
    
    if(jpeg){
        for(i in 1:2){
            plot(pca_out[,c(i, i+1)], col=locus_colour, pch='.', main = 'Tocky Locus (per each cell)')

            plot(pca_out[,c(i, i+1)], col=tempcol, pch='.', main = 'Clusters')
            
            lg <- is.na(x@Data$Angle)
            tmp <- tapply(x@Data$Angle[!lg], clusters[!lg], mean)
            tmp <- Locus_to_colour(tmp)
            palette <- as.factor(clusters); levels(palette) <- tmp
            palette <- as.character(palette)

            mean_x <- tapply(pca_out[,i], clusters, mean)
            mean_y <- tapply(pca_out[,i+1], clusters, mean)
            plot(pca_out[,c(i, i+1)], col = rgb(0,0,0,alpha=0.01), pch='.', main = 'Clusters')
            text(mean_x, mean_y, labels = names(mean_x), col = 2, cex= 1.5)
            
        }
        TockyLocusLegend(cex = 2)

        plot(x = c(-1,1), y = c(-1,1),  ann = FALSE, xaxt = 'n', yaxt= 'n',  bty = "n",col = 0)
        text(-0.75, 1, labels = "Clusters", cex = 2)
        clusters <- as.factor(clusters)
        collevel <- color_palette(as.numeric(levels(clusters)))
        df <- data.frame(Cluster = levels(clusters))
        legend(-0.9,0.9, legend = levels(clusters), col = collevel, pch = 19, cex = 1)
        
        dev.off()
        
    }
    
    return(invisible(x))
}



