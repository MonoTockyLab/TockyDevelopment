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

#' Perform Dimensional Reduction using Canonical Correspondence Analysis (CCA) for TockyPrepData.
#' @param x A TockyPrepData after running the function Tocky.
#' @param variables Variables (markers) for CCA analysis. If NULL, variables are to be chosed interactively.
#' @param marker_neg_gate Whether autofluorescence values are all considered to be zero or not. This approach is recommended. Perform DefineNegative for the same variables to be used by TockyCCA in advance.
#' @param select Whether to interactively select markers to be processed. If FALSE, all the log-transformed markers apart from Timer fluorescence (stored in the TockyPrepData) will be used.
#' @return The slot Reduction will contain CCA results.
#' @export
#' @examples
#' \dontrun{
#' x <- TockyCCA(x)
#'}
#' @importFrom vegan cca

TockyCCA <- function(x,  variables = NULL, marker_neg_gate = TRUE,  select = FALSE){

    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
    }
    data <- x@Data[!is.na(x@Data$Angle),]
    used_indices <- which(!is.na(x@Data$Angle))
    
    choices <- colnames(data)
    if(is.null(variables)){
        if(!select){
            var <- x@Transformation$logdata_parameters$logged.channel
            var <- paste(var, 'logdata', sep='.')
        }else{
            var <- select.list(choices, graphics = TRUE, title = "Data to be included for CCA", multiple =TRUE)
        }

    }else{
        var <- variables
    }


    tmpdata  <-  data[,var]
    

    if(marker_neg_gate){
        for(i in var){

            df <- x@QCdata$negative_gate_def
            tp_threshold <- df[df$variable == i,]$negative.gate

            tmp <- tmpdata[,i]
            lg <- tmp <= tp_threshold
            tmp[lg] <- 0
            tmpdata[,i] <- tmp
        }
    }
    
    lg <- rowSums(tmpdata) > 0
    env <- as.matrix(data[lg,c("Angle","Intensity")])
    tmpdata <- tmpdata[lg,]
    used_indices <- used_indices[lg]
    cca_out  <-  cca(tmpdata~ env)
     
    x@Tocky$CCA <- cca_out
    x@Tocky$CCA$variables <- var
    
    x@Tocky$CCA$used_indices <- used_indices
    
    return(x)
}


#' Generate CCA heatmap plots for Tocky data (Timer-Blue vs Timer-Red 2d plots)
#' @param x A TockyPrepData object produced by the function Tocky
#' @param jpeg A logical arguement. If FALSE, it will open a device window in which plots are generated.
#' @param select Whether to interactively select markers to be processed. If FALSE, all the log-transformed markers apart from Timer fluorescence (stored in the TockyPrepData) will be used.
#' @param ncol The number of columns in plot
#' @param nrow The number of rows in plot
#' @export
#' @examples
#' \dontrun{
#' PlotTockyCCA(x)
#'}

PlotTockyCCA  <-  function(x, ncol = 4, nrow  = 3, jpeg = FALSE, select = FALSE){
        if(!inherits(x, "TockyPrepData")){
            stop("Use a TockyPrepData for x. \n")
            
        }
        
        if(is.null(x@Tocky$CCA)){
            stop("Perform TockyCCA before using this function. \n")
            
        }
        
        data <- x@Data
        
        choices <- colnames(data)
        if(!select){
                var <- x@Transformation$logdata_parameters$logged.channel
                var <- paste(var, 'logdata', sep='.')
                
            }else{
                var <- select.list(choices, graphics = TRUE, title = "Data to be included for CCA plot", multiple =TRUE)

        }
            var <- c(var, c("Angle", "Intensity"))
            var <- unique(var)
        
    if(jpeg){
        jpeg(filename = 'PlotTockyCCA.jpeg', width = 480*3, height = 480*3)
        par(mfrow = c(nrow, ncol))
        cex.lab = 2
        
    }else{

        par(mfrow = c(nrow, ncol))
        par(mar = c(4,4,1,1))
        par(cex.axis = 0.6)
        par(mgp = c(1, 0.001, 0))
        par(tck = -0.025)
        cex.lab = 1
        
    }
    
    n <- length(var)
    
    data <- x@Data[rownames(x@Tocky$CCA$CCA$wa),]
    
    for(i in 1:n){
           expression <- c()
           tpexprs  <-  as.data.frame(x@Tocky$CCA$CCA$wa)
           expression <-  data[rownames(x@Tocky$CCA$CCA$wa),var[i]]
           colnames(tpexprs)[1:2] <- paste("CCA", 1:2, sep='')
           
           if(var[i]=='Intensity'){
               expression[expression==0] <- NA
           }
           
           col <- cut(expression, breaks = 49)
           lg <- is.na(col)
           col <- as.factor(col)
           levels(col) <- colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)
           col <- as.character(col)
           col[is.na(col)] <- rgb(0,0,0,alpha = 0.01)

           
           tptitle <-  sub(var[i], pattern='.logdata', replacement = '')
           plot(tpexprs[,1:2], col = col, main = tptitle, pch = 19, cex = 0.2, cex.lab = cex.lab)
       }

    if(jpeg){
        dev.off()
    }
}




#' Generate CCA heatmap plots for Tocky data (Timer-Blue vs Timer-Red 2d plots)
#'
#' @param x A TockyPrepData object produced by the function Tocky
#' @param ncol The number of columns in plot
#' @param nrow The number of rows in plot
#' @param jpeg A logical argument. If FALSE, it will open a device window in which plots are generated.
#' @param select Whether to interactively select markers to be processed. If FALSE, all the log-transformed markers apart from Timer fluorescence (stored in the TockyPrepData) will be used.
#' @param colour Either 'Spectral' or 'BlueRed' for Angle colour key.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#' @param max_cells_displayed The maximum number of cells to display in the plots
#' @param verbose Logical indicating whether to print progress messages and outputs.
#' @return Generates a CCA heatmap plot and optionally saves it as a jpeg file.
#' @examples
#' \dontrun{
#' BiplotCCA(x)
#'}
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette dev.new dev.off rgb
#' @importFrom methods show
#' @import graphics utils
#' @importClassesFrom TockyPrep TockyPrepData
BiplotCCA  <-  function(x, ncol = 2, nrow  = 1, jpeg = FALSE, select = FALSE, colour = 'Spectral', xlim = NULL, ylim = NULL, max_cells_displayed = 30000, verbose = FALSE){
    if(!(inherits(x, "TockyPrepData"))){
        stop("Use a TockyPrepData for x. \n")
    }
    
    if(is.null(x@Tocky$CCA)){
        stop("Perform TockyCCA before using this function. \n")
        
    }
    
    data <- x@Data
    WA_scores <- x@Tocky$CCA$CCA$wa
    
    choices <- colnames(data)
    if(select){
        var <- select.list(choices, graphics = TRUE, title = "Data to be included for CCA plot", multiple =TRUE)
        var <- c(var, c("Angle", "Intensity"))
        var <- unique(var)
    }else{
        var <- c("Angle","Intensity")
    }
    
    
    
    if(jpeg){
        jpeg(filename = file.path('CCAplot.jpeg'), width = 480*3, height = 480*3)
        par(mfrow = c(nrow, ncol))
        cex.lab = 2
        
    }else{
        
        par(mfrow = c(nrow, ncol))
        par(mar = c(4,4,4,4))
        par(cex.axis = 0.8)
        par(mgp = c(2, 0.5, 0))
        par(tck = -0.025)
        par(xpd=TRUE)
        cex.lab = 1
        
    }
    
    data <- x@Data[rownames(WA_scores),]
    
    if(!is.null(max_cells_displayed)){
        max_cells_displayed <- min(max_cells_displayed, nrow(data))
        lg <- sample(1:nrow(data), max_cells_displayed)
        data <- data[lg,]
        WA_scores <- WA_scores[lg,]
    }
    
    for(i in 1:2){
        expression <- c()
        tpexprs  <-  as.data.frame(WA_scores)
        expression <-  data[rownames(WA_scores),var[i]]
        colnames(tpexprs)[1:2] <- paste("CCA", 1:2, sep='')
        
        if(var[i]=='Intensity'){
            expression[expression==0] <- NA
        }
        
        if(var[i]=='Intensity'){
            expression[expression==0] <- NA
            col <- cut(expression, breaks = 49)
            lg <- is.na(col)
            col <- as.factor(col)
            levels(col) <- colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)
            col <- as.character(col)
            col[is.na(col)] <- rgb(0,0,0,alpha = 0.01)
        }
        
        if(var[i]=="Angle"){
            if(colour == 'BlueRed'){
                col <- angle_to_colour(expression)
            }else{
                col <- cut(expression, breaks = 49)
                lg <- is.na(col)
                col <- as.factor(col)
                levels(col) <- colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)
                col <- as.character(col)
                col[is.na(col)] <- rgb(0,0,0,alpha = 0.01)
                
            }
            
        }
                
        tptitle <-  sub(var[i], pattern='.logdata', replacement = '')
        if(!is.null(xlim)&!is.null(ylim)){
            plot(tpexprs[,1:2], col = col, main = tptitle, pch = 19, cex = 0.2, cex.lab = cex.lab, xlim = xlim, ylim = ylim)
        }else{
            plot(tpexprs[,1:2], col = col, main = tptitle, pch = 19, cex = 0.2, cex.lab = cex.lab)
        }
        
        biplot_values <- x@Tocky$CCA$CCA$biplot
        biplot_names <- sub(rownames(biplot_values), pattern='env', replacement='')
        rownames(biplot_values) <- biplot_names
        
        normalized_biplot_vectors <- biplot_values / sqrt(rowSums(biplot_values^2))
        
        plot_min_x <- min(WA_scores[,1])
        plot_max_x <- max(WA_scores[,1])
        plot_min_y <- min(WA_scores[,2])
        plot_max_y <- max(WA_scores[,2])
        
        compute_scaling_factor <- function(vector, plot_min_x, plot_max_x, plot_min_y, plot_max_y) {
            x_component <- vector["CCA1"]
            y_component <- vector["CCA2"]
            
            scaling_x <- if (x_component > 0) (plot_max_x - 0) / x_component else if (x_component < 0) (plot_min_x - 0) / x_component else Inf
            scaling_y <- if (y_component > 0) (plot_max_y - 0) / y_component else if (y_component < 0) (plot_min_y - 0) / y_component else Inf
            
            scaling_factor <- min(abs(scaling_x), abs(scaling_y))
            return(scaling_factor)
        }
        
        num_vectors <- nrow(normalized_biplot_vectors)
        scaling_factors <- sapply(1:num_vectors, function(j) {
            vector <- normalized_biplot_vectors[j,]
            compute_scaling_factor(vector, plot_min_x, plot_max_x, plot_min_y, plot_max_y)
        })
        
        global_scaling_factor <- min(scaling_factors) * 0.8
        
        scaled_biplot_vectors <- normalized_biplot_vectors * global_scaling_factor
        
        arrows(0, 0, scaled_biplot_vectors[, "CCA1"], scaled_biplot_vectors[, "CCA2"], col=1, length=0.1, lwd=1.2)
        text(scaled_biplot_vectors[, c("CCA1", "CCA2")], labels = biplot_names, col = 1, pos = 1, cex=1.5)
        
        tx <- (min(tpexprs[,1]) + max(tpexprs[,1])*2)/3
        ty <- max(tpexprs[,2])*1.1
        height <- max(tpexprs[,2])*0.06
        width <- (max(tpexprs[,1]) - min(tpexprs[,1]))/5
        express <- data[, var[i]]
        if(var[i]=='Intensity'){
            
            if(verbose){
                cat("Intensity \n")
                show(summary(express))
            }

            plot_color_code(expression = express, x= tx, y = ty, width = width, height = height, colour = 'Spectral', method = 'Expression')
            
        }
        if(var[i]=='Angle'){
            
            if(verbose){
                cat("Angle \n")
                show(summary(express))
                
                }
                

            plot_color_code(expression = express, x = tx, y = ty, width = width, height = height,  colour = 'Spectral', method = 'Expression')
            
        }
        
        
    }
    
    if(jpeg){
        dev.off()
    }
}



