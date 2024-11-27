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

#' Perform Principal Component Analysis for TockyPrepData.
#' @param x A TockyPrepData after running the function Tocky.
#' @param variables Variables (markers) for PCA analysis. If NULL, variables are to be chosed interactively.
#' @param marker_neg_gate Whether autofluorescence values are all considered to be zero. This approach is recommended and requires DefineNegative.
#' @param cleaning Whether data cleaning is performed to remove cells with a zero value for a marker. This is not recommended unless you have a reason that you cannot use DefineNegative to collapse negative data.
#' @param select Whether to interactively select markers to be processed. If FALSE, all the log-transformed markers apart from Timer fluorescence (stored in the TockyPrepData) will be used.
#' @param Timer Whether Timer data (normalised or Angle/Intensity) are used.
#' @export
#' @examples
#' \dontrun{
#' x <- TockyPCA(x, variables = NULL, cleaning = FALSE, marker_neg_gate = TRUE)
#'}
#' @importFrom stats prcomp
TockyPCA <- function(x, variables = NULL, marker_neg_gate = TRUE, cleaning = FALSE, select = TRUE, Timer = FALSE){

    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
    }
    data <- x@Data
    
    if(Timer){
        choices <- c(colnames(data)[grepl(pattern = 'logdata', colnames(data))|grepl(pattern = 'normalised', colnames(data))], 'Angle', 'Intensity')
    }else{
        choices <- c(colnames(data)[grepl(pattern = 'logdata', colnames(data))])
    }


    if(is.null(variables)){
        if(!select){
            var <- x@Transformation$logdata_parameters$logged.channel
            var <- paste(var, 'logdata', sep='.')
        }else{
            var <- select.list(choices, graphics = TRUE, title = "Data to be included for PCA", multiple =TRUE)
        }

    }else{
        var <- variables
    }
    
    if(Timer){
        var2 <- setdiff(var, c(colnames(data)[grepl(pattern = 'normalised', colnames(data))], 'Angle', 'Intensity'))
    }else{
        var2 <- var
    }
    
    if(any(var %in% c("Angle", "Intensity"))){
        data <- data[!is.na(data$Angle),]
    }

    if(marker_neg_gate){
        if(!is.null(x@QCdata$negative_gate_def)){
            neg_df_var <- x@QCdata$negative_gate_def[,"variable"]
            logic <- all(var2 %in% neg_df_var)
            if(!logic){
                stop("Perform DefineNegative for the same set of variables. \n")
            }
        }else{
            stop("Perform DefineNegative before runnng this function. \n")

        }
    }
    tmpdata  <-  data[,var]

    if(marker_neg_gate){
        

        
        for(i in var2){
            df <- x@QCdata$negative_gate_def
            lgdf <- df[,"variable"] == i
            tp_threshold <- df[lgdf, "negative.gate"]
            tmp <- tmpdata[,i]
            tmp[tmp <= tp_threshold] = 0
            tmpdata[,i] <- tmp
        }
    }
    pca_out  <-  prcomp(tmpdata)
    x@Tocky$PCA <- pca_out
    x@Tocky$PCA$variables <- var
    
    return(x)
}



#' Produce PCA plot
#' @param x A TockyPrepData after running the function Tocky.
#' @param filename A character string for file name
#' @param jpeg A logical arguement. If FALSE, it will open a device window in which plots are generated.
#' @param cluster A logical arguement. If TRUE, clusters are coloured in output plots.
#' @return A TockyPrepData simply to avoid null return.
#' @export
#' @examples
#' \dontrun{
#' PlotTockyPCA(x)
#'}
#' @importClassesFrom TockyPrep TockyPrepData
PlotTockyPCA <- function(x, jpeg = FALSE, cluster = FALSE, filename = NULL){

    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
    }
    
    if(is.null(x@Tocky$PCA)){
        stop("Perform TockyPCA before running this function. \n")
        
    }
    data <- x@Data

    pca_out  <-  x@Tocky$PCA
    var <- x@Tocky$PCA$variables
    if(cluster){
        clusters <- x@Clustering$PCAclusters
        tempcol <- as.factor(clusters); levels(tempcol) <- 1:length(levels(tempcol))
        tempcol <- as.vector(tempcol); tempcol <- as.numeric(tempcol)
        
        col <- color_palette(unique(tempcol))
        tempcol <- as.factor(tempcol)
        levels(tempcol) <- col
        tempcol <- as.character(tempcol)
    }
    if(length(var) > 2){
        n <- round(length(var)^0.5, digits = 0)+1
        if(jpeg){
            if(is.null(filename)){
                filename <- 'PCA.jpeg'
            }
            jpeg(filename = filename, width = 480*3, height = 480*3)
        }
        
            par(mfrow = c(n,n))

            par(mar = c(4,3,1,1))
            par(cex.axis = 0.6)
            par(mgp = c(1.2, 0.2, 0))
            par(tck = -0.025)
            cex.lab = 1
        for(j in 1:(length(var)-1)){
            tmp <- pca_out$x[,j:(j+1)]
            if(nrow(tmp) > 30000){
                tmp <- tmp[sample(1:nrow(tmp), 30000), ]
            }
            
            if(cluster){

                plot(tmp, pch ='.', col = tempcol)
                
            }else{
                plot(tmp, pch ='.', col = 4)
            }
            
        }
        if(jpeg){
            dev.off()
        }

        
    }
    
    scree <- pca_out$sdev
    names(scree) <- paste("PC", 1:length(scree), sep='')
    
    pdf(file = 'scree_plot.pdf')
    barplot(scree)
    dev.off()
    
    return(invisible(x))
}


#' Set a gate to define negative (and positive) for each marker expression
#' @param x A TockyPrepData.
#' @param reduction Choose whether to use PCA or CCA as a reduction method.
#' @return A TockyPrepData (unchanged) for safety.
#' @export
#' @examples
#' \dontrun{
#' PlotDimRedLoading(x)
#'}
#' @importClassesFrom TockyPrep TockyPrepData

PlotDimRedLoading <- function(x, reduction = 'PCA'){
    
    if(reduction == 'CCA'){
        tmp <- x@Tocky$CCA$CCA$wa[,1:2]
        xlim <- c(min(tmp[,1]), max(tmp[,1]));ylim <- c(min(tmp[,2]), max(tmp[,2]));
        loading<- x@Tocky$CCA$CCA$v[,1:2]
        plot(loading[,1], loading[,2], pch = 19, cex = 0.2, xlab="CCA1", ylab="CCA2", main="CCA Variables", xlim = xlim, ylim = ylim)
        tmptxt <- sub(rownames(loading), pattern = '.logdata', replacement = '')
        text(loading[,1], loading[,2], labels=tmptxt, cex=1, pos = 1)
        abline(v= 0, h= 0 , col = 8)
        biplot <- x@Tocky$CCA$CCA$biplot[,1:2]*5
        arrows(0, 0, biplot[,1], biplot[,2],col="blue", length=0.1)
        tmptxt2 <- sub(rownames(biplot), pattern = 'env', replacement = '')
        text(biplot[,1], biplot[,2], labels=tmptxt2, cex=1, pos=3, col = 4)
                              
             
    }
    if(reduction == 'PCA'){
        loading <- x@Tocky$PCA$rotation[,1:2]
        # Plot
        plot(loading[,1], loading[,2], type="n",
             xlab="PC1", ylab="PC2",
             main="PCA Loadings")
        arrows(0, 0, loading[,1], loading[,2],col="blue", length=0.1)
        tmptxt <- sub(rownames(loading), pattern = '.logdata', replacement = '')
                    
        abline(v= 0, h= 0 , col = 8)


        text(loading[,1], loading[,2], labels=tmptxt,
                         cex=1, pos=3)
    }
    
     return(invisible(x))

}





#' Generate PCA heatmap plots for Tocky data (Timer-Blue vs Timer-Red 2d plots)
#'
#' @param x A TockyPrepData object produced by the function Tocky
#' @param ncol The number of columns in plot
#' @param nrow The number of rows in plot
#' @param jpeg A logical argument. If FALSE, it will open a device window in which plots are generated.
#' @param select Whether to interactively select markers to be processed. If FALSE, all the log-transformed markers apart from Timer fluorescence (stored in the TockyPrepData) will be used.
#' @param biplot_scaling A number for multiplying biplot values for visibility. Default is 3.
#' @param colour Either 'Spectral' or 'BlueRed' for Angle colour key.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#' @return Generates a CCA heatmap plot and optionally saves it as a jpeg file.
#' @examples
#' \dontrun{
#' PlotPCAHeatmap(x)
#'}
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette dev.new dev.off rgb
#' @import graphics utils
#' @importClassesFrom TockyPrep TockyPrepData
PlotPCAHeatmap  <-  function(x, ncol = 2, nrow  = 1, jpeg = FALSE, select = FALSE, biplot_scaling = 3, colour = 'Spectral', xlim = NULL, ylim = NULL){
    if(!(inherits(x, "TockyPrepData"))){
        stop("Use a TockyPrepData for x. \n")
    }
    
    if(is.null(x@Tocky$PCA)){
        stop("Perform TockyPCA before using this function. \n")
        
    }
    
    data <- x@Data
    PCAscores <- x@Tocky$PCA$x
    
    choices <- colnames(data)
    if(select){
        var <- select.list(choices, graphics = TRUE, title = "Data to be included for PCA plot", multiple =TRUE)
        var <- c(var, c("Angle", "Intensity"))
        var <- unique(var)
    }else{
        var <- c("Angle","Intensity")
    }
    
    
    
    if(jpeg){
        jpeg(filename = file.path('PCAplot.jpeg'), width = 480*3, height = 480*3)
        par(mfrow = c(nrow, ncol))
        cex.lab = 2
        
    }else{
        dev.new()
        par(mfrow = c(nrow, ncol))
        par(mar = c(4,4,1,1))
        par(cex.axis = 0.6)
        par(mgp = c(1, 0.001, 0))
        par(tck = -0.025)
        cex.lab = 1
        
    }
    
    for(i in 1:2){
        expression <- c()
        tpexprs  <-  as.data.frame(PCAscores)
        expression <-  data[,var[i]]
        colnames(tpexprs)[1:2] <- paste("PC", 1:2, sep='')
        
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
      
        
        tx <- (min(tpexprs[,1]) + max(tpexprs[,1]))/2
        ty <- max(tpexprs[,2])*0.95
        height <- max(tpexprs[,2])*0.05
        width <- (max(tpexprs[,1]) - min(tpexprs[,1]))/5
        
        if(var[i]=='Intensity'){
            plot_color_code(expression = expression, x= tx, y = ty, width = width, height = height, colour = 'Spectral', method = 'Expression')
            
        }
        if(var[i]=='Angle'){
            plot_color_code(x = tx, y = ty, width = width, height = height,  colour = 'Spectral', method = 'Angle')
            
        }
        
        
    }
    
    if(jpeg){
        dev.off()
    }
}



