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

#' Plot Cluster Percentage by Bar plot
#'
#' @param x A TockyPrepData after running the function ClusteringPCA.
#' @param p_adjust_method A method for p-value adjustment in statistical tests.
#' @param colours (optional) A vector specifying colours for different groups in the plot. For example, colours = c("purple", "black").
#' @param verbose Logical indicating whether to print progress messages and outputs.
#' @return The TockyPrepData with added statistics and plots for significant clusters.
#' @export
#' @examples
#' \dontrun{
#' x <- PlotClusterPercentage(x)
#'}
#' @import ggplot2 stats
#' @importFrom methods show
#' @importClassesFrom TockyPrep TockyPrepData
PlotClusterPercentage <- function(x, p_adjust_method = 'fdr', colours = NULL, verbose = TRUE){
    if(!(inherits(x, "TockyPrepData"))){
        stop("Use a TockyPrepData. \n")
    }

    sampledef <- x@sampledef$sampledef
    if(is.null(sampledef)){
        stop("Perform SampleDef. \n")
    }

    
    if(is.null(x@Tocky$PCAclusters)){
        stop("Perform ClusteringPCA before using this function. \n")
    }
    
    x <- cluster_cell_num(x)
    df_res <- x@Tocky[['cluster_cellnumbers']]
    df_percentage <- 100 * apply(df_res, 2, function(x) x/sum(x))
    
    df_long <- as.data.frame(df_percentage)
    df_long$cluster <- rownames(df_long)
    
    df_long <- reshape(df_long,
    direction = "long",
    varying = list(names(df_long)[-which(names(df_long) == "cluster")]),
    v.names = "percentage",
    idvar = "cluster",
    timevar = "file",
    times = names(df_long)[-which(names(df_long) == "cluster")])
    df_long <- df_long[order(df_long$cluster), ]
    
    df_long <- merge(df_long, sampledef, by = "file", all.x = TRUE)
    
    df_group_averages <- aggregate(percentage ~ group + cluster, data = df_long, function(x) {
        c(average = mean(x), sd = sd(x))
    })
    df_group_averages <- do.call(data.frame, df_group_averages)
    names(df_group_averages)[-(1:2)] <- c("average_percentage", "sd_percentage")
    df_group_averages$average_percentage <- sapply(df_group_averages$average_percentage, '[', 1)
    df_group_averages$sd_percentage <- sapply(df_group_averages$sd_percentage, '[', 1)
    
    p_values_list <- split(df_long, df_long$cluster)
    p_values_data <- lapply(p_values_list, function(data) {
        wilcox_test <- wilcox.test(percentage ~ group, data = data)
        data.frame(cluster = unique(data$cluster), p_value = wilcox_test$p.value)
    })
    
    p_values <- do.call("rbind", p_values_data)
    
    
    p_values$p_adjusted <- p.adjust(p_values$p_value, method = p_adjust_method)
    p_values$significance <- p_values$p_adjusted < 0.05
    df_group_averages <- merge(df_group_averages, p_values, by = "cluster")
    df_group_averages$cluster <- factor(df_group_averages$cluster, levels = as.character(1:max(as.numeric(df_group_averages$cluster))))
    df_group_averages <- df_group_averages[order(as.numeric(df_group_averages$cluster)), ]

    
    significant_clusters <- df_group_averages[df_group_averages$significance == TRUE, ]
    
    
    if(nrow(significant_clusters) > 0) {
        significant_clusters <- aggregate((average_percentage + sd_percentage) ~ cluster,
                                          data = significant_clusters,
                                          max)
        names(significant_clusters)[2] <- "y_position"
    } else {
        significant_clusters <- data.frame(cluster = character(0), y_position = numeric(0))
    }

    
    if(verbose){
        show(p_values)
        
    }
    p1 <- ggplot(df_group_averages, aes(x = !!sym("cluster"), y = !!sym("average_percentage"), fill = !!sym("group"))) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = !!sym("average_percentage") - !!sym("sd_percentage"),
    ymax = !!sym("average_percentage") + !!sym("sd_percentage")),
    position = position_dodge(width = 0.9), width = 0.25) +
    theme_minimal() +
    labs(x = "Cluster", y = "Percentage (%)", title = "Cluster Percentage") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

    
    p1 <-    p1 + geom_text(data = significant_clusters,
    aes(x = !!sym("cluster"), y = !!sym("y_position"), label = "*"),
    position = position_dodge(width = 0.45),
    vjust = -0.5,
    size = 6,
    inherit.aes = FALSE)
    
    if(!is.null(colours)){
        
        p1 <- p1 + scale_fill_manual(values = colours)
    }
    
    print(p1)
    
    out <- list(cluster_percentage = df_long,
    group_statistics = df_group_averages,
    p_values = p_values)
    
    x@Tocky[['PlotClusterPercentage']] <- out

    return(x)
}




#' Retrieve Cluster Percentage Data and Stats
#'
#' This function retrieves the cluster percentage data from a TockyPrepData that has already been processed with the PlotClusterPercentage function. It can display the statistics in the Terminal or write them to CSV files.
#'
#' @param x A TockyPrepData that has been processed with the PlotClusterPercentage function.
#' @param writeResults A logical value. If TRUE, two files will be generated containing group statistics and p-values, respectively. If FALSE, these statistical results are displayed in the Terminal.
#' @param filename (optional) Base name for the output files when writeResults is TRUE.
#' @return The same TockyPrepData passed as input, for consistency in function design, though the function primarily focuses on data retrieval and display or file writing.
#' @export
#' @examples
#' \dontrun{
#'   GetStatsClusterPercentage(x)
#'}
#' @importFrom utils write.table
#' @importFrom methods show
#' @importClassesFrom TockyPrep TockyPrepData
GetStatsClusterPercentage <- function(x, writeResults = FALSE, filename = "cluster_percentage"){
    if(!(inherits(x, "TockyPrepData"))){
        stop("Use a TockyPrepData. \n")
    }

    if(is.null(x@Tocky$PlotClusterPercentage)){
        stop("Apply PlotClusterPercentage. \n")
    }
    
    df <- x@Tocky$PlotClusterPercentage
    
    if(!writeResults){
        show(df$group_statistics)
        show(df$p_values)
    } else {
        write.table(df$group_statistics, file = paste0(filename, '_group_statistics.csv'), sep=',', quote = FALSE, row.names = FALSE, col.names = TRUE)
        write.table(df$p_values, file = paste0(filename, '_pvalues.csv'), sep=',', quote = FALSE, row.names = FALSE, col.names = TRUE)
    }

    return(x)
}




