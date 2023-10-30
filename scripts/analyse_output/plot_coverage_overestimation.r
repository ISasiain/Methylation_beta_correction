#!/usr/bin/Rscript

# Script Name: plot_coverage_overestimation.r

# Description:
# This script generates a multi-line coverage plot from input data. It accepts an RDS file containing coverage data for different samples and creates a plot that displays the coverage of each sample on the same graph. The resulting plot can be customized with a specified title and saved as a PNG image.

# Required R Packages:
# - optparse: For parsing command-line arguments.
# - ggplot2: For creating and customizing the plot.

# Configuration of Command Line Arguments:
# - The script uses the 'optparse' package to configure and parse command-line arguments.
# - The available options include:
#   - "--input": Path to the input RDS file containing coverage data for multiple samples.
#   - "--title": The title of the plot.
#   - "--output": Name of the output plot file (without the file extension).

# Loading Required R Packages:
# - The script loads the 'optparse' package for command-line argument parsing.
# - It also loads the 'ggplot2' package for plot creation and customization.

# Parsing Command Line Arguments:
# - The 'optparse' package is used to parse the command-line arguments, including the input file, title, and output file name.

# Reading Input File and Configuring Plotting:
# - The script reads the input RDS file and creates a multi-line plot.
# - It iterates through the data for different samples and adds them as separate lines to the plot.
# - The color of each line is automatically chosen from a rainbow of colors.
# - The x-axis represents "1-Purity," and the y-axis represents "Coverage."
# - The plot style is customized with a classic theme and dotted grid lines.

# Saving Plot:
# - After generating the plot, it is saved as a PNG image with the filename provided in the '--output' argument, and the ".png" extension is added.

# Example Usage:
# - Users can run the script from the command line with various options to create a customized multi-line coverage plot.

# Author: Iñaki Sasiain Casado
# Affiliation: Johan Staaf lab


# ==================================
# INSTALLING AND/OR LOADING PACKAGES
# ==================================

if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", quietly = TRUE) }

suppressPackageStartupMessages(library(optparse))

if(!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", quietly = TRUE) }

suppressPackageStartupMessages(library(ggplot2))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

argument_list <- list(

  make_option(c("-i", "--input"), type="character",   
              help="Name (path) of the input file",
              metavar= "[path to data]"),

  make_option(c("-t", "--title"), type="character",
              help="The title of the plot must be entered here",
              metavar= "[Plot's title]"),
  
  make_option(c("-o", "--output"), type="character", 
              help="Name of the output plot (without extension)",
              metavar= "[name of plot]")  

)

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program plots the coverage of different samples together"))


# ===========================================
# READING INPUT FILE AND CONFIGURING PLOTTING
# ===========================================

file_to_plot <- readRDS(arguments$input)
plot_created <- FALSE #Create boolean to check if the plot has been created
my_col <- 0 #Add counter to check the number of coverages that have been added, to get the coluour of the line
all_cols <- rainbow(length(file_to_plot)) #Calculate total number of coverages to plot



for (sample in names(file_to_plot)) {

  #Change color
  my_col <- my_col + 1

  #Tranformr data to plot into dataframe
  sample_to_plot <- data.frame(as.numeric(names(file_to_plot[[sample]][["Coverage_plot"]])), unname(file_to_plot[[sample]][["Coverage_plot"]]))
  colnames(sample_to_plot) <- c("x", "y")

  #Create plot if it has not been created
  if (!plot_created)  {
    
    my_plot <- ggplot(data=sample_to_plot,
                      aes(x=x, y=y, color = all_cols[my_col])) +
                      geom_line(color = all_cols[my_col],linewidth=0.5, alpha=0.5) +
                      xlab("1-Purity") +
                      ylab("Coverage") +
                      theme_classic() +
                      theme(plot.title = element_text(size = 16),
                        axis.text = element_text(size = 14),
                        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
                        panel.grid.minor = element_blank())
                      

    
    plot_created <- TRUE
    
  } else { #Append to plot if it has alreday been created
    
    my_plot <- my_plot + geom_line(data=sample_to_plot,
                                   aes(x=x, y=y), 
                                   color = all_cols[my_col],
                                   linewidth=0.5, alpha=0.5)
    
  }
  
}


# ============
# SAVING PLOT
# ============

ggsave(paste(arguments$output, ".png", sep=""))