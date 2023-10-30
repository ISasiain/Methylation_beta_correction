#!/usr/bin/Rscript

# Script Name: heatmap_script.r

# Description:
# This script generates a heatmap to visualize the methylation status of CpGs in different samples. It uses three sets of beta values (original, cancer-corrected, and microenvironment-corrected) and an annotation file to create the heatmap. The generated heatmap provides insights into methylation changes in various samples.

# Required R Packages:
# - ComplexHeatmap: For creating and customizing the heatmap.
# - optparse: For parsing command-line arguments.

# Loading Required Packages:
# - The script sets a specific CRAN mirror and installs the 'ComplexHeatmap' package if not already installed.
# - It also loads the 'ComplexHeatmap' and 'optparse' packages for heatmap generation and command-line argument parsing, respectively.

# Configuration of Command Line Arguments:
# - The script uses the 'optparse' package to configure and parse command-line arguments.
# - The available options include:
#   - "--original_betas": Path to the RObject file containing the original uncorrected beta values.
#   - "--cancer_corrected_betas": Path to the RObject file containing the corrected tumor beta values.
#   - "--microenvironment_corrected_betas": Path to the RObject file containing the corrected microenvironment beta values.
#   - "--annotation_file": Path to the annotation CSV file with the methylation status of each sample.
#   - "--output_prefix": Prefix for the output PNG heatmap file.

# Loading Data:
# - The script reads the beta values and annotation data from the specified files using the 'readRDS' and 'read.table' functions, respectively.

# Creating the Heatmap:
# - The script generates a heatmap with three panels, one for each set of beta values (original, cancer-corrected, microenvironment-corrected).
# - The heatmap is customized with various options, including colors, labels, and clustering parameters.
# - The resulting heatmap is saved as a PNG file with the user-defined prefix.

# Example Usage:
# - Users can run the script from the command line with various options to create a customized methylation heatmap.

# Author: Iñaki Sasiain Casado
# Affiliation: Johan Staaf lab @ Lund University / Oncology & Pathology


# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

# Installing pachage to create heat-maps
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

#Loading R package
suppressPackageStartupMessages(library("ComplexHeatmap"))

if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))



# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

argument_list <- list(

  make_option(c("-o", "--original_betas"), type="character", 
              help="The path to the RObject file containing the original uncorrected beta values must be entered here.",
              metavar = "[path]"),

  make_option(c("-c", "--cancer_corrected_betas"), type="character",
              help="The path to the RObject file containing the corrected tumor beta values must be entered here.",
              metavar="[path]"),

  make_option(c("-m", "--microenvironment_corrected_betas"), type="character",
              help="The path to the RObject file containing the corrected microenvironment beta values must be entered here.",
              metavar="[path]"),

  make_option(c("-a", "--annotation_file"), type="character",
              help="The path to the annotation csv file with the methylation status of each sample must be entered here.",
              metavar="[path]"),

  make_option(c("-p", "--output_prefix"), type="character", default="Heatmap",
              help="The prefix pf the png file in which the created heatmap will be saved can be entered here",
              metavar="[prefix]")         

)

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program plots the results of the methylation beta value correction based on the beta correction for samples without known purity values."))


# ================
# LOADING THE DATA
# ================

# Loading beta values
original_values <- readRDS(arguments$original_betas)  #Uncorrected beta values
corrected_values <- readRDS(arguments$cancer_corrected_betas)  #Corrected beta values
microenvironment_values <- readRDS(arguments$microenvironment_corrected_betas)


# Loading annotation
annotation <- read.table(arguments$annotation_file) #Annotation file
annotation <- annotation[colnames(original_values),] #Removing non used samples



# ====================
# CREATING THE HEATMAP
# ====================


# Plotting heat map
png(paste(arguments$output_prefix, ".png", sep=""), width=400, height=300)

a <- HeatmapAnnotation(
        Pyroseq=annotation,
        col=list(Pyroseq=c("Hypermethylated"="black", "BRCA1_null"="grey"))
)

Heatmap(original_values, 
        column_title="ORIG. β",
        row_title="CpGs",
        name="β value", 
        show_row_dend = FALSE,
        show_column_names = FALSE,
        show_row_name = FALSE,
        top_annotation=a,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2") +

Heatmap(corrected_values, 
        column_title="TUMOR β",
        row_title="CpGs",
        name="β value",
        show_row_dend = FALSE,
        show_column_names = FALSE, 
        show_row_name = FALSE,
        top_annotation=a,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2") +

Heatmap(microenvironment_values, 
        column_title="MICROENV. β",
        row_title="CpGs",
        name="β value",
        show_row_dend = FALSE,
        show_column_names = FALSE,
        show_row_name = FALSE,
        top_annotation=a,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2")

dev.off()
