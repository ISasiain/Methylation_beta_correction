#!/usr/bin/Rscript

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
uppressPackageStartupMessages(library("ComplexHeatmap"))

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

  make_option(c("-p", "--output_prefix"), type="character", default="Heatmap"
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
png("original_and_corrected_heatmap.png", width=800, height=800)

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
        top_annotation=a,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2")

dev.off()
