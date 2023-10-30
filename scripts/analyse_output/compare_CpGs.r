#!/usr/bin/Rscript

# Script Name: analyze_cpg_context.r

# Description:
# This script analyzes the context of CpGs in different datasets and generates pie charts to visualize the distribution of CpG contexts. It uses input R objects containing CpG vectors and an R object containing the annotated CpG context to create these visualizations.

# Required R Packages:
# - optparse: For parsing command-line arguments.
# - ggplot2: For creating and customizing pie charts.
# - VennDiagram: For generating Venn diagrams.

# Loading Required Packages:
# - The script sets a specific CRAN mirror and installs the required packages if not already installed.
# - It also loads the 'optparse', 'ggplot2', and 'VennDiagram' packages for argument parsing, pie chart creation, and Venn diagram generation, respectively.

# Configuration of Command Line Arguments:
# - The script uses the 'optparse' package to configure and parse command-line arguments.
# - The available options include:
#   - "--cpgs_to_compare": Comma-separated paths to R objects containing CpG vectors for analysis.
#   - "--context_annotation": Path to the R object file containing the annotated CpG context.
#   - "--prefix": Prefix for output file names (optional, default is 'output').

# Parsing Command Line Arguments:
# - The 'optparse' package is used to parse the command-line arguments, including input files and the output file prefix.

# Loading Data:
# - The script reads CpG context annotation data from the specified R object file.
# - It also loads multiple sets of CpG vectors from separate R object files.

# Checking Common CpGs:
# - The script creates Venn diagrams to visualize common CpGs among different datasets and saves them as PNG files.

# Checking Context of the CpGs:
# - The script generates pie charts to show the distribution of CpG contexts for each dataset.
# - It differentiates between promoter, proximal, distal, CGI, shore, ocean, and chromatin accessibility contexts.
# - The generated pie charts are saved as PNG files with appropriate names based on the dataset type.

# Example Usage:
# - Users can run the script from the command line with various options to analyze CpG context and create visualizations.

# Author: Iñaki Sasiain Casado
# Affiliation: Johan Staaf lab @ Lund University / Oncology & Pathology



# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

#Loading optparse to parse command line arguments
if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))

#Loading ggplot2 for plotting
if(!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2") }

suppressPackageStartupMessages(library(ggplot2))

#Loading ggplot2 for plotting
if(!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram") }

suppressPackageStartupMessages(library(VennDiagram))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

#Creating argument list
argument_list <- list(

  make_option(c("-c", "--cpgs_to_compare"), type="character",  
              help="The path of all the R objects containing vectors with CpGs to analyse have to be specified here separated by a comma",
              metavar = "[COMMA_SEPARATED_PATHS]"),

  make_option(c("-a", "--context_annotation"), type="character",  
              help="The path to the RObject file containing the annotated CpG context must be entered here.",
              metavar = "[PATH]"),

  make_option(c("-p", "--prefix"), type="character",  
              help="A prefix for the output plots can be entered here. The default prefix is [%output]",
              default="output",
              metavar = "[PREFIX]")

)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This programm gets purity and original beta Robjects from a R object containing those or more elements. It also allow the user to cerate training anjd validation dataset"))


# ============
# LOADING DATA
# ============

cat("\n\nLoading the data...\n")

# Loading CpG context csv
cpg_context <- readRDS(arguments$context_annotation)

# Loading cpg vectors

#Splitting argument by the comma
path_vec <- strsplit(arguments$cpgs_to_compare, ",")[[1]]

#Creating a list to append cpgs to analyse
cpgs_ls <- list()

#Filling the list
for(file in path_vec) {

    #Generate file name
    name <- strsplit(tail(strsplit(file, "\\/")[[1]], n=1), "_")[[1]][1]

    cpgs_ls[[name]] <- readRDS(file)

}



# ====================
# CHECKING COMMON CPGS
# ====================

cat("\n\nGenerating Venn diagram...\n")

#Creating and saving Venn plot
venn.diagram(cpgs_ls,
             filename=paste(arguments$prefix, ".VennDiagram.png", sep=""))




# ============================
# CHECKING CONTEXT OF THE CpGs
# ============================

# Create dataframe to check the cpgs included
context_df <- as.data.frame(matrix(NA, nrow=length(c(names(cpgs_ls), "ALL")), ncol=length(colnames(cpg_context))))
rownames(context_df) <- c(names(cpgs_ls), "ALL")
colnames(context_df) <- colnames(cpg_context)


cat("\n\nGenerating dataframe...\n")

#Filling the dataframe with the number of CpGs in each context

for (type in rownames(context_df)) {

  #Determine the cpgs whose context has to be checked
  if (type == "ALL") {

    #Using sapply to get a single vector
    cpgs <- rownames(cpg_context)

  } else {

    #Using sapply to get a single vector
    cpgs <- cpgs_ls[[type]]

  }

  # Adding context counts to dataframe
  for (context in colnames(context_df)) {
    context_df[type, context] <- sum(cpg_context[cpgs, context]==1, na.rm=TRUE)
  }

}



cat("\n\nPlotting...\n")

for (type in rownames(context_df)) {

  # Create a dataframe per each type to be plotted 
  type_df <- as.data.frame(cbind(as.numeric(context_df[type,]), as.character(colnames(context_df))))
  colnames(type_df) <- c("count", "context")
  rownames(type_df) <- colnames(context_df)
  
  type_df$"count" <- as.numeric(type_df$"count")



  # Create pie chart for promoter, proximal, and distal
  prom_prox_dis <- type_df[c("promoter","proximal","distal"),]

  gg_prom_prox_dis <- ggplot(data=prom_prox_dis, aes(x="", y=count, fill=context)) +
         geom_bar(stat="identity") +
         theme_void() +
         coord_polar("y", start=0) +
         geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

  ggsave(filename=paste(type, "_prom_prox_dis.CpG_context_piechart.png", sep=""))
  
  # Create pie chart for cgi, shore, and ocean
  cgi_shore_ocean <- type_df[c("cgi","shore","ocean"),]

  gg_cgi_shore_ocean <- ggplot(data=cgi_shore_ocean, aes(x="", y=count, fill=context)) +
         geom_bar(stat="identity") +
         theme_void() +
         coord_polar("y", start=0) +
         geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

  ggsave(filename=paste(type, "_cgi_shore_ocean.CpG_context_piechart.png", sep=""))


  if (type!="ALL") {
    # Create pie chart for chromatin accessibility
    atac <- type_df[c(paste("atac", type, sep=""),paste("nonAtac", type, sep="")),]

    gg_atac <- ggplot(data=atac, aes(x="", y=count, fill=context)) +
          geom_bar(stat="identity") +
          theme_void() +
          coord_polar("y", start=0) +
          geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

    ggsave(filename=paste(type, "_atac.CpG_context_piechart.png", sep=""))

  } else {

    #All BRCA
    atac <- type_df[c("atacBRCA","nonAtacBRCA"),]

    gg_atac <- ggplot(data=atac, aes(x="", y=count, fill=context)) +
          geom_bar(stat="identity") +
          theme_void() +
          coord_polar("y", start=0) +
          geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

    ggsave(filename=paste("ALL_BRCA_atac.CpG_context_piechart.png", sep=""))


    #All LUSC
    atac <- type_df[c("atacLUSC","nonAtacLUSC"),]

    gg_atac <- ggplot(data=atac, aes(x="", y=count, fill=context)) +
          geom_bar(stat="identity") +
          theme_void() +
          coord_polar("y", start=0) +
          geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

    ggsave(filename=paste("ALL_LUSC_atac.CpG_context_piechart.png", sep=""))



    #All LUAC
    atac <- type_df[c("atacLUAD","nonAtacLUAD"),]

    gg_atac <- ggplot(data=atac, aes(x="", y=count, fill=context)) +
          geom_bar(stat="identity") +
          theme_void() +
          coord_polar("y", start=0) +
          geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

    ggsave(filename=paste("ALL_LUAC_atac.CpG_context_piechart.png", sep=""))

  }
}

cat("\n=================\n")
cat ("PROCESS FINISHED")
cat("\n=================\n")