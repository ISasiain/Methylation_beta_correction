#!/usr/bin/Rscript

# Script Name: get_most_variable_cpgs.r

# Description:
# This script allows users to select the most variable CpGs from a reference dataset based on variance or by specifying a number. The selected CpGs are extracted from the dataset to be analyzed.

# Required R Packages:
# - optparse: For parsing command-line arguments.

# Configuration of Command Line Arguments:
# - The script uses the 'optparse' package to configure and parse command-line arguments.
# - The available options include:
#   - "--reference_dataset": Path to the RObject containing the reference data (betas of each CpG).
#   - "--data_to_analyse": Path to the RObject containing the data to analyze (betas of each CpG).
#   - "--filter_by_number": Set to TRUE to select a specific number of the most variable CpGs.
#   - "--number_of_CpGs": Number of most variable CpGs to keep.
#   - "--filter_by_variance": Set to TRUE to select CpGs over a specific variance threshold.
#   - "--variance": Variance threshold for filtering CpGs.
#   - "--include_CpG_list": Set to TRUE to generate a file containing the included CpG names.
#   - "--prefix": Prefix for the output file names.

# Loading Required R Packages:
# - The script loads the 'optparse' package for command-line argument parsing.

# Parsing Command Line Arguments:
# - The 'optparse' package is used to parse the command-line arguments, including file paths and options for filtering CpGs based on variance or a specific number.

# Loading Data:
# - The script reads the reference and data to analyze from RObject files using the 'readRDS' function.

# Getting Most Variable CpGs:
# - The script calculates the variance for all CpGs in the reference dataset.
# - If filtering by a specific number of CpGs is selected, it sorts the CpGs based on variance and selects the top 'number_of_CpGs'.
# - If filtering by variance threshold is selected, it chooses CpGs with variances greater than or equal to the specified 'variance' threshold.

# Generating Output:
# - The selected CpGs are used to extract relevant data from the dataset to be analyzed.
# - The extracted data is saved as an RObject file with a user-defined prefix and file extension "_most_variable_CpGs.RData."
# - If the 'include_CpG_list' option is set to TRUE, a separate RObject file is created with the names of the included CpGs.

# Example Usage:
# - Users can run the script from the command line with various options to select the most variable CpGs based on their preferences.

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


# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

#Creating argument list
argument_list <- list(

  make_option(c("-r", "--reference_dataset"), type="character",  
              help="The path to the RObject containing the reference data (betas of each CpG) must be entered here.",
              metavar = "[PATH_TO_REFERNCE]"),

  make_option(c("-a", "--data_to_analyse"), type="character",  
              help="The path to the RObject containing the data to analyse (betas of each CpG) must be entered here.",
              metavar = "[PATH_TO_DATA_TO_ANALYSE]"),

  make_option(c("-N", "--filter_by_number"), type="numeric", default=FALSE,
              help="This argument must be set to TRUE if the user wants to get a certain number of the most variant CpG",
              metavar="[TRUE/FALSE]"),

  make_option(c("-n", "--number_of_CpGs"), type="numeric",  
              help="The number of the most variable CpGs to keep must be entered here.",
              metavar = "[CpG_NUM]"),

  make_option(c("-V", "--filter_by_variance"), type="numeric", default=FALSE,
              help="This argument must be set to TRUE if the user wants to get teh CpGs over a specific variance threshold",
              metavar="[TRUE/FALSE]"),

  make_option(c("-v", "--variance"), type="numeric",  
              help="The variance threshold to filter the CpGs must be entered here CpGs to keep must be entered here.",
              metavar = "[Var]"),

  make_option(c("-l", "--include_CpG_list"), type="logical",  
              help="This argument must be set to TRUE to generate a file containing the included CpG names as an output. The default argument is [%default]",
              default=FALSE,
              metavar = "[TRUE/FALSE]"),

  make_option(c("-p", "--prefix"), type="character",  
              help="The prefix of the output files can be selected here. The default argument is [%default]",
              default="output",
              metavar = "[PREFIX]")
)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This programm gets purity and original beta Robjects from a R object containing those or more elements. It also allow the user to cerate training anjd validation dataset"))


# ===========
# LODING DATA
# ===========

# Reading the reference betas
reference <- readRDS(arguments$reference_dataset)

# Reading the data to analyse
to_analyse <- readRDS(arguments$data_to_analyse)


# =================================
#    GETTING MOST VARIABLE CPGs
# =================================

if (arguments$filter_by_number) {

# Detemine the variance of all the rows (CpGs)
 cpgs_variance <- apply(
                        reference,
                        MARGIN=1,
                        FUN=var
                   )

# Create a vector to sort the rows. Get only the number of rows containing the CpGs that want to be included
sorting_vec <- order(cpgs_variance, decreasing=TRUE)[1:arguments$number_of_CpGs]

# Sorting the betas dataframe and getting only the CpGs to include
cpgs <- rownames(reference[sorting_vec,])

} else if (arguments$filter_by_variance) {

# Detemine the variance of all the rows (CpGs)
cpgs_variance <- apply(
                        reference,
                        MARGIN=1,
                        FUN=var
                   )

# Get the CpGs over the variance threshold

cpgs <- rownames(refernce)[cpgs_variance >= arguments$variance]
}

# ===========================
#     GENERATING OUTPUT
# ===========================

to_analyse <- to_analyse[cpgs,]

saveRDS(to_analyse, file=paste(arguments$prefix, "_most_variable_CpGs.RData", sep=""))

if (arguments$include_CpG_list) {

    saveRDS(cpgs, file=paste(arguments$prefix, "_CpG_vector.RData", sep="")) 

}