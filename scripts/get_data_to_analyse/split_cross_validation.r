#!/usr/bin/Rscript

# Script Name: split_cross_validation.r

# Description:
# This script processes and splits a dataset of beta values and purities for cross-validation. 
# It can work with data where beta values and purities are in separate R objects or in the same object.
# It allows users to filter CpGs from sexual chromosomes, select the number of most variable CpGs to include, and perform cross-validation.

# Required R Packages:
# - optparse: Used for parsing command-line arguments.
# - caret: Required for data splitting and cross-validation.
# - Other standard packages are used for data manipulation and analysis.

# Loading Required Packages:
# The script sets a specific CRAN mirror and installs the required packages if not already installed.
# It also loads the 'optparse' and 'caret' packages for argument parsing and data splitting.

# Configuration of Command Line Arguments:
# The script uses the 'optparse' package to configure and parse command-line arguments.
# The available options include paths to data files, argument flags for various data processing steps, and more.

# Parsing Command Line Arguments:
# The script parses the command-line arguments to configure the data processing tasks, such as data file paths, cross-validation settings, and filtering options.

# Processing Command Line Arguments:
# The provided arguments are processed, and relevant data is loaded. It determines whether beta values and purities are in separate objects or the same object.

# Data Loading and Preprocessing:
# Depending on the configuration, the script loads beta values and purities from the provided R objects. It adapts column names for beta values.

# Filtering CpGs from Sexual Chromosomes:
# If the user opts to filter CpGs from sexual chromosomes, the script loads an annotation file, matches it with the beta values, and removes relevant rows.

# Creating Groups for Cross-Validation:
# The script uses the 'caret' package to split the data into groups for cross-validation. The number of groups can be customized.

# Getting Most Variable CpGs:
# If specified, the script selects the most variable CpGs based on variance for inclusion in the output data.

# Creating Output Files for Each Fold:
# The script saves training and test datasets for each fold of the cross-validation as separate R objects in the specified output directory.

# Example Usage:
# Users can run the script from the command line with various options to preprocess and split beta values and purities for cross-validation.

# Author: [Your Name]
# Affiliation: [Your Affiliation]


# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

#Loading optparse to parse command line arguments
if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))

#Loading catools to split the dataset in training and test
if(!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret") }

suppressPackageStartupMessages(library(caret))


# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

#Creating argument list
argument_list <- list(

  make_option(c("-s", "--in_same_object"), type="logical",  
              help="The user must especify if the purity and betas are in the same R object.",
              metavar = "[TRUE/FALSE]"),

  make_option(c("-B", "--path_to_directory"), type="character",  
              help="The path to the R object with the data to process must be entered here. If the purities and betas are in different R object files the one containing the betas should be entered here",
              metavar = "[path]"),

  make_option(c("-P", "--path_to_purities"), type="character",  
              help="The path to the R object with the purities must be entered here. This should only be used if the purities and betas are in different R objects.",
              metavar = "[path]"),

  make_option(c("-b", "--betas"), type="character",  
              help="The name of the object that contains the beta values inside the R object must be entered here",
              metavar = "[var_name]"),

  make_option(c("-u", "--purity"), type="character",  
              help="The name of the object that contains the beta values inside the R object must be entered here",
              metavar = "[var_name]")

  make_option(c("-k", "--groups_in_k_fold_cross_validation"), type="numeric",
              help="The number of groups in which the data will be splitted to perform the cross validation can be entered here",
              metavar = "[number]"),

  make_option(c("-o", "--output_directory"), type = "character", default =".",
              help = "The directory in which the output files will be stored can be. The default value is the current directory",
              metavar ="[path]"),

  make_option(c("-f", "--filter_non_autosomes"), type="logical", default=FALSE,
              help="Set this argument to TRUE to discard CpGs from sexual chromosomes", 
              metavar="[TRUE/FALSE]"),

  make_option(c("-A", "--path_to_annotation"), type="character",
              help="Specify the name of the path containing the annotation file.",
              metavar="[path]"),

  make_option(c("-a", "--annotation_dataframe_name"), type="character",
              help="Specify the name of the dataframe containing the annotation information.",
              metavar="[variable name]"),

  make_option(c("-c", "--column_chromosomes"), type="character",
              help="Select the name of column of the dataframe containig the information about the chromosomes.",
              metavar="[column name]"),

  make_option(c("-N", "--get_only_CpGs_with_highest_variance"), type="logical", default=FALSE,
              help="This argument should be sent to true if the user wants to get only the betas of the CpGs with the highest variance",
              metavar="[TRUE/FALSE]"),

  make_option(c("-n", "--number_of_most_variable_CpGs_to_include"), type="numeric",
              help="The number of the most variant CpGs to be included in the output should be included here",
              metavar="[number]")       

)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program gets purity and original beta Robjects from a R object containing those or more elements. It also allow the user to cerate training anjd test dataset"))


# ===========
# LODING DATA
# ===========

cat("\n\nLoading the data...\n")

#If the betas and purities are inside the same R object load a single file
if (arguments$in_same_object) { 

  load(arguments$path_to_directory) #Loding the R object with the purities and betas

#If the betas and purities are not inside the same R object load two files
} else {

  load(arguments$path_to_directory) #Loding the R object with the betas
  load(arguments$path_to_purities) #Loding the R object with the purities

}

# Getting betas and purities from the loaded R objects
betas <- get(arguments$betas)
purities <- get(arguments$purity)

#Adapting the names of the colnames of the beta value dataframe
colnames(betas) <- lapply(colnames(betas), 
                         function (name) {strsplit(name, "-01")[[1]][1]})

# ======================================
# FILTERING CpGs FROM SEXUAL CHROMOSOMES
# ======================================

# If the user has selected to filter sexual chromosomes the following code will be run
if (arguments$filter_non_autosomes) {
  
  cat("\n\nFiltering CpGs from sexual chromosomes...\n")

  #Loading annotation file with the chromosome info and getting the data
  load(arguments$path_to_annotation) 
  annotation <- get(arguments$annotation_dataframe_name)

  #Ordering the annotation file based on the order of the CpG beta dataframe
  annotation <- annotation[rownames(betas),]

  #Get a vector of CpGs which are not in the sexual chromosomes
  included_CpGs <- rownames(annotation)[!annotation[,arguments$column_chromosomes] %in% c("chrX", "chrY")]

  #Removing the CpGs belonging to the sexual chromosomes (the chromosomes that are not chrX a d chrY)
  betas <- betas[included_CpGs,]

}


# ====================================
# CREATING GROUPS FOR CROSS VALIDATION
# ====================================

    cat("\n\nSplitting the data for the cross validation...\n")

    # Using the createFoldas from the caret package to create groups for the cross validation
    # The number of groups was set to 6 generate have around 100 samples per group
    groups <- createFolds(y=colnames(betas),
                          k=6,
                          list=TRUE,
                          returnTrain=TRUE)


    # Generating a list with the training and test data for each created fold
    folds_list <- lapply(X=names(groups), 
                         FUN=function(fold) {
                            
                            # Getting the training and test beta values
                            training <- betas[,groups[[fold]]]
                            test <- betas[,-groups[[fold]]]

                            # Getting the training and test purity values
                            purity_training <- purities[colnames(training)]
                            purity_test <- purities[colnames(test)]

                            # Adding the created dataframes to a list and returning it
                            return(list("Training"=training, 
                                        "Test"=test, 
                                        "Purity_training"=purity_training,
                                        "Purity_test"=purity_test))               
                            }
                        )

    # Using the name of ecah fold as the id of the elements of fold_list
    names(folds_list) <- names(groups)



# =====================================================
# GETTING MOST VARIABLE CPGs BASED ON THE TRAINING DATA
# =====================================================

# If the user has selected to get only the most variable CpGs the following code will be run
if(arguments$get_only_CpGs_with_highest_variance) {

        for (fold in names(folds_list)) {

            # Detemine the variance of all the rows (CpGs)
            cpgs_variance <- apply(
                                folds_list[[fold]][["Training"]],
                                MARGIN=1,
                                FUN=var
                        )

            # Create a vector to sort the rows. Get only the number of rows containing the CpGs that want to be included
            sorting_vec <- order(cpgs_variance, decreasing=TRUE)
            sorting_vec <- sorting_vec[1:arguments$number_of_most_variable_CpGs_to_include]

            # Getting only the cpgs selected in the training and test datasets
            folds_list[[fold]][["Training"]] <- folds_list[[fold]][["Training"]][sorting_vec,]
            folds_list[[fold]][["Test"]] <- folds_list[[fold]][["Test"]][sorting_vec,]

        }
    }


# ===========================================================
# CREATING OUTPUT FILES FOR EACH FOLD OF THE CROSS VALIDATION
# ===========================================================

    cat("\n\nSaving output files...\n")

    for (fold in names(folds_list)) {

        # Saving training and test datasets
        saveRDS(folds_list[[fold]][["Training"]], paste(arguments$output_directory, "/", fold, "_BetasTraining.RData", sep=""))
        saveRDS(folds_list[[fold]][["Test"]], paste(arguments$output_directory, "/", fold, "_BetasTest.RData", sep=""))

        saveRDS(folds_list[[fold]][["Purity_training"]], paste(arguments$output_directory, "/", fold, "_PurityTraining.RData", sep=""))
        saveRDS(folds_list[[fold]][["Purity_test"]], paste(arguments$output_directory, "/", fold, "_PurityTest.RData", sep=""))
    }


cat("\n\n****************\n")
cat("PROCESS FINISHED\n")
cat("****************\n")