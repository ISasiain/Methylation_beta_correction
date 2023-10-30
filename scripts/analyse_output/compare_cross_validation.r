#!/usr/bin/Rscript

# Script Name: compare_cross_validation.r

# Description:
# This script compares and analyzes sample purity predictions from different datasets. It takes as input R objects containing predicted purity values for various folds and cancer types and generates boxplots and quality metrics to assess the predictions.

# Required R Packages:
# - ggplot2: For creating data visualizations.
# - stringr: For string manipulation.
# - tibble: For data manipulation.
# - tidyr: For data tidying.
# - dplyr: For data manipulation.
# - optparse: For parsing command-line arguments.

# Loading Required Packages:
# - The script sets a specific CRAN mirror and installs the required packages if not already installed.
# - It also loads the 'ggplot2', 'stringr', 'tibble', 'tidyr', 'dplyr', and 'optparse' packages for data analysis and argument parsing.

# Configuration of Command Line Arguments:
# - The script uses the 'optparse' package to configure and parse command-line arguments.
# - The available options include:
#   - "--names_of_the_folds": Names of the different folds used, separated by a comma.
#   - "--names_of_the_predictions_to_compare": Names of the different predictions to compare, separated by a comma.
#   - "--plot_multiple_cancer_types": Set to TRUE if multiple cancer types will be plotted.
#   - "--directory/es": The name of the directory, or directories separated by a comma, containing the different predictions and folds.
#   - "--reference_purity/es": The path of the R object, or R objects separated by a comma, containing the reference purity values.
#   - "--output_prefix": Prefix for the output plots.

# Parsing Command Line Arguments:
# - The 'optparse' package is used to parse the command-line arguments, including names of folds, predictions, directories, reference purities, and output prefix.

# Processing Command Line Arguments:
# - The script processes the provided arguments, splitting them into vectors and ordering them if necessary.

# Getting Values to Compare:
# - The script loads actual purities and predicted values for each prediction, fold, and cancer type.
# - It calculates a quality metric, the distance to the estimate, for each prediction.

# Determining Quality Metrics:
# - The script determines quality metrics for each prediction and fold and stores the results in a data frame.

# Plotting Results:
# - The script creates boxplots and quality metric visualizations for the predictions.
# - The generated plots are saved as PNG files with appropriate names based on the dataset type.

# Example Usage:
# - Users can run the script from the command line with various options to compare and analyze sample purity predictions.

# Author: Iñaki Sasiain Casado
# Affiliation: Johan Staaf lab @ Lund University / Oncology & Pathology



# ==========================================
#      IMPORTING THE REQUIRED PACKAGES 
# ==========================================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

if(!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2") }

library(ggplot2)

if(!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr") }

library(stringr)

if(!requireNamespace("tibble", quietly = TRUE)) {
  install.packages("tibble") }

library(tibble)

if(!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr") }

library(tidyr)

if(!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr") }

library(dplyr)

if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

library(optparse)


# ==========================================
#      PARSING COMMAND LINE ARGUMENTS
# ==========================================

argument_list <- list(

    make_option(c("-f", "--names_of_the_folds"), type="character", 
    help="The names that identify the different folds used must be entered here separated by a comma",
    metavar="[comma_separated_list_of_folds]"),

    make_option(c("-c", "--names_of_the_predictions_to_compare"), type="character", 
    help="The names that identify the different predictions must be entered here separated by a comma",
    metavar="[comma_separated_list_of_prediction_names]"),

    make_option(c("-m","--plot_multiple_cancer_types"), type="logical", default=FALSE,
    help="The user can specify if more than one cancer types have to be included in the plot",
    metavar="[TRUE/FALSE]"),

    make_option(c("-d", "--directory/es"), type="character",
    help="The name of the directory, or directories separacted by a comma if the user has selected to plot different cancer types, containing the different predictions made and different folds for each must be entered here",
    metavar="[path_to_dorectory]"),

    make_option(c("-p", "--reference_purity/es"), type="character", 
    help="The path of the R object, or R objects separated by a comma is the user has selected to plot different cancer types, containing the refernce purity values to compare the predictions with must be entered here",
    metavar="[path_to_purities]"),

    make_option(c("-o", "--output_prefix"), type="character", default="result_analysis",
    help="The user can especify a prefix for the output plots. Default [%default]",
    metavar = "[prefix]")
)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This script compares the output of different sample purity prtedictions made based on beta values"))

#Code to use if a single cancer type will be plotted
if (arguments$plot_multiple_cancer_types==FALSE) {

# ===================================================
#         PROCESSING COMMAND LINE ARGUMENTS
# ===================================================

# Splitting the provided argument in a vector of folds
vec_of_folds <- strsplit(arguments$names_of_the_folds, split=",")[[1]]

# Splitting the provided argument in a vector of predictions to analyse
vec_of_preds <- strsplit(arguments$names_of_the_predictions_to_compare, split=",")[[1]]

# Apply strsplit to each element of the Prediction column
to_order <- as.numeric(sapply(vec_of_preds, function(x) strsplit(x, split = "_")[[1]][2]))

# Sorting the vector
vec_of_preds <- vec_of_preds[order(to_order)]


# Creating an empty list for the predictied values per prediction and fold

# Creating a list containing each prediction made as an element
prediction_list <- vector("list", length(vec_of_preds))
names(prediction_list) <- vec_of_preds

# Assign a list with the fold name as id to each prediction made
for (prediction in names(prediction_list)) {

    prediction_list[[prediction]] <- vector("list", length(vec_of_folds))
    names(prediction_list[[prediction]]) <- vec_of_folds

}

# ===============================================
#            GETTING VALUES TO COMPARE
# ===============================================


# ACTUAL PURITIES

# Loading the actual 1-Purity vector into a variable
actual_1_minus_P <- readRDS(arguments$reference_purity)

# Iterate through predictions and folds and append the predictions
for (pred in vec_of_preds) {
  for (fold in vec_of_folds) {

    # Getting the path to each prediction
    path_to_prediction <- paste(arguments$directory, "/", pred, "/", fold, sep="")
    pattern <- paste("*_", gsub("_", "", pred), "\\.RData", sep="")


    # List files in the directory that match the pattern
    matching_files <- list.files(path_to_prediction, pattern = pattern, full.names = TRUE)


    # Raise an error if more than one paths are identified
    if (length(matching_files) != 1) {

      stop("More than one or 0 files containing the estimated purity per prediction and fold has been detected. Execution halted.")
    }

    #Append to list
    prediction_list[[pred]][[fold]] <- readRDS(matching_files[1])

  }
}


# =====================================================================
#          DETERMINING QUALITY METRICS (distance to estimate)
# =====================================================================

# Creating an empty dataframe to store the results
dti_df <- data.frame(matrix(nrow = length(vec_of_preds), ncol = length(vec_of_folds)))

# Adding column names and row names
rownames(dti_df) <- vec_of_preds
colnames(dti_df) <- vec_of_folds

# Iterate through the predictions and folds to determine the dist_to_est (quality metric) per each
for (pred in vec_of_preds) {

  for (fold in vec_of_folds) {

    dti_df[pred, fold] <- mean(
      sapply(
        names(prediction_list[[pred]][[fold]]),
        FUN = function (sample) {
        

        # Determining the distance to the estimate (the mean is calculated in case there were more than one maximums. THis can be removed after softening the coverage plot)
        mean(abs(1 - actual_1_minus_P[sample] - prediction_list[[pred]][[fold]][[sample]][["1-Pur_estimates"]]),4)
        }
      )
    )
  }
}

print(dti_df)

# ================================================
#                PLOTTING RESULTS
# ================================================


# Transforming the row names into a column (using tibble)
dti_df <- dti_df %>% rownames_to_column("Prediction")

# Transforming the dataframe to long format to be plotted
dti_df <- dti_df %>% pivot_longer(cols=all_of(vec_of_folds), 
                                  names_to="Fold", 
                                  values_to="Distance_to_estimate")

print(dti_df)

# Plotting the results 
ggplot(dti_df, aes(x = factor(Prediction, levels = vec_of_preds), y = Distance_to_estimate)) +
  geom_boxplot(fill="palegreen1", ) +
  geom_smooth(aes(group=1), linewidth=2, se=FALSE, color="palegreen4") +
  labs(x = "Prediction", y = "Mean distance to estimate") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Saving the plot
ggsave(paste(arguments$output_prefix, "cross_validation.boxplot.png", sep="."))


#Code to use if multiple cancer types will be plotted
} else {


# ===================================================
#         PROCESSING COMMAND LINE ARGUMENTS
# ===================================================

# Splitting the provided argument in a vector of folds
vec_of_folds <- strsplit(arguments$names_of_the_folds, split=",")[[1]]

# Splitting the provided argument in a vector of predictions to analyse
vec_of_preds <- strsplit(arguments$names_of_the_predictions_to_compare, split=",")[[1]]

# Apply strsplit to each element of the Prediction column
to_order <- as.numeric(sapply(vec_of_preds, function(x) strsplit(x, split = "_")[[1]][2]))

# Sorting the vector
vec_of_preds <- vec_of_preds[order(to_order)]


# Defining vector directories with predictions
directories_preds <- strsplit(arguments$"directory/es", ",")[[1]]

#Defining vector with path to actual purities
pur_paths <- strsplit(arguments$"reference_purity/es", ",")[[1]]


#Creating vector with cancer type names
types <- lapply(directories_preds, function(x) tail(strsplit(x,"/")[[1]], n=1))

# Creating a list containing each prediction made as an element
prediction_list <- vector("list", length(types))
names(prediction_list) <- types

# Creating an empty list for type, prediction and fold
for (type in types) {

  prediction_list[[type]] <- vector("list", length(vec_of_preds))
  names(prediction_list[[type]]) <- vec_of_preds

  for (prediction in vec_of_preds)

    prediction_list[[type]][[prediction]] <- vector("list", length(vec_of_folds))
    names(prediction_list[[type]][[prediction]]) <- vec_of_folds
  }


# Creating an empty to store actual purities
actual_p_list <- vector("list", length(types))
names(actual_p_list) <- types


# ===============================================
#            GETTING VALUES TO COMPARE
# ===============================================


# FILLING LISTS; actual and predicted

for (type in types) {

  dir <- grep(type, directories_preds, value = TRUE)
  path_actual_p <- grep(type, pur_paths, value = TRUE)


  # Getting the actual purities for the current cancer type
  actual_p_list[[type]] <- readRDS(path_actual_p)

  for (pred in vec_of_preds) {
    for (fold in vec_of_folds) {

        # Getting the path to each prediction
    path_to_prediction <- paste(dir, "/", pred, "/", fold, sep="")
    pattern <- paste("*_", gsub("_", "", pred), ".RData", sep="")

    # List files in the directory that match the pattern
    matching_files <- list.files(path_to_prediction, pattern = pattern, full.names = TRUE)

    print(matching_files)


    # Raise an error if more than one paths are identified
    if (length(matching_files) != 1) {

      stop("More than one or 0 files containing the estimated purity per prediction and fold has been detected. Execution halted.")
    }

    #Append to list
    prediction_list[[type]][[pred]][[fold]] <- readRDS(matching_files[1])

    }
  }
}


# ===============================================
#            GETTING QUALITY METRICS
# ===============================================

# Results dataframe
to_plot_df <- data.frame(matrix(ncol=4, nrow=0))


#Iterate through cancer types
for (type in types){

  # Iterate through the predictions and folds to determine the dist_to_est (quality metric) per each
  for (pred in vec_of_preds) {

    for (fold in vec_of_folds) {

      to_plot_df <- rbind(to_plot_df, c(type, pred, fold, 
        mean(sapply(
          names(prediction_list[[type]][[pred]][[fold]]),
          FUN = function (sample) {
          

          # Determining the distance to the estimate (the mean is calculated in case there were more than one maximums. THis can be removed after softening the coverage plot)
          mean(abs(1 - actual_p_list[[type]][sample] - prediction_list[[type]][[pred]][[fold]][[sample]][["1-Pur_estimates"]]),4)
            }
          )
        ) 
      )
      )
    }
  }
}

#Adding colnames to the new plot
colnames(to_plot_df) <- c("Type", "Prediction", "Fold", "Distance_to_estimate")

print(to_plot_df)

# ===============================================
#              PLOTTING THE RESULTS
# ===============================================

# Create a boxplot
ggplot(to_plot_df, aes(x = factor(Prediction, levels = vec_of_preds), y = as.numeric(Distance_to_estimate) , fill = Type)) +
  geom_boxplot() +
  geom_smooth(aes(group=Type, color=Type), linewidth=2, se=FALSE) +
  labs(x = "Prediction", y = "Mean distance to estimate", fill = "Cancer Type", color="Mean dis. to est.") +
  scale_fill_manual(values = c("palegreen1", "tomato", "steelblue1")) +
  scale_color_manual(values = c("palegreen4", "tomato4", "steelblue4")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Saving the plot
ggsave(paste(arguments$output_prefix, ".cancer_types_combined.cross_validation.boxplot.png", sep="."))

}

cat("\n=================\n")
cat ("PROCESS FINISHED")
cat("\n=================\n")