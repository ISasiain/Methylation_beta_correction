#!/usr/bin/Rscript

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

    make_option(c("-d", "--directory"), type="character",
    help="The name of the directory containing the different predictions made and different folds for each must be entered here",
    metavar="[path_to_dorectory]"),

    make_option(c("-p", "--reference_purity"), type="character", 
    help="The path of the R object containing the refernce purity values to compare the predictions with must be entered here",
    metavar="[path_to_purities]"),

    make_option(c("-o", "--output_prefix"), type="character", default="result_analysis",
    help="The user can especify a prefix for the output plots. Default [%default]",
    metavar = "[prefix]")
)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This script compares the output of different sample purity prtedictions made based on beta values"))


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

# Create an empty list to store predictions
prediction_list <- vector("list", length(vec_of_preds))

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
      stop("More than one file containing the estimated purity per prediction and fold has been detected. Execution halted.")
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

print(vec_of_preds)
print(vec_of_folds)

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

# Plotting the results 
ggplot(dti_df, aes(x = factor(Prediction, levels = vec_of_preds), y = Distance_to_estimate)) +
  geom_boxplot(fill="olivedrab4") +  # Add a boxplot for each prediction
  labs(x = "Prediction", y = "Mean distance to estimate (% points)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Saving the plot
ggsave(paste(arguments$output_prefix, "cross_validation.boxplot.png", sep="."))