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

#Adding this to adapt the sample names
names(actual_1_minus_P) <- paste(names(actual_1_minus_P), "-01A", sep="")


# ESTIMATED PURITIES

# Get all the combinations of predictions and folds
combinations <- expand.grid(vec_of_preds, vec_of_folds)

# Creating a function to append each prediction to the list

append_prediction <- function(prediction, fold, directory) {

    #Getting the path to each prediction (SLIGHTLY HARDCODED!!!  WE ARE ASSUMING THAT THE FILE ENDS WITH cpgXXX.RData)
    path_to_prediction <- paste(directory, "/", prediction, "/", fold, sep="")
    pattern <- paste("*_", gsub("_", "", prediction), "\\.RData", sep="")

    # List files in the directory that match the pattern
    matching_files <- list.files(path_to_prediction, pattern = pattern, full.names = TRUE)

    #PROBLEM!!! SOME FILES HAVE NOT BEEN GENERATED!!!!!!!1
    print(prediction)
    # Print the matching files
    print(matching_files)

    if (length(matching_files) != 1) {

      stop("More than one file containing the estimated purity per prediction and fold has been detected. Execution halted.")

    }

    #Appending each prediction
    prediction_list[[prediction]][[fold]] <- readRDS(matching_files[1])
}

#Applying the function to each row of the combinations dataframe using lapply

lapply(X=seq_len(nrow(combinations)), 
       FUN=function(i) {append_prediction(prediction=combinations$Var1[i],
                                          fold=combinations$Var2[i], 
                                          directory=arguments$directory)
                        }
)
                
# =====================================================================
#          DETERMINING QUALITY METRICS (distance to estimate)
# =====================================================================

# Creating an empty dataframe to store the results
dti_df <- data.frame(matrix(nrow = length(vec_of_preds), ncol = length(vec_of_folds)))

# Adding column names and row names
colnames(dti_df) <- vec_of_folds
rownames(dti_df) <- vec_of_preds

for (prediction in vec_of_preds) {

  for (fold in vec_of_folds) {

    dti_df[prediction, fold] <- mean(
      sapply(
        names(prediction_ls[[prediction]][[fold]]),
        FUN = function (sample) {
        # Determining the distance to the estimate (the mean is calculated in case there were more than one maximums. THis can be removed after softening the coverage plot)
        dis_to_estimate <- mean(abs(1-actual_1_minus_P[sample] - prediction_ls[[prediction]][[fold]][[sample]][["1-Pur_estimates"]]),4)
        }
      )
    )
  }
}


# ================================================
#                PLOTTING RESULTS
# ================================================


# Transforming the row names into a column (using tibble)
dti_df <- dti_df %>% rownames_to_column("Prediction")


# Transforming the dataframe to long format to be plotted
dti_df_long <- pivot_longer(dti_df, cols="Prediction", names_to="Fold", values_to="Distance_to_estimate")

# Plotting the results 

ggplot(dti_df_long, aes(x = Prediction, y = Quality, color = Fold)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
  labs(x = "Prediction", y = "Quality Metric") +
  theme_classic()

# Saving the plot
ggsave(paste(arguments$output_prefix, "cross_validation.boxplot.png", sep="."))
