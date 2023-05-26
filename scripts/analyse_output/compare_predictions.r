#!/usr/bin/Rscript

# ==========================================
#      IMPORTING THE REQUIRED PACKAGES 
# ==========================================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

if(!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2") }

library(ggplot2)

if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

library(optparse)


# ==========================================
#      PARSING COMMAND LINE ARGUMENTS
# ==========================================

argument_list <- list(

    make_option(c("-c", "--predictions_to_compare"), type="character", 
    help="The paths to the R objects cotaining different predictions must be entered here separated by a coma",
    metavar="[comma_separated_list_of_paths]"),

    make_option(c("-p", "--reference_purity"), type="character", 
    help="The path of the R object containing the refernce purity values to compare the predictions with must be entered here",
    metavar="[path_to_purities]"),

  make_option(c("-o", "--output_prefix"), type="character", default="result_analysis",
              help="The user can especify a prefox for the output plots. Default [%default]",
              metavar = "[prefix]")
)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This script compares the output of different sample purity prtedictions made based on beta values"))

# =============================================
#      GETTING THE PREDICTIONS TO COMPARE
# =============================================

# Splitting the provided argument in a vector of paths
vec_of_paths <- strsplit(arguments$predictions_to_compare, split=",")[[1]]

# Creating a list to store the predictions to compare
prediction_ls <- list()

# Appending each prediction to the prediction list
for (path in vec_of_paths) {

    # Getting the variable name to store the predicted purities
    var_name <- strsplit(path, split="/")[[1]][length(strsplit(path, split="/")[[1]])]
    var_name <- strsplit(var_name, split="\\.")[[1]][1]


    # Append elements to the list using the value stored in var_name as id
    prediction_ls[[var_name]] <- readRDS(path)
}

# ======================================
#      GETTING THE ACTUAL PURITIES
# ======================================

# Loading the actual 1-Purity vector into a variable
actual_1_minus_P <- readRDS(arguments$reference_purity)

#Adding this to adapt the sample names
names(actual_1_minus_P) <- paste(names(actual_1_minus_P), "-01A", sep="")

# =======================================
#      CALCULATING QUALITY METRICS
# =======================================

# Create an empty dataframe to store the calculated data
metrics_df <- data.frame()

#Calculating metrics of each prediction and appending them to the dataframe as new lines
for (prediction in names(prediction_ls)) {

  mean_dis_to_est <- mean(
    sapply(
      names(prediction_ls[[prediction]]),
      FUN = function (sample) {
      # Determining the distance to the estimate (the mean is calculated in case there were more than one maximums. THis must be removed after softening the coverage plot)
      dis_to_estimate <- mean(abs(1-actual_1_minus_P[sample] - prediction_ls[[prediction]][[sample]][["1-Pur_estimates"]]),4)
      }
    )
  )

  # Calculating correlation between actual and estimate purities. The estimated values are obtained using an lapply command
  # and sorted as the vector containing the actual 1-Purity values

  correlation <- cor(actual_1_minus_P, sapply(names(prediction_ls[[prediction]]), FUN = function(sample) {1-prediction_ls[[prediction]][[sample]][["1-Pur_estimates"]]})[names(actual_1_minus_P)])

  metrics_df <- rbind(metrics_df, c(prediction, mean_dis_to_est, correlation))

}

# Renaming columns and row of the dataframe
colnames(metrics_df) <- c("Prediction", "Mean_dis_to_est", "Correlation")
rownames(metrics_df) <- metrics_df$Prediction

#Transform metrics to numeric values
metrics_df$Correlation <- as.numeric(metrics_df$Correlation)
metrics_df$Mean_dis_to_est <- as.numeric(metrics_df$Mean_dis_to_est)

# ====================================================================
#      SORTING DATAFRAME (only if numeric parameters are compared)
# ====================================================================

# Apply strsplit to each element of the Prediction column
to_order <- as.numeric(sapply(metrics_df$Prediction, function(x) strsplit(x, split = "_")[[1]][2]))

# Sorting the dataframe
metrics_df <- metrics_df[order(to_order),]

# ===========================
#      PLOTTING RESULTS
# ===========================

#Setting the prediction id as a factor
metrics_df$Prediction <- factor(metrics_df$Prediction, levels=metrics_df$Prediction)

print(metrics_df)

# Plotting and saving a barplot of the distance to estimate

ggplot(metrics_df, aes(x = Prediction, y = Mean_dis_to_est)) +
  geom_bar(stat = "identity", fill = "blue") +
  xlab("Prediction") +
  ylab("Mean Distance to Estimate") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(paste(arguments$output_prefix, "dis_to_int.barplot_comp.png", sep="."))

# Plotting and saving a barplot of the correlation

ggplot(metrics_df, aes(x = Prediction, y = Correlation)) +
  geom_bar(stat = "identity", fill = "blue") +
  xlab("Prediction") +
  ylab("Correlation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(paste(arguments$output_prefix, "correlation.barplot_comp.png", sep="."))