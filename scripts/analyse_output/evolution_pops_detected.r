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

if(!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr") }

library(stringr)

# =============================================================
#      PARSING COMMAND LINE ARGUMENTS AND LOADING THE DATA
# =============================================================

argument_list <- list(

    make_option(c("-f", "--list_of_files"), type="character",
                help="The user must especify the path of the R objects containing the intercepts of the identified populations to analyse separated by commas.",
                metavar="[comma_separated_paths]")

)

# =============================================
#      GETTING THE INTERCEPT DATA TO COMPARE
# =============================================

# Splitting the provided argument in a vector of paths
vec_of_paths <- strsplit(arguments$predictions_to_compare, split=",")[[1]]

# Creating a list to store the predictions to compare
pop_num_ls <- list()

# Appending each prediction to the prediction list
for (path in vec_of_paths) {

    # Getting the variable name to store the predicted purities
    var_name <- strsplit(path, split="/")[[1]][length(strsplit(path, split="/")[[1]])]
    var_name <- str_remove(var_name, "output_reg.intercepts.RData")

    print(var_name)

    # Append elements to the list using the value stored in var_name as id
    prediction_ls[[var_name]] <- rowSums(!is.na(readRDS(path)))
}


# ==========================================================
#      TRANSFORMING DATA INTO A DATAFRAME TO BE PLOTTED
# ==========================================================

# Creating a dataframe to append the proportions of cpgs with one, two or three populations
pops_proportions <- data.frame()

#Adding one column per prtediction analysed
for (prediction in names(prediction_ls)) {

    pops_proportions <- cbind(pops_proportions, c(length(prediction_ls[predction][prediction_ls[predction]==3])/length(prediction_ls[predction]),
                                                  length(prediction_ls[predction][prediction_ls[predction]==2])/length(prediction_ls[predction]),
                                                  length(prediction_ls[predction][prediction_ls[predction]==1])/length(prediction_ls[predction])))

}

pops_proportions <- cbind(pops_proportions,names(prediction_ls))

# Adding rownames and colnames
colnames(pops_proportions) <- c("prop_3", "prop_2", "prop_1", "prediction")
rownames(pops_proportions) <- names(prediction_ls)


# ==================================
#       PLOTTING THE DATAFRAME
# ==================================

ggplot(data=pops_proportions) +
    geom_bar(x=predction)
