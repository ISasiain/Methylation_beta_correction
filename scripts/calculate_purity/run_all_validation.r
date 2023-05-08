#!/usr/bin/Rscript

# ==================================
# INSTALLING AND/OR LOADING PACKAGES
# ==================================

#Install and load optparse to parse command line arguments
if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", quietly = TRUE) }

suppressPackageStartupMessages(library(optparse))

#Install and load parallel to run the script in parallel
if(!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel", quietly = TRUE) }

suppressPackageStartupMessages(library(parallel))

#Install and load doSnow to run the script in parallel
if(!requireNamespace("doSNOW", quietly = TRUE)) {
  install.packages("doSNOW", quietly = TRUE) }

suppressPackageStartupMessages(library(doSNOW))

#Install and load foreach to run the script in parallel
if(!requireNamespace("foreach", quietly = TRUE)) {
  install.packages("foreach", quietly = TRUE) }

suppressPackageStartupMessages(library(foreach))

#Install and load Kendall to create a progress bar
if(!requireNamespace("Kendall", quietly = TRUE)) {
  install.packages("Kendall", quietly = TRUE) }

suppressPackageStartupMessages(library(Kendall))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

#Defining command line arguments
argument_list <- list(

  make_option(c("-c", "--cores"), type="integer", default=1,  
              help="Number of cores to be used to run the program [default %default]",
              metavar = "[number]"),
  
  make_option(c("-a", "--alpha"), type="double", default=0.7, 
              help="Alpha value to determine the prediction intervals of each CpG [default %default]",
              metavar= "[floating number]"),

  make_option(c("-r", "--max_RSE"), type="double", default=0.5,
              help="Maximum Residual Standard Error allowed per CpG regression [default %default]",
              metavar="[floating number]"),

  make_option(c("-s", "--min_slope"), type="double", default=0.3,
              help="Minimum slope allowed per CpG regression [default %default]", 
              metavar="[floating number]"),

  make_option(c("-p", "--percentage_to_interval"), type="double", default=4.0,
              help="Percentage of the maximum coverage to include in the 1-Purity interval [default %default]",
              metavar="[floating number]"),

  make_option(c("-d", "--regression_data"), type="character",
              help="The directory containing the regression parameters must be entered here.",
              metavar="[directory]"),

  make_option(c("-b", "--betas_to_analyse"), type="character",
              help="The path to the R object contaoining the betas to analyse must be entered here.",
              metavar="[directory]"),

  make_option(c("-o", "--output_filename"), type="character", default="prediected purities",
              help="The name of the output R object containing the predicted values must be entered here. Default [%default]",
              metavar="[filename]"),

  make_option(c("-l", "--output_location"), type="character", default="./",
              help="The name of the location output R object containing the predicted values must be entered here. Default [%default]",
              metavar="[path_to_directory]")  

)

#Parsing arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program estimates the Purity values of samples based on the beta values of each sample's CpGs."))


# ==============================
# LOADING THE REQUIRED FUNCTIONS
# ==============================

# Getting the directory of the functions sourced on this script 
# (They must be located on the same directory than the main script)
dir <- commandArgs()[4]
dir <- gsub("--file=", "", dir)

#Getting the complete path of the functions to be called
fun1 <- gsub("run_all_validation.r", "predicting_purity.r", dir)
fun2 <- gsub("run_all_validation.r", "purity_coverage.r", dir)

#Sourcing the functions
source(fun1)
source(fun2)

# ========================
# PREPARING THE INPUT DATA
# ========================


#Reading the R objects containing the regression data as dataframes
my_slopes <- readRDS(list.files(arguments$regression_data, pattern="*reg.slopes.RData", full.names=TRUE))
my_intercepts <- readRDS(list.files(arguments$regression_data, pattern="*reg.intercepts.RData", full.names=TRUE))
my_RSE <- readRDS(list.files(arguments$regression_data, pattern="*reg.RSE.RData", full.names=TRUE))
my_df <- readRDS(list.files(arguments$regression_data, pattern="*reg.df.RData", full.names=TRUE))


#Reading beta values
unadj_validation <- readRDS(arguments$betas_to_analyse)

#Create a list to append all the predicted purity intervals
list_of_predicted_intervals <- list()

# ===========================
# CONFIGURING PARALLELIZATION
# ===========================

# Printing the number of cores to be used
cat("\nUsing", arguments$cores, "cores\n")

# Creatinhg clusters to run the script in  parallel
cl <- makeCluster(arguments$cores)

#Registering the clusters
registerDoSNOW(cl)

# ========================================
# RUNNING THE ANALYSIS WITH A PROGRESS BAR
# ========================================

# Printing command line message
cat("\nRunning the analysis...\n\n")

# Getting the names of the samples to analyse
samples <- colnames(unadj_validation)

# Defining the progress bar
p_bar <- txtProgressBar(min = 0, 
                        max = length(samples), 
                        style = 3)

# Creating a function to follow the execution of the script
progress <- function(n) setTxtProgressBar(p_bar, n)
opts <- list(progress = progress)

# Running the sourced functions in parallel for each sample. The execution level will be followed through a progress bar
list_of_predicted_intervals <- foreach(s = samples, .packages = "Kendall", .options.snow = opts) %dopar% {

  # Defining an empty dataframe to add the all the 1-Purity predicted intervals for all the CpGs of a sample
  interval_df <- data.frame()

  # Predicting all the 1-Purity intervals for each CpG of each sample and append them to the empty interval_df
  for (cpg in rownames(unadj_validation)) {
    interval_df <- rbind(interval_df, predicting_purity(beta=unadj_validation[cpg, s],
                                                        slopes=my_slopes[cpg, ],
                                                        intercepts=my_intercepts[cpg, ],
                                                        RSE=my_RSE[cpg, ],
                                                        degrees_of_freedom=my_df[cpg, ],
                                                        slope_threshold=arguments$min_slope,
                                                        RSE_threshold=arguments$max_RSE,
                                                        alpha=arguments$alpha
                                                        ))
  }

  # Calculate the 1-Purity estimate and interval for the sample analysed.
  # The results with be shown in list named with the sample id
  list(name = s, value = purity_value_per_sample(pred_purity_confidence=interval_df,
                                                 interval_threshold=arguments$percentage_to_interval))
}

# Stop the defined clusters
stopCluster(cl)

# Append the list defined for each sample to the list containing the predicted values for all the samples.
# The sample id is used to identify each element of the list
list_of_predicted_intervals <- setNames(lapply(list_of_predicted_intervals, function(x) x$value), sapply(list_of_predicted_intervals, function(x) x$name))

# Save the list of predicted values as an R object
saveRDS(list_of_predicted_intervals, file=paste(arguments$output_location, arguments$output_filename, ".RData", sep=""))

#Print message to show the end of the execution
cat("\n\n**********************\n")
cat("   PROCESS FINISHED\n")
cat("**********************\n")
