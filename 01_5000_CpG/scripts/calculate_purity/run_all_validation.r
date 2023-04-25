#!/usr/bin/Rscript

# ==================================
# INSTALLING AND/OR LOADING PACKAGES
# ==================================

if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", quietly = TRUE) }

suppressPackageStartupMessages(library(optparse))

if(!requireNamespace("caTools", quietly = TRUE)) {
  install.packages("caTools", quietly = TRUE) }

suppressPackageStartupMessages(library(caTools))

if(!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel", quietly = TRUE) }

suppressPackageStartupMessages(library(parallel))

if(!requireNamespace("doSNOW", quietly = TRUE)) {
  install.packages("doSNOW", quietly = TRUE) }

suppressPackageStartupMessages(library(doSNOW))

if(!requireNamespace("foreach", quietly = TRUE)) {
  install.packages("foreach", quietly = TRUE) }

suppressPackageStartupMessages(library(foreach))

if(!requireNamespace("Kendall", quietly = TRUE)) {
  install.packages("Kendall", quietly = TRUE) }

suppressPackageStartupMessages(library(Kendall))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

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
              metavar="[floating number]")
)

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program estimates the Purity values of samples based on the beta values of each sample's CpGs."))

# =====================================================================
# SETTING WORKING DIRECTORY AND LOADING THE REQUIRED DATA AND FUNCTIONS
# =====================================================================

#Setting the working directory
setwd("~/Methylation/adjustBetas")

#Sourcing the functions to be used
source("./predicting_purity.R")
source("./purity_coverage.R")

#Loading data to workspace
load("workspace_tcgaBrca_top5000.RData")

# ========================
# PREPARING THE INPUT DATA
# ========================

#Reading csv files as dataframes
total_betas <- read.csv("output_training_betas.original.csv", sep=",", row.names=1)
my_slopes <- read.csv("output_training_reg.slopes.csv", sep=",", row.names=1)
my_intercepts <- read.csv("output_training_reg.intercepts.csv", sep=",", row.names=1)
my_RSE <- read.csv("output_training_reg.RSE.csv", sep=",", row.names=1)
my_df <- read.csv("output_training_reg.df.csv", sep=",", row.names=1)

#Get the validation data
#Setting the seed for spliting the datainterval_
set.seed(1)

#Splitting the data
train_samples <- sample.split(colnames(betaUnadjusted), SplitRatio=0.8)

#Creating validation datasets for betas and purity
unadj_validation <- betaUnadjusted[,!train_samples]
purity_validation <- purityVector[!train_samples]

#Create a list to append all the predicted purity intervals
list_of_predicted_intervals <- list()

# ===========================
# CONFIGURING PARALLELIZATION
# ===========================

cat("\nUsing", arguments$cores, "cores\n")

cl <- makeCluster(arguments$cores)
registerDoSNOW(cl)

# ========================================
# RUNNING THE ANALYSIS WITH A PROGRESS BAR
# ========================================

cat("\nRunning the analysis...\n\n")

samples <- colnames(unadj_validation)

p_bar <- txtProgressBar(min = 0, 
                        max = length(samples), 
                        style = 3)

progress <- function(n) setTxtProgressBar(p_bar, n)

opts <- list(progress = progress)

list_of_predicted_intervals <- foreach(s = samples, .packages = "Kendall", .options.snow = opts) %dopar% {

  interval_df <- data.frame()
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

  # Return a named list directly
  list(name = s, value = purity_value_per_sample(pred_purity_confidence=interval_df,
                                                 interval_threshold=arguments$percentage_to_interval))
}

stopCluster(cl)

list_of_predicted_intervals <- setNames(lapply(list_of_predicted_intervals, function(x) x$value), sapply(list_of_predicted_intervals, function(x) x$name))

saveRDS(list_of_predicted_intervals, file="out_with_reg_test2.RData")

cat("\n\n**********************\n")
cat("   PROCESS FINISHED\n")
cat("**********************\n")
