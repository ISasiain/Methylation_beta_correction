#!/usr/bin/Rscript

cat("\nPreparing the working envirnonemnt...\n")

# ==================================
# INSTALLING AND/OR LOADING PACKAGES
# ==================================

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

num_of_clusters <- detectCores() - 1

cat("\nUsing", num_of_clusters, "cores\n")

cl <- makeCluster(num_of_clusters)
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
    interval_df <- rbind(interval_df, predicting_purity(unadj_validation[cpg, s], my_slopes[cpg, ], my_intercepts[cpg, ], my_RSE[cpg, ], my_df[cpg, ]))
  }

  # Return a named list directly
  list(name = s, value = purity_value_per_sample(interval_df))
}

stopCluster(cl)

list_of_predicted_intervals <- setNames(lapply(list_of_predicted_intervals, function(x) x$value), sapply(list_of_predicted_intervals, function(x) x$name))

saveRDS(list_of_predicted_intervals, file="out_with_reg_test2.RData")

cat("\n\n**********************\n")
cat("   PROCESS FINISHED\n")
cat("**********************\n")
