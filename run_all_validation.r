#Installing packages if necessary
if(!requireNamespace("caTools", quietly = TRUE)) {
  install.packages("caTools") }

library(caTools)

if(!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel") }

library(doParallel)

if(!requireNamespace("foreach", quietly = TRUE)) {
  install.packages("foreach") }

library(foreach)

#Setting the working directory
setwd("~/Methylation/adjustBetas")

#Sourcing the functions to be used
source("./predicting_purity.R")
source("./purity_coverage.R")

##data set reduced to top 5000 CpGs in set of 630 TCGA BRCA tumors
load("workspace_tcgaBrca_top5000.RData")


#Reading csv files as dataframes
total_betas <- read.csv("output_betas.original.csv", sep=",", row.names=1)
my_slopes <- read.csv("output_training_reg.slopes.csv", sep=",", row.names=1)
my_intercepts <- read.csv("output_training_reg.intercepts.csv", sep=",", row.names=1)
my_RSE <- read.csv("output_training_reg.RSE.csv", sep=",", row.names=1)
my_SSx <- read.csv("output_training_reg.SSx.csv", sep=",", row.names=1)
my_SSy <- read.csv("output_training_reg.SSy.csv", sep=",", row.names=1)
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

for (sample in colnames(unadj_validation)) {

    cat("\nSamples analysed;",which(colnames(unadj_validation)==sample), "/", ncol(unadj_validation), "\n")
    interval_df <- data.frame()

    for (cpg in rownames(unadj_validation)) {

        interval_df <- rbind(interval_df, predicting_purity(unadj_validation[cpg,sample], my_slopes[cpg,], my_intercepts[cpg,], my_RSE[cpg,], my_SSx[cpg,], my_SSy[cpg,], my_df[cpg,]))
    }

    list_of_predicted_intervals[[sample]] <- purity_value_per_sample(interval_df)

}

saveRDS(list_of_predicted_intervals, file="out_without_reg.RData")

cat("\n**********************\n")
cat("   PROCESS FINISHED\n")
cat("**********************\n")
