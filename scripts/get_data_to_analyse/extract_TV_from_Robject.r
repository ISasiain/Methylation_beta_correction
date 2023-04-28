#!/usr/bin/Rscript

# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")


if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))

if(!requireNamespace("caTools", quietly = TRUE)) {
  install.packages("caTools") }

suppressPackageStartupMessages(library(caTools))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

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

  make_option(c("-p", "--purity"), type="character",  
              help="The name of the object that contains the beta values inside the R object must be entered here",
              metavar = "[var_name]"),

  make_option(c("-S", "--split_in_training_and_validation"), type="logical", default=FALSE,
              help="This argument must be entered to specify if the data has to be splitted to create training and validation datasets.",
              metavar = "[boolean]"),

  make_option(c("-v", "--percentage_to_validation"), type="double", default=20.0,
              help="The percentage of the samples to be included in the validation data must be entered here.",
              metavar = "[number]"),

  make_option(c("-o", "--output_directory"), type = "character", default =".",
              help = "The directory in which the output files will be stored can be. The default value is the current directory",
              metavar ="[path]"),

  make_option(c("-f", "--filter_non_autosomes"), type="logical", default=FALSE,
              help="Set this argument to TRUE to discard CpGs from sexual chromosomes", 
              metavar="[TRUE/FALSE]"),

  make_option(c("-A", "--path_to_annotation"), type="character",
              help="Select the name of the path containing the annotation file.",
              metavar="[path]"),

  make_option(c("-a", "--annotation_dataframe_name"), type="character",
              help="Select the name of the dataframe containing the annotation information.",
              metavar="[variable name]"),

  make_option(c("-c", "--column_chromosomes"), type="character",
              help="Select the name of column of the dataframe containig the information about the chromosomes.",
              metavar="[column name]")            

)

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This programm gets purity and original beta Robjects from a R object containing those or more elements. It also allow the user to cerate training anjd validation dataset"))


# ================
# LODING R OBJECT 
# ================

cat("\n\nLoading the data...\n")

if (arguments$in_same_object) {

  load(arguments$path_to_directory) #Loding the R object with the purities and betas

} else {

  load(arguments$path_to_directory) #Loding the R object with the betas
  load(arguments$path_to_purities) #Loding the R object with the purities

}

betas <- get(arguments$betas)
purities <- get(arguments$purity)


# ======================================
# FILTERING CpGs FROM SEXUAL CHROMOSOMES
# ======================================

if (arguments$filter_non_autosomes) {
  
  cat("\n\nFiltering CpGs from sexual chromosomes...\n")

  #Loading the file
  load(arguments$path_to_annotation) 
  annotation <- get(arguments$annotation_dataframe_name)

  #Ordering the annotation file based on the order of the CpG beta dataframe
  annotation <- annotation[rownames(betas),]

  #Get a vector of CpGs which are not in the sexual chromosomes
  included_CpGs <- rownames(annotation)[!annotation[,arguments$column_chromosomes] %in% c("chrX", "chrY")]

  #Removing the CpGs belonging to the sexual chromosomes (the chromosomes that are not chrX a d chrY)
  betas <- betas[included_CpGs,]

}

# ========================================
# CREATING TRAINING AND VALIDATION SUBSETS
# ========================================

if (arguments$split_in_training_and_validation) {

cat("\n\nCreating training and validation datasets...\n")

#Setting the seed for spliting the data
set.seed(1)

#Splitting the data
split_ratio <- 1 - arguments$percentage_to_validation/100
train_samples <- sample.split(colnames(betas), SplitRatio=split_ratio)

#Creating training and validation datasets for the beta values
unadj_training <- betas[,train_samples]
unadj_validation <- betas[,!train_samples]

#Creating training and validation datasets for the purity_vector
purity_training <- purities[train_samples]
purity_validation <- purities[!train_samples]

}

# ========================================
# CREATING TRAINING AND VALIDATION SUBSETS
# ========================================

cat("\n\nSaving the data...\n")

if (arguments$split_in_training_and_validation) {

    saveRDS(unadj_training, paste(arguments$output_directory, "/betas_training.RData", sep=""))
    saveRDS(unadj_validation, paste(arguments$output_directory, "/betas_validation.RData", sep=""))

    saveRDS(purity_training, paste(arguments$output_directory, "/purity_training.RData", sep=""))
    saveRDS(purity_validation, paste(arguments$output_directory, "/purity_validation.RData", sep=""))

} else {

    saveRDS(arguments$betas, paste(arguments$output_directory, "/betas.RData", sep=""))
    saveRDS(arguments$purity, paste(arguments$output_directory, "/purity.RData", sep=""))

}

cat("\n\n****************\n")
cat("PROCESS FINISHED\n")
cat("****************\n")
