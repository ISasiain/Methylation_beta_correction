#!/usr/bin/Rscript

# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

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

  make_option(c("-d", "--path_to_directory"), type="character",  
              help="The path to the R object with the data to process must be entered here.",
              metavar = "[path]"),

  make_option(c("-b", "--betas"), type="character",  
              help="The name of the object that contains the beta values inside the R object must be entered here",
              metavar = "[var_name]"),

  make_option(c("-p", "--purity"), type="character",  
              help="The name of the object that contains the beta values inside the R object must be entered here",
              metavar = "[var_name]"),

  make_option(c("-s", "--split_in_training_and_validation"), type="logical", default=FALSE,
              help="This argument must be entered to specify if the data has to be splitted to create training and validation datasets.",
              metavar = "[boolean]"),

  make_option(c("-v", "--percentage_to_validation"), type="double", default=20.0,
              help="The percentage of the samples to be included in the validation data must be entered here.",
              metavar = "[number]"),

  make_option(c("-o", "--output_directory"), type="character", default=".",
              help="The directory in which the output files will be stored can be. The default value is the current directory",
              metavar = "[path]")

)

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This programm gets purity and original beta Robjects from a R object containing those or more elements. It also allow the user to cerate training anjd validation dataset"))


# ================
# LODING R OBJECT 
# ================

load(arguments$path_to_directory)

# ========================================
# CREATING TRAINING AND VALIDATION SUBSETS
# ========================================

if (arguments$split_in_training_and_validation) {

#Setting the seed for spliting the data
set.seed(1)

#Splitting the data
split_ratio <- 1 - arguments$percentage_to_validation/100
train_samples <- sample.split(colnames(get(arguments$betas)), SplitRatio=split_ratio)

#Creating training and validation datasets for the beta values
unadj_training <- get(arguments$betas)[,train_samples]
unadj_validation <- get(arguments$betas)[,!train_samples]

#Creating training and validation datasets for the purity_vector
purity_training <- get(arguments$purity)[train_samples]
purity_validation <- get(arguments$purity)[!train_samples]

}

# ========================================
# CREATING TRAINING AND VALIDATION SUBSETS
# ========================================

if (arguments$split_in_training_and_validation) {

    saveRDS(unadj_training, paste(arguments$output_directory, "/betas_training.RData", sep=""))
    saveRDS(unadj_validation, paste(arguments$output_directory, "/betas_validation.RData", sep=""))

    saveRDS(purity_training, paste(arguments$output_directory, "/purity_training.RData", sep=""))
    saveRDS(purity_validation, paste(arguments$output_directory, "/purity_validation.RData", sep=""))

} else {

    saveRDS(arguments$betas, paste(arguments$output_directory, "/betas.RData", sep=""))
    saveRDS(arguments$purity, paste(arguments$output_directory, "/purity.RData", sep=""))

}
