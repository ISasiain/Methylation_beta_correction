#!/usr/bin/Rscript

## -SCRIPT'S NAME: preprocessing_data.r
#
## - DESCRIPTION: 
#
#   This script preprocesses methylation beta values of different samples, allowing to split
# the data in training and validation datasets, filter cpgs from sexual chromosomes, split
# the cpgs based on the infinium design and get only the number of the most variable cpgs.
# 
## - USED R PACKAGES:
#
#   *OPTPARSE. Parsing command line arguments
#   *CATOOLS. Splitting dataset in training and validation
#
## - USER DEFINED FUNCTIONS:
#   
#   None
#
## - PROCEDURE:
#
#   1. Installing (if necessary) and loading packages and configuring command line arguments.
#
#   2. Loading the data. A single R object will be loaded if the betas and purities are stored
#      in the same R object, if they weren't, the two R objects would be loaded independently.
# 
#   3. Filter cpgs from sexual chromosomes if the user had selected to do so. Only the cpgs not 
#      belonging to ChrX or ChrY will be kept in the dataframe.
#
#   4. Keepeing the specified number of most variable cpgs if the user had selected to do so. The
#      remaining cpgs will be deleted from the dataframe.
#
#   5. Dividing dataset in function of the infinium design if the user has selected to do so. The
#      dataframe containing the beta values will be splitted in two, one for each infinium design.
#
#   6. Creating training and validation datasets based on the specified proportion if the user has
#      selected to do so. If the beta dataset had previously been splitted regarding the infinium 
#      design one training and validation dataset will be generated for each infinium type.
#
#   7. Saving filtered and splitted datasets as R objects. The number of produced files will depend
#      on the times in which the beta values have been splitted (infinium design and 
#      training/validation).
#
## - INPUT FILES:
#
#     -> File stored as an R object containing the beta values used for the purity prediction. This
#        file can be an independen .RData file or included in the same object than the purity values.
#
#     -> File stored as an R object containing the actual purity values of the predicted samples. This
#        file can be an independen .RData file or included in the same object than the beta values.
#
#     -> If the option to filter sexual chromosomes has been activated an annotation file containing the
#        chromosomes to which each cpg belong should be provided stored as an R object.
#
#     -> If the option to split the data based on the infinium technique has been activated an annotation 
#        file containing the design used for the determination of the beta value of each cpg should be
#        entered as an R object.
#
## - OUTPUT FILES:
#
#    The number of output files will depend on the options chosen;
#
#     * If the user has not selected to split the data on training and validation datasets and in function 
#       of the infinium design;
#
#       -> An R object file containing a dataframe with the betas of each sample and cpg.
#       -> An R object file containing a named vector with the actual purity values of each sample.
#
#     * If the user has selected to split the data on training and validation datasets but not in function 
#       of the infinium design;
#
#       -> Two R object files, training and validation, containing dataframe with the betas of ecah sample and cpg. 
#       -> Two R object files, training and validation, containing named vectors with the actual purity values of each sample.
#
#     * If the user has selected not to split the data on training and validation datasets but splitting it in function 
#       of the infinium design;
#
#       -> Two R object files, with cpgs determined using infiniumI and infiniumII, containing dataframe with the betas of ecah sample and cpg. 
#       -> An R object file containing a named vector with the actual purity values of each sample.
#
#     * If the user has selected to split the data on training and validation datasets qnd also in function 
#       of the infinium design;
#
#       -> Four R object files, training and validation with cpgs determined using infiniumI and infiniumII, containing dataframe with the betas
#          of ecah sample and cpg. 
#       -> Two R object files, training and validation, containing named vectors with the actual purity values of each sample.
#
## - USAGE:
#
#     The script must be run on the command line using the following flags. 
#
#     """
#     Rscript path_to_script/preprocessing_data.r -s [TRUE/FALSE] -B [path_to_betas] -P [path_to_purities] -b [betas_name] -p [purities_name]
#     -S [TRUE/FALSE] -v [float/int] -o [path_to_output] -f [TRUE/FALSE] -a [path_to_annotation_file] -a [name_of_annotation_df] -c [chr_column_name]
#     -d [TRUE/FALSE] -D [path_to_design_file] -N [TRUE/FALSE] -n [integer]
#     """
#     
#     *The function of the command line options are the following; 
#
#       -a: The user must especify if the purity and betas are in the same R object.
#       -B: The path to the R object with the data to process must be entered here. If the purities and betas are in different R object files the one containing the betas should be entered here
#       -P: The path to the R object with the purities must be entered here. This should only be used if the purities and betas are in different R objects.
#       -b: The name of the object that contains the beta values inside the R object must be entered here
#       -p: The name of the object that contains the beta values inside the R object must be entered here
#       -S: This argument must be entered to specify if the data has to be splitted to create training and validation datasets.
#       -v: The percentage of the samples to be included in the validation data must be entered here. Default = 20 %.
#       -o: The directory in which the output files will be stored can be. The default value is the current directory.
#       -f: Set this argument to TRUE to filter CpGs from sexual chromosomes
#       -A: Specify the name of the path containing the annotation file.
#       -a: Specify the name of the dataframe containing the annotation information.
#       -c: Select the name of column of the dataframe containig the information about the chromosomes.
#       -d: If the user wants to split the data based on the infinium design (infinium I or II) this argument should be set to TRUE
#       -D: he path to the r object containing a list with the CpGs quantified with each infinium method.
#       -N: This argument should be sent to true if the user wants to get only the betas of the CpGs with the highest variance
#       -n: The number of the most variant CpGs to be included in the output should be included here
#
## - VERSION: 1.0
#
## - DATE: 18/05/2023
#
## - AUTHOR: Iñaki Sasiain Casado


# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

#Loading optparse to parse command line arguments
if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))

#Loading catools to split the dataset in training and validation
if(!requireNamespace("caTools", quietly = TRUE)) {
  install.packages("caTools") }

suppressPackageStartupMessages(library(caTools))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

#Creating argument list
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
              help="The percentage of the samples to be included in the validation data must be entered here. Default [%default]",
              metavar = "[number]"),

  make_option(c("-o", "--output_directory"), type = "character", default =".",
              help = "The directory in which the output files will be stored can be. The default value is the current directory",
              metavar ="[path]"),

  make_option(c("-f", "--filter_non_autosomes"), type="logical", default=FALSE,
              help="Set this argument to TRUE to discard CpGs from sexual chromosomes", 
              metavar="[TRUE/FALSE]"),

  make_option(c("-A", "--path_to_annotation"), type="character",
              help="Specify the name of the path containing the annotation file.",
              metavar="[path]"),

  make_option(c("-a", "--annotation_dataframe_name"), type="character",
              help="Specify the name of the dataframe containing the annotation information.",
              metavar="[variable name]"),

  make_option(c("-c", "--column_chromosomes"), type="character",
              help="Select the name of column of the dataframe containig the information about the chromosomes.",
              metavar="[column name]"),

  make_option(c("-d", "--split_based_on_design"), type="logical", default=FALSE,
              help="If the user wants to split the data based on the beta value detremination methods (infinium I or II) this argument should be set to TRUE",
              metavar="[TRUE/FALSE]"),

  make_option(c("-D", "--path_to_design_list"), type="character",
              help="The path to the r object containing a list with the CpGs quantified with each infinium method.",
              metavar="[path]"),

  make_option(c("-N", "--get_only_CpGs_with_highest_variance"), type="logical", default=FALSE,
              help="This argument should be sent to true if the user wants to get only the betas of the CpGs with the highest variance",
              metavar="[TRUE/FALSE]"),

  make_option(c("-n", "--number_of_most_variable_CpGs_to_include"), type="numeric",
              help="The number of the most variant CpGs to be included in the output should be included here",
              metavar="[number]")       

)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This programm gets purity and original beta Robjects from a R object containing those or more elements. It also allow the user to cerate training anjd validation dataset"))


# ===========
# LODING DATA
# ===========

cat("\n\nLoading the data...\n")

#If the betas and purities are inside the same R object load a single file
if (arguments$in_same_object) { 

  load(arguments$path_to_directory) #Loding the R object with the purities and betas

#If the betas and purities are not inside the same R object load two files
} else {

  load(arguments$path_to_directory) #Loding the R object with the betas
  load(arguments$path_to_purities) #Loding the R object with the purities

}

# Getting betas and purities from the loaded R objects
betas <- get(arguments$betas)
purities <- get(arguments$purity)

# ======================================
# FILTERING CpGs FROM SEXUAL CHROMOSOMES
# ======================================

# If the user has selected to filter sexual chromosomes the following code will be run
if (arguments$filter_non_autosomes) {
  
  cat("\n\nFiltering CpGs from sexual chromosomes...\n")

  #Loading annotation file with the chromosome info and getting the data
  load(arguments$path_to_annotation) 
  annotation <- get(arguments$annotation_dataframe_name)

  #Ordering the annotation file based on the order of the CpG beta dataframe
  annotation <- annotation[rownames(betas),]

  #Get a vector of CpGs which are not in the sexual chromosomes
  included_CpGs <- rownames(annotation)[!annotation[,arguments$column_chromosomes] %in% c("chrX", "chrY")]

  #Removing the CpGs belonging to the sexual chromosomes (the chromosomes that are not chrX a d chrY)
  betas <- betas[included_CpGs,]

}


# =================================
#    GETTING MOST VARIABLE CPGs
# =================================

# If the user has selected to get only the most variable CpGs the following code will be run
if(arguments$get_only_CpGs_with_highest_variance) {

  # Detemine the variance of all the rows (CpGs)
  cpgs_variance <- apply(
                        betas,
                        MARGIN=1,
                        FUN=var
                   )

  # Create a vector to sort the rows. Get only the number of rows containing the CpGs that want to be included
  sorting_vec <- order(cpgs_variance, decreasing=TRUE)[1:arguments$number_of_most_variable_CpGs_to_include]

  # Sorting the betas dataframe and getting only the number of CpGs selected
  betas <- betas[sorting_vec,]
}


# ===================================================
# DIVIDING DATASET IN FUNCTION OF THE INFINIUM METHOD
# ===================================================

if (arguments$split_based_on_design) {

  designs <- readRDS(arguments$path_to_design_list)


  betas_I <- betas[rownames(betas) %in% designs[["I"]],]
  betas_II <- betas[rownames(betas) %in% designs[["II"]],]

}

# ========================================
# CREATING TRAINING AND VALIDATION SUBSETS
# ========================================

#If the user has selected to filter sexual chromosomes the following code will be run
if (arguments$split_in_training_and_validation) {

  if (!arguments$split_based_on_design) {
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
    
  } else {

    cat("\n\nCreating training and validation datasets for each of the beta estimation method...\n")

    #Setting the seed for spliting the data
    set.seed(1)

    #Splitting the data
    split_ratio <- 1 - arguments$percentage_to_validation/100
    train_samples <- sample.split(colnames(betas_I), SplitRatio=split_ratio)
    

    #Creating training and validation datasets for the beta values
    unadj_training_I <- betas_I[,train_samples]
    unadj_training_II <- betas_II[,train_samples]
    unadj_validation_I <- betas_I[,!train_samples]
    unadj_validation_II <- betas_II[,!train_samples]

    #Creating training and validation datasets for the purity_vector
    purity_training <- purities[train_samples]
    purity_validation <- purities[!train_samples]

  }
}

# =======================
# SAVING OUTPUT R OBJECTS
# =======================

cat("\n\nSaving the data...\n")

#Run if the user has selected to split the data in training and validation subsets
if (arguments$split_in_training_and_validation) {

   if (!arguments$split_based_on_design) {  

    #Saving beta values
    saveRDS(unadj_training, paste(arguments$output_directory, "/betas_training.RData", sep=""))
    saveRDS(unadj_validation, paste(arguments$output_directory, "/betas_validation.RData", sep=""))

    #Saving purity values
    saveRDS(purity_training, paste(arguments$output_directory, "/purity_training.RData", sep=""))
    saveRDS(purity_validation, paste(arguments$output_directory, "/purity_validation.RData", sep=""))

   } else {

    #Saving beta values
    saveRDS(unadj_training_I, paste(arguments$output_directory, "/betas_training_I.RData", sep=""))
    saveRDS(unadj_training_II, paste(arguments$output_directory, "/betas_training_II.RData", sep=""))
    saveRDS(unadj_validation_I, paste(arguments$output_directory, "/betas_validation_I.RData", sep=""))
    saveRDS(unadj_validation_II, paste(arguments$output_directory, "/betas_validation_II.RData", sep=""))

    #Saving purity value
    saveRDS(purity_training, paste(arguments$output_directory, "/purity_training.RData", sep=""))
    saveRDS(purity_validation, paste(arguments$output_directory, "/purity_validation.RData", sep=""))

   }

#Run if the user has not selected to split the data in training and validation subsets
} else {

  if (!arguments$split_based_on_design) {  

    #Saving beta values
    saveRDS(arguments$betas, paste(arguments$output_directory, "/betas.RData", sep=""))
    
    #Saving purity values
    saveRDS(arguments$purity, paste(arguments$output_directory, "/purity.RData", sep=""))

  } else {

    #Saving beta values
    saveRDS(arguments$betas_I, paste(arguments$output_directory, "/betas_I.RData", sep=""))
    saveRDS(arguments$betas_I, paste(arguments$output_directory, "/betas_II.RData", sep=""))
    
    #Saving purity values
    saveRDS(arguments$purity, paste(arguments$output_directory, "/purity.RData", sep=""))

  }

}

cat("\n\n****************\n")
cat("PROCESS FINISHED\n")
cat("****************\n")
