#!/usr/bin/Rscript

# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

#Loading optparse to parse command line arguments
if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))

#Loading catools to split the dataset in training and test
if(!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret") }

suppressPackageStartupMessages(library(caret))


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

  make_option(c("-u", "--purity"), type="character",  
              help="The name of the object that contains the beta values inside the R object must be entered here",
              metavar = "[var_name]"),

  make_option(c("-S", "--split_in_training_and_validation"), type="logical", default=FALSE,
              help="This argument must be set to true to specify if the data has to be splitted to create training and test datasets.",
              metavar = "[boolean]"),

  make_option(c("-v", "--percentage_to_validation"), type="double", default=20.0,
              help="The percentage of the samples to be included in the test data must be entered here. Default [%default]",
              metavar = "[number]"),

  make_option(c("-C", "--split_for_cross_validation"), type="logical", default=FALSE,
              help="This argument must be set to true to specify if the data has to be splitted to create different training and test datasets to perform a cross-test.",
              metavar = "[boolean]"),

  make_option(c("-k", "--groups_in_k_fold_cross_validation"), type="numeric",
              help="The number of groups in which the data will be splitted to perform the cross validation can be entered here",
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

  make_option(c("-N", "--get_only_CpGs_with_highest_variance"), type="logical", default=FALSE,
              help="This argument should be sent to true if the user wants to get only the betas of the CpGs with the highest variance",
              metavar="[TRUE/FALSE]"),

  make_option(c("-n", "--number_of_most_variable_CpGs_to_include"), type="numeric",
              help="The number of the most variant CpGs to be included in the output should be included here",
              metavar="[number]")       

)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This programm gets purity and original beta Robjects from a R object containing those or more elements. It also allow the user to cerate training anjd test dataset"))


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

#Adapting the names of the vector
names(purities) <- lapply(names(purities), 
                         function (name) {paste(name, "-01A", sep="")})


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


# ====================================
# CREATING GROUPS FOR CROSS VALIDATION
# ====================================

if (arguments$split_for_cross_validation) {

    cat("\n\nSplitting the data for the cross validation...\n")

    # Using the createFoldas from the caret package to create groups for the cross validation
    # The number of groups was set to 6 generate have around 100 samples per group
    groups <- createFolds(y=colnames(betas),
                          k=6,
                          list=TRUE,
                          returnTrain=TRUE)


    # Generating a list with the training and test data for each created fold
    folds_list <- lapply(X=names(groups), 
                         FUN=function(fold) {
                            
                            # Getting the training and test beta values
                            training <- betas[,groups[[fold]]]
                            test <- betas[,-groups[[fold]]]

                            # Getting the training and test purity values
                            purity_training <- purities[colnames(training)]
                            purity_test <- purities[colnames(test)]

                            # Adding the created dataframes to a list and returning it
                            return(list("Training"=training, 
                                        "Test"=test, 
                                        "Purity_training"=purity_training,
                                        "Purity_test"=purity_test))               
                            }
                        )

    # Using the name of ecah fold as the id of the elements of fold_list
    names(folds_list) <- names(groups)

}

# ========================================
# CREATING TRAINING AND TEST SUBSETS
# ========================================

#If the user has selected to filter sexual chromosomes the following code will be run
if (arguments$split_in_training_and_validation) {

    cat("\n\nCreating training and test datasets...\n")

    #Setting the seed for spliting the data
    set.seed(1)

    #Splitting the data
    split_ratio <- 1 - arguments$percentage_to_validation/100
    train_samples <- sample.split(colnames(betas), SplitRatio=split_ratio)

    #Creating training and test datasets for the beta values
    unadj_training <- betas[,train_samples]
    unadj_validation <- betas[,!train_samples]

    #Creating training and test datasets for the purity_vector
    purity_training <- purities[train_samples]
    purity_validation <- purities[!train_samples]
    
}


# =====================================================
# GETTING MOST VARIABLE CPGs BASED ON THE TRAINING DATA
# =====================================================

# If the user has selected to get only the most variable CpGs the following code will be run
if(arguments$get_only_CpGs_with_highest_variance) {
    
    if (arguments$split_in_training_and_validation) {

        # Detemine the variance of all the rows (CpGs)
        cpgs_variance <- apply(
                        unadj_training,
                        MARGIN=1,
                        FUN=var)

        # Create a vector to sort the rows. Get only the number of rows containing the CpGs that want to be included
        sorting_vec <- order(cpgs_variance)[1:arguments$number_of_most_variable_CpGs_to_include]

        # Sorting the betas dataframe and getting only the number of CpGs selected
        unadj_training <- unadj_training[sorting_vec,]
    
    } else if (arguments$split_for_cross_validation) {


        for (fold in folds_list) {

            # Detemine the variance of all the rows (CpGs)
            cpgs_variance <- apply(
                                fold$"Training",
                                MARGIN=1,
                                FUN=var
                        )

            # Create a vector to sort the rows. Get only the number of rows containing the CpGs that want to be included
            sorting_vec <- order(cpgs_variance)[1:arguments$number_of_most_variable_CpGs_to_include]

            # Getting only the cpgs selected in the training and test datasets
            fold$"Training" <- fold$"Training"[sorting_vec,]
            fold$"Test" <- fold$"Test"[sorting_vec,]

        }
    }
}


# ===========================================================
# CREATING OUTPUT FILES FOR EACH FOLD OF THE CROSS VALIDATION
# ===========================================================

if (arguments$split_in_training_and_validation) {

    cat("\n\nSaving output files...\n")

    #Saving beta values
    saveRDS(unadj_training, paste(arguments$output_directory, "/BetasTraining.RData", sep=""))
    saveRDS(unadj_validation, paste(arguments$output_directory, "/BetasTest.RData", sep=""))

    #Saving purity values
    saveRDS(purity_training, paste(arguments$output_directory, "/PurityTraining.RData", sep=""))
    saveRDS(purity_validation, paste(arguments$output_directory, "/PurityTest.RData", sep=""))

} else if (arguments$split_for_cross_validation) {

    cat("\n\nSaving output files...\n")

    for (fold in names(folds_list)) {

        # Saving training and test datasets
        saveRDS(folds_list[[fold]][["Training"]], paste(arguments$output_directory, "/", fold, "_BetasTraining.RData", sep=""))
        saveRDS(folds_list[[fold]][["Test"]], paste(arguments$output_directory, "/", fold, "_BetasTest.RData", sep=""))

        saveRDS(folds_list[[fold]][["Purity_training"]], paste(arguments$output_directory, "/", fold, "_PurityTraining.RData", sep=""))
        saveRDS(folds_list[[fold]][["Purity_test"]], paste(arguments$output_directory, "/", fold, "_PurityTest.RData", sep=""))
    }

}

cat("\n\n****************\n")
cat("PROCESS FINISHED\n")
cat("****************\n")