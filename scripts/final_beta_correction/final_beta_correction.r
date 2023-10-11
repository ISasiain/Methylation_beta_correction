#!/usr/bin/Rscript

## -SCRIPT'S NAME: final_beta_correction.r
#
## - DESCRIPTION: 
#
#   This script correct beta values of samples whose purity has been estimated based on
#   a reference cohort
# 
## - USED R PACKAGES:
#
#   *OPTPARSE. Parsing command line arguments
#   *DOPARALLEL. Parallelizing execution
#   *PARALLEL. Parallelizing script
#
## - USER DEFINED FUNCTIONS:
#   
#   *adjustBeta(). This function corrects beta values based on the sample purity generating
#    corrected betas for the actual cancer cells and the tumor microenvironment. It also outputs
#    the parameters of the regressions used for the correction
#
## - PROCEDURE:
#
#   1. Installing (if necessary) and loading packages, configuring command line arguments and sourcing 
#      user defined functions.
#
#   2. Configuring parallelization.
#
#   3. Loading the data, adding the seed to run the analysis with and running the adjustBeta() function
#      per each CpG using a parallelized apply function.
#
#   4. Adding the produced results to a result list.
#
#   5. Saving each element of the result list as an independent R object.
#
## - INPUT FILES:
#
#    -> Dataframe stored as an R object containig the original beta values of the cpgs and samples
#       to be corrected.
#
#    -> Named vector stored as an R object containing the purity values of the samples used for the analysis
#
#
## - OUTPUT FILES:
#
#    -> R object file containinga dataframe with the original beta values
#
#    -> R object file containing a dataframe with the corrected tumor beta values
#
#    -> R object file containing a dataframe with the corrected microenvironment beta values
#
#    -> R object file containing a dataframe with the slopes of the regressions used for the beta correction
#
#    -> R object containing a dataframe with the intercepts of the regressions used for the beta correction.
#
#    -> R object containing a dataframe with the Residual Standard Error of the regressions used for the beta correction.
#
#    -> R object containing a dataframe with the degrees of freedom of the regressions used for the beta correction.
#
#    -> R object containing a dataframe with the methylation patterns (populations) detected during the correction.
#
## - USAGE:
#
#     The script must be run on the command line using the following flags. 
#
#     """
#     Rscript path_to_script/new_purity_corrector.r -c [num_of_cores] -b [path_to_betas] 
#     -p [path_to_puritied] -o [path_to_save_output_files] -n [prefix_output_files]
#     """
#     
#     *The function of the command line options are the following; 
#
#       -c: Number of cores to be used to run the program. Default: 1.
#       -B: The path to the file with the beta values of the refernce cohort must be entered here. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.
#       -P: The path to the file with the purity values of the refrence cohort must be entered here. The file must be an R object containing a dictionary vector.
#       -b: Path to the file with the beta values to be corrected whose sample purity has been estimated. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.
#       -p: Path to the tsv file with the predicted sample purity values of the samples whose betas have to be corrected. The file must be the tsv text file generated as an output of run_all_validation.r.
#       -o: The path to the location where the output files will be saved must be entered here. The output is an R object. Default: working directory.
#       -n: The prefix to be used to name the output files. Default: output.
#
## - VERSION: 1.0
#
## - DATE: 17/05/2023
#
## - AUTHOR: Mattias Aine  (mattias.aine@med.lu.se)
## - ADAPTED BY: Iñaki Sasiain
## - AFFILIATION: Johan Staaf lab @ Lund University / Oncology & Pathology

# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

options(repos = "https://cran.r-project.org/")

if(!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel") }

suppressPackageStartupMessages(library(doParallel))

if(!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel") }

suppressPackageStartupMessages(library(parallel))

if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

argument_list <- list(

  make_option(c("-c", "--cores"), type="integer", default=1,  
              help="Number of cores to be used to run the program [default %default]",
              metavar = "[number]"),

  make_option(c("-P", "--ref_cohort_purity"), type="character",  
              help="Path to the file with the purity values of of the refernce cohort. The file must be an R object containing a dictionary vector.",
              metavar = "[file path]"),

  make_option(c("-B", "--ref_cohort_betas"), type="character",  
              help="Path to the file with the beta values of the refernce cohort. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.",
              metavar = "[file path]"),

  make_option(c("-p", "--est_purity"), type="character",  
              help="Path to the tsv file with the predicted sample purity values of the samples whose betas have to be corrected. The file must be the tsv text file generated as an output of run_all_validation.r.",
              metavar = "[file path]"),

  make_option(c("-b", "--betas_to_correct"), type="character",  
              help="Path to the file with the beta values to be corrected whose sample purity has been estimated. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.",
              metavar = "[file path]"),

  make_option(c("-o", "--output"), type="character", default="./",
              help="The path to the location where the output files will be saved must be entered here. The output is an R object. Default [%default]",
              metavar = "[file path]"),

  make_option(c("-n", "--output_name"), type="character", default="output",
              help="The prefix to be used to name the output files. Default [%default]",
              metavar = "[file path]")

)

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program corrects methylation beta values providing parameters of the regressions used for the correction."))



# =====================================
#   SOURCING THE CORRECT BETAS FUNCTION
# =====================================

dir <- commandArgs()[4]
dir <- gsub("--file=", "", dir)
dir <- gsub("new_purity_corrector.r", "new_function_correctBetas.r", dir)

source(dir)

# ===========================
# CONFIGURING PARALLELIZATION
# ===========================

cat("\nUsing", arguments$cores,"cores\n\n")

#Creating the cluster to run the process in parallel
cl <- makeCluster(arguments$cores)  
registerDoParallel(cl)  

#Making sure that all packages have access to the flexmix package. Using invisible()
#to avoid proning anything in the terminal

invisible(clusterEvalQ(cl, {library("flexmix")}))

# ================
# LOADING THE DATA
# ================

cat("\nLoading the data...\n\n")

# Loading refernce files (cohort)
cohort_betas <- readRDS(arguments$ref_cohort_betas)
cohort_purities <- readRDS(arguments$ref_cohort_purity)

# Loading betas to correct
to_correct_betas <- readRDS(arguments$betas_to_correct)

# Loading predicted purities
predicted_purities <- read.table(arguments$est_purity, 
                                header=TRUE,
                                sep="\t")

# Removing samples with more than one estimates (if any)
predicted_purities <- predicted_purities[which(predicted_purities[,"num_of_est"]==1),]

# Transforming the predicted_purities dataframe into a vector
predicted_purities_vec <- predicted_purities[,"estimate"]
names(predicted_purities_vec) <- predicted_purities[,"sample"]

# ================
# FILTERING CPGS
# ================

cat("\nChecking cpgs...\n\n")

if (length(!(rownames(to_correct_betas) %in% rownames(cohort_betas))) != 0) {

    # Printing warning message
    cat("\n",  length(!(rownames(to_correct_betas) %in% rownames(cohort_betas))), "CpGs are not included into the refernce cohort, so they can not be corrected.\n\n")

    # Filtering not included CpGs
    to_correct_betas <- to_correct_betas[rownames(to_correct_betas) %in% rownames(cohort_betas),]
}

# Remove CpGs from the cohort dataset that are not included into the data to correct to speed up the process.
cohort_betas <- cohort_betas[rownames(cohort_betas) %in% rownames(to_correct_betas),]

# ============
# MERGING DATA
# ============

# Creating a single purities vector
purities <- c(cohort_purities, predicted_purities_vec)

# Creating a single betas dataframe
betas <- merge(cohort_betas, to_correct_betas, by="row.name")

#Removing sample purities not contained into the beta dataset. It generates errors
purities <- purities[colnames(betas)]


# ====================================
# PREPROCESSING AND ANALYSING THE DATA
# ====================================


cat("\nRunning the correction...\n\n")

#Adding seed to each row of the beta value dataframe
betaRun <- cbind(seed=1:nrow(betas),betas)

#Storing sample names
betaNames <- colnames(betas)

#Running the analysis in parallel
res <- parRapply(cl = cl, #ClusterS to run the process
                 betaRun, #Beta values+the added seed
                 adjustBeta, #Function to correct betas
                 purity=purities, #Purity values
                 snames=betaNames, #Sample names
                 seed=TRUE) #The seed has been added in the data


# ====================
# CREATING RESULT LIST
# ====================

# Creating a list to add the results (Only of the samples to be corrected) ADAPT THIS!!!!
result_list <- list(

  betas.original = do.call("rbind",lapply(res,function(x) x$y.orig)), #Original beta values
  betas.tumor = do.call("rbind",lapply(res,function(x) x$y.tum)), #Corrected tumor beta values
  betas.microenvironment = do.call("rbind",lapply(res,function(x) x$y.norm)), #Corrected microenvironment beta values

)

# Creating a list to add the parameters of the correction regressions
reg_list <- list(

  cpg.populations =  do.call("rbind",lapply(res,function(x) x$groups)), #Methylation patterns (populations) of each CpG
  reg.slopes = do.call("rbind",lapply(res,function(x) x$res.slopes)), #Slopes of the populations
  reg.intercepts = do.call("rbind",lapply(res,function(x) x$res.int)), #Intercepts of the populations
  reg.RSE = do.call("rbind",lapply(res,function(x) x$res.rse)), #Residual standard error
  reg.df = do.call("rbind",lapply(res,function(x) x$res.df)) #Degrees of freedom of the reversed regressions

)

# =====================
# CREATING OUTPUT FILES
# =====================

#Defining a function to store the elements of the result list to RData files
df_to_RObj <- function(df, filename) {
  saveRDS(df, filename)
}

#Defining a function to store the elements of the result list as tsv files
df_to_tsv <- function(df,filename) {
  # Saving text file (tsv)
  write.table(df, 
            file=filename,
            col.names=NA, 
            sep="\t")
}

#Creating output files per each dataframe of the result_list list (getting only the results for the CpGs of the predicted samples)
lapply(names(result_list), function(n) {
  df_to_RObj(result_list[[n]][,names(predicted_purities_vec)],filename=paste(arguments$output, arguments$output_name,"_",n,".samples_to_correct.RData",sep=""))
  df_to_tsv(result_list[[n]][,names(predicted_purities_vec)],filename=paste(arguments$output, arguments$output_name,"_",n,".samples_to_correct.tsv",sep=""))
})

#Creating output files per each dataframe of the reg_list list
lapply(names(reg_list), function(n) {
  df_to_RObj(result_list[[n]],filename=paste(arguments$output, arguments$output_name,"_",n,".RData",sep=""))
})

# Stop clusters used in parallelization
stopCluster(cl)

#cat("\n=================\n")
#cat ("PROCESS FINISHED")
#cat("\n=================\n")