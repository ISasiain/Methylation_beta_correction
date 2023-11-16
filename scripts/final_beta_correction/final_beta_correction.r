#!/usr/bin/Rscript

## -SCRIPT'S NAME: final_beta_correction.r
#
## - DESCRIPTION: 
#
#   This script correct beta values of samples whose purity has been estimated based on
#   a reference cohort refitting the reference regressions to include both, the reference data 
#   points and the new ones (betas to correct + estimated purity).
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
#   3. Loading the data, if one of the samples has more than one purity estimates remove it sample, 
#      as no reliable correction is possible Adding the seed to run the analysis with and running 
#      the adjustBeta() function per each CpG using a parallelized apply function.
#
#   4. Filtering CpGs. The CpG to correct not included in the reference data set must be removed, as
#      they can not be corrected. Also, the CpGs in the reference data set not included in the data to
#      correct are removed to speed up the process.
#
#   5. Mergining reference and data to analyse to refit the regresisions
#
#   6. Generating results and adding the them to a result list.
#
#   7. Saving each element of the result list as an independent R object or TSV file.
#
## - INPUT FILES:
#
#    -> Dataframe stored as an R object containig the reference beta values.
#
#    -> Named vector stored as an R object containing the purity of the reference samples.
#
#    -> TSV file containing the purity estimates of the samples whose betas are to be corrected. This 
#       must be the output of purity_estimator.r
#
#    -> Dataframe stored as an R object containig the beta values with the betas of the samples to correct.
#
#    -> Vector of CpGs to correct stored as an R object. This must be only provided if the user chooses to limit
#       the beta correction to certain CpGs.
#
#
## - OUTPUT FILES:
#
#    -> R object file containing a dataframe with the original beta values
#
#    -> TSV file containing the original beta values of each sample
#
#    -> R object file containing a dataframe with the corrected tumor beta values
#
#    -> TSV file containing the corrected cancer cell beta values of each sample
#
#    -> R object file containing a dataframe with the corrected microenvironment beta values
#
#    -> TSV file containing the corrected microenvironment beta values of each sample
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
#     Rscript path_to_script/new_purity_corrector.r -c [num_of_cores] -B [path_to_ref_betas] 
#     -P [path_to_ref:purities] -b [path_to_betas_to_correct] -p [path_to_estimated purities]
#     -F [correct_only_certain_CpGs: TRUE/FALSE] -f [vec_CpGs_to_correct]
#     -o [path_to_save_output_files] -n [prefix_output_files]
#     """
#     
#     *The function of the command line options are the following; 
#
#       -c: Number of cores to be used to run the program. Default: 1.
#       -B: The path to the file with the beta values of the reference cohort must be entered here. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.
#       -P: The path to the file with the purity values of the refrence cohort must be entered here. The file must be an R object containing a dictionary vector.
#       -b: Path to the file with the beta values to be corrected whose sample purity has been estimated. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.
#       -p: Path to the tsv file with the predicted sample purity values of the samples whose betas have to be corrected. The file must be the tsv text file generated as an output of run_all_validation.r.
#       -F: This argument should be set TRUE if a list with the CpGs to correct wants to be provided.
#       -f: The path of the Robject containing a vector with the CpGs to correct should be entered here.
#       -o: The path to the location where the output files will be saved must be entered here. The output is an R object. Default: working directory.
#       -n: The prefix to be used to name the output files. Default: output.
#
## - VERSION: 1.0
#
## - DATE: 30/10/2023
#
## - AUTHOR: Iñaki Sasiain
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

  make_option(c("-F", "--only_certain_CpGs"), type="logical", default=FALSE,   
              help="This argument should be set TRUE if a list with the CpGs to correct wants to be provided.",
              metavar = "[TRUE/FALSE]"),

  make_option(c("-f", "--CpGs_to_correct_vec"), type="character",  
              help="The path of the Robject containing a vector with the CpGs to correct should be entered here.",
              metavar = "[file path]"),

  make_option(c("-o", "--output"), type="character", default="./",
              help="The path to the location where the output files will be saved must be entered here. The output is an R object. Default [%default]",
              metavar = "[file path]"),

  make_option(c("-n", "--output_name"), type="character", default="output",
              help="The prefix to be used to name the output files. Default [%default]",
              metavar = "[file path]"),

  make_option(c("-r", "--refitting"), type="logical",
              help="If the user wants to refit the refrence regressions used for the beta correction this argument should be set to TRUE",
              metavar="[TRUE/FALSE]"),

  make_option(c("-R", "--ref_regressions"), type="character",  
              help="Path of the directory containing the RObjects with the parameters of the refernce regressions",
              metavar = "[file path]"),

  make_option(c("-P", "--ref_cohort_purity"), type="character",  
              help="Path to the file with the purity values of of the reference cohort. The file must be an R object containing a dictionary vector.",
              metavar = "[file path]"),

  make_option(c("-B", "--ref_cohort_betas"), type="character",  
              help="Path to the file with the beta values of the reference cohort. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.",
              metavar = "[file path]"),

  make_option(c("-p", "--est_purity"), type="character",  
              help="Path to the tsv file with the predicted sample purity values of the samples whose betas have to be corrected. The file must be the tsv text file generated as an output of run_all_validation.r.",
              metavar = "[file path]"),

  make_option(c("-b", "--betas_to_correct"), type="character",  
              help="Path to the file with the beta values to be corrected whose sample purity has been estimated. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.",
              metavar = "[file path]")
)

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program corrects methylation beta values providing parameters of the regressions used for the correction."))

if (arguments$refitting == FALSE) {

# =====================================
#   SOURCING THE CORRECT BETAS FUNCTION
# =====================================

dir <- commandArgs()[4]

dir <- gsub("--file=", "", dir)
dir <- gsub("final_beta_correction.r", "new_function_correctBetas.r", dir)

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

# Loading reference files (cohort)
cohort_betas <- readRDS(arguments$ref_cohort_betas)
cohort_purities <- readRDS(arguments$ref_cohort_purity)

# Loading betas to correct
to_correct_betas <- readRDS(arguments$betas_to_correct)

# Loading predicted purities
predicted_purities <- read.table(arguments$est_purity, 
                                sep="\t")


# Removing samples with more than one estimates (if any)
predicted_purities <- predicted_purities[which(predicted_purities[,2]==1),]

# Transforming the predicted_purities dataframe into a vector
predicted_purities_vec <- 1-predicted_purities[,3]
names(predicted_purities_vec) <- predicted_purities[,1]

# ================
# FILTERING CPGS
# ================

cat("\nChecking cpgs...\n\n")
# Use only the specified CpGs if that option has been selected
if (arguments$only_certain_CpGs) {

  #Getting the vector
  vec_of_cpgs <- readRDS(arguments$CpGs_to_correct_vec)

  #Keeping only CpGs of interest
  to_correct_betas <- to_correct_betas[vec_of_cpgs,]

}

# Checking if the CpGs are included in the reference data
if (sum(!(rownames(to_correct_betas) %in% rownames(cohort_betas))) != 0) {

    # Printing warning message
    cat("\n",  sum(!(rownames(to_correct_betas) %in% rownames(cohort_betas))), "CpG(s) is/are not included into the reference cohort, so it/they can not be corrected.\n\n")

    # Filtering not included CpGs
    to_correct_betas <- to_correct_betas[rownames(to_correct_betas) %in% rownames(cohort_betas),]
}

# Remove CpGs from the cohort dataset that are not included into the data to correct to speed up the process.
cohort_betas <- cohort_betas[rownames(cohort_betas) %in% rownames(to_correct_betas),]

#Sorting the cohort betas dataframe based on the rownames of to_correct_betas
cohort_betas <- cohort_betas[rownames(to_correct_betas),]


# ============
# MERGING DATA
# ============

# Creating a single purities vector
purities <- c(cohort_purities, predicted_purities_vec)

# Creating a single betas dataframe
betas <- cbind(cohort_betas, to_correct_betas)

#Adapt name of betas
colnames(betas) <- lapply(colnames(betas), 
                         function (name) {strsplit(name, "-01")[[1]][1]})

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
  betas.microenvironment = do.call("rbind",lapply(res,function(x) x$y.norm)) #Corrected microenvironment beta values
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
  saveRDS(df, file=filename)
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
  df_to_RObj(reg_list[[n]],filename=paste(arguments$output, arguments$output_name,"_",n,".RData",sep=""))
})

# Stop clusters used in parallelization
stopCluster(cl)

cat("\n=================\n")
cat ("PROCESS FINISHED")
cat("\n=================\n")

} else {

# ============
# LOADING DATA
# ============

cat("\nLoading the data...\n\n")

#Reading the R objects containing the regression data as dataframes
my_slopes <- readRDS(list.files(arguments$ref_regressions, pattern="*reg.slopes.RData", full.names=TRUE))
my_intercepts <- readRDS(list.files(arguments$ref_regressions, pattern="*reg.intercepts.RData", full.names=TRUE))


# Loading betas to correct
to_correct_betas <- readRDS(arguments$betas_to_correct)


# Loading predicted 1 - purities
predicted_1mPurities <- read.table(arguments$est_purity, 
                                sep="\t")


# Removing samples with more than one estimates (if any)
predicted_1mPurities <- predicted_1mPurities[which(predicted_1mPurities[,2]==1),]

# Transforming the predicted_1mPurities dataframe into a vector
predicted_1mPurities_vec <- predicted_1mPurities[,3]
names(predicted_1mPurities_vec) <- predicted_1mPurities[,1]


# ================
# FILTERING CPGS
# ================

cat("\nChecking cpgs...\n\n")

# Use only the specified CpGs if that option has been selected
if (arguments$only_certain_CpGs) {

  #Getting the vector
  vec_of_cpgs <- readRDS(arguments$CpGs_to_correct_vec)

  #Keeping only CpGs of interest
  to_correct_betas <- to_correct_betas[vec_of_cpgs,]

}

if (sum(!(rownames(to_correct_betas) %in% rownames(my_slopes))) != 0) {

    # Printing warning message
    cat("\n",  sum(!(rownames(to_correct_betas) %in% rownames(my_slopes))), "CpG(s) is/are not included into the refernce cohort, so it/they can not be corrected.\n\n")
    
    # Filtering not included CpGs
    to_correct_betas <- to_correct_betas[rownames(to_correct_betas) %in% rownames(my_slopes),]
}

# Remove CpGs from the regressions that are not included into the data to correct to speed up the process.
my_slopes <- my_slopes[rownames(my_slopes) %in% rownames(to_correct_betas),]
my_intercepts <- my_intercepts[rownames(my_intercepts) %in% rownames(to_correct_betas),]

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


# ==============================================
# CORRECTING BETAS BASED ON REFERNCE REGRESSIONS
# ==============================================

cat("\nCorrecting betas...\n\n")

#Generating function to identify the refernce regression to which each CpG of each sample belongs
identify_regression <- function(beta, estimated_1mPurity, vec_slopes, vec_intercepts) {


    # Generating a vector of distances
    vec_of_dis <- unname(sapply(c(1:length(vec_slopes)), 
                         FUN = function(pop) {

                                    #Ignoring NA values that appear when less than 3 populations are detected
                                    #in the refernce rrgressions
                                    if (!is.na(vec_slopes[pop])) {
                                        beta - (vec_slopes[pop] * estimated_1mPurity + vec_intercepts[pop])

                                    } else {
                                        NA

                                    }
                                }
                            ))


    #Determining the population (vector index) with the lowest absolute distance. If the distances are equal the first
    #population with be chosen by default
    pop_identified <- which.min(abs(vec_of_dis))

    #Generating and returning output named vector
    return(
        c("Slope"=vec_slopes[pop_identified],
          "Intercept"=vec_intercepts[pop_identified],
          "Distance"=vec_of_dis[pop_identified])
        )
    }


#Generate function to correct betas based on the identified regression. The parametres specified must be 
#from the beta VS 1-P regressions.
correcting_betas <- function(slope, intercept, distance, to_correct) {


    if (to_correct=="Tumor") {
        #The tumor beta value will be obtained using the intercept and the calculated distance.
        #The maximum possible value will allways be kept below or equal to 1
        tum_beta <- if(intercept + distance <= 1) {
                        intercept + distance
                    } else {
                        1
                    }

        return(tum_beta)


    } else if (to_correct=="Microenvironment") {
        #The microenvironment beta value will be obtained using the intercept and slope when 1-P=1 and the calculated distance.
        #The minimum possible value will allways be kept below or equal to 1
        env_beta <- if (intercept + slope + distance >= 0) {
                        intercept + slope + distance
                    } else {
                        0
                    }
        
        return(env_beta)

    }
}

 # Creating dataframes to store the corrected values
 corrected_tumor <- data.frame()
 corrected_microenvironment <- data.frame()

 # Correcting betas through a parallelized for loop
 foreach(cpg=rownames(to_correct_betas)) %dopar% {

    corrected_tumor[cpg,colnames(to_correct_betas)] <- sapply(names(predicted_1mPurities_vec), 

        function(sample) {

        identified_reg <- identify_regression(
            beta=to_correct_betas[cpg, sample],
            estimated_1mPurity=predicted_1mPurities_vec[sample],
            vec_slopes=my_slopes[cpg,],
            vec_intercepts=my_intercepts[cpg,]
            )


        corrected_tum <- correcting_betas(
            slope=identified_reg["Slope"],
            intercept=identified_reg["Intercept"],
            distance=identified_reg["Distance"],
            to_correct="Tumor"
            )


        return(corrected_tum)

        }
    )

    corrected_microenvironment[cpg,colnames(to_correct_betas)] <- sapply(names(predicted_1mPurities_vec), 

        function(sample) {

        identified_reg <- identify_regression(
            beta=to_correct_betas[cpg, sample],
            estimated_1mPurity=predicted_1mPurities_vec[sample],
            vec_slopes=my_slopes[cpg,],
            vec_intercepts=my_intercepts[cpg,]
            )


        corrected_env <- correcting_betas(
            slope=identified_reg["Slope"],
            intercept=identified_reg["Intercept"],
            distance=identified_reg["Distance"],
            to_correct="Microenvironment"
            )

        return(corrected_env)

        }
    )
 }

 # Stop clusters used in parallelization
 stopCluster(cl)

# =======================
# GENERATING OUTPUT FILES
# =======================

cat("\nGenerating output files...\n\n")

# Generating RObject files
saveRDS(corrected_tumor, file=paste(arguments$output, arguments$output_name, ".tumor.samples_to_correct.RData", sep=""))
saveRDS(corrected_microenvironment, file=paste(arguments$output, arguments$output_name, ".microenvironment.samples_to_correct.RData", sep=""))
saveRDS(to_correct_betas, file=paste(arguments$output, arguments$output_name, ".original.samples_to_correct.RData", sep=""))


# Generating tsv files
write.table(corrected_tumor, 
            file=paste(arguments$output, arguments$output_name, ".tumor.samples_to_correct.tsv", sep=""),
            col.names=NA, 
            sep="\t")
write.table(corrected_microenvironment, 
            file=paste(arguments$output, arguments$output_name, ".microenvironment.samples_to_correct.tsv", sep=""),
            col.names=NA, 
            sep="\t")
write.table(to_correct_betas, 
            file=paste(arguments$output, arguments$output_name, ".original.samples_to_correct.tsv", sep=""),
            col.names=NA, 
            sep="\t")

cat("\n=================\n")
cat ("PROCESS FINISHED")
cat("\n=================\n")

}