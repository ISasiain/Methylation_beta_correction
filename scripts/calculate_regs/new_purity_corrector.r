#!/usr/bin/Rscript

#####======================================================================#####
###    Methylation beta value corrector based on the actual sample purity    ###
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Adapted by: Iñaki Sasiain
##Affiliation: Johan Staaf lab @ Lund University / Oncology & Pathology



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

  make_option(c("-b", "--input_betas"), type="character",  
              help="The path to the file with the beta values to be analysed must be entered here. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.",
              metavar = "[file path]"),

  make_option(c("-p", "--input_purity"), type="character",  
              help="The path to the file with the purity values of the samples to be analysed must be entered here. The file must be an R object containing a dictionary vector.",
              metavar = "[file path]"),

  make_option(c("-o", "--output"), type="character", default="./",
              help="The path to the location where the output file will be saved must be entered here. The output is an R object. Default [%default]",
              metavar = "[file path]"),

  make_option(c("-n", "--output_name"), type="character", default="output",
              help="The prefix to be used to name the output files. Default [%default]",
              metavar = "[file path]")

)

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program corrects methylation beta values providing parameters of the regressions used for the correction."))



# =============================
# SOURCING THE FUNCTION
# =============================

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



# =============================================
# LOADING, PREPROCESSING AND ANALYSING THE DATA
# =============================================

unadjusted_betas <- readRDS(arguments$input_beta)
purities <- readRDS(arguments$input_purity)

#Adding seed to each row of the beta value dataframe
betaRun <- cbind(seed=1:nrow(unadjusted_betas),unadjusted_betas)

#Storing sample names
betaNames <- colnames(unadjusted_betas)

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

result_list <- list(

  betas.original = do.call("rbind",lapply(res,function(x) x$y.orig)), #Original beta values
  betas.tumor = do.call("rbind",lapply(res,function(x) x$y.tum)), #Corrected tumor beta values
  betas.microenvironment = do.call("rbind",lapply(res,function(x) x$y.norm)), #Corrected microenvironment beta values
  cpg.populations =  do.call("rbind",lapply(res,function(x) x$groups)), #Methylation patterns (populations) of each CpG
  reg.slopes = do.call("rbind",lapply(res,function(x) x$res.slopes)), #Slopes of the populations
  reg.intercepts = do.call("rbind",lapply(res,function(x) x$res.int)), #Intercepts of the populations
  reg.RSE = do.call("rbind",lapply(res,function(x) x$res.rse)), #Residual standard error
  reg.df = do.call("rbind",lapply(res,function(x) x$res.df)) #Degrees of freedom of the reversed regressions

)

# =====================
# CREATING OUTPUT FILES
# =====================

##Creating a prefix for the output files
out_prefix <- arguments$output_name

#Defining a function to export the dataframes as csv files
df_to_RObj <- function(df, filename) {
  saveRDS(df, filename)
}

#Creating output files per each dataframe of the output_list list
lapply(names(result_list), function(n) {
  df_to_RObj(result_list[[n]],filename=paste(arguments$output, out_prefix,"_",n,".RData",sep=""))
})

stopCluster(cl)

#cat("\n=================\n")
#cat ("PROCESS FINISHED")
#cat("\n=================\n")