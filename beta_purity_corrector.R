##METHYLATION PURITY CORRECTOR (adapted from Mattias Aines script)
# ===============================================
# SETTING WORKING DIRECTORY AND LOCATION OF FILES
# ===============================================

##set/create own home directory below:
setwd("~/Methylation/adjustBetas")
HOME <- "~/Methylation/adjustBetas"

# ============================================================
# LOADING THE REQUIRED PACKAGES AND SOURCING THE REFERNCE DATA
# ============================================================

##load required packages

if(!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel") }

library(doParallel)

if(!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel") }

library(parallel)

if(!requireNamespace("caTools", quietly = TRUE)) {
  install.packages("caTools") }

library(caTools)

#Loading and sourcing R data and functions

##data set reduced to top 5000 CpGs in set of 630 TCGA BRCA tumors
load("workspace_tcgaBrca_top5000.RData")

##source - flexmix loaded on source
source(paste0(HOME,"/function_correctBetas.r"))

# ========================================
# CREATING TRAINING AND VALIDATION SUBSETS
# ========================================

#Setting the seed for spliting the data
set.seed(1)

#Splitting the data
train_samples <- sample.split(colnames(betaUnadjusted), SplitRatio=0.8)

#Creating training and validation datasets for the beta values
unadj_training <- betaUnadjusted[,train_samples]
unadj_validation <- betaUnadjusted[,!train_samples]

#Creating training and validation datasets for the purity_vector
purity_training <- purityVector[train_samples]
purity_validation <- purityVector[!train_samples]

# ============================================
# CORRECTING BETAS FOR THE TRAINING POPULATION
# ============================================


#need to feed seed to function so fully reproducible

no_cores <- detectCores(logical = TRUE)

cat("Using", no_cores-1,"cores","\n")

cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)  

##estimated runtime is ~500 sec on 7 cores

clusterEvalQ(cl, {
  library("flexmix")
})

##add rng seed to each row and pass to each paralell instance
betaRun<-cbind(seed=1:nrow(unadj_training),unadj_training)
betaNames<-colnames(unadj_training)

#clusterSetRNGStream(cl, 20200918) ##using this way to pass seed will not make exactly replicable..
res<-parRapply(cl = cl, betaRun, adjustBeta,purity=purity_training,snames=betaNames,seed=TRUE)

# =======================
# STORING THE OUTPUT DATA
# =======================

##getting parameters of the created 
temp1<-do.call("rbind",lapply(res,function(x) x$y.tum)) #Corrected tumor betas
temp2<-do.call("rbind",lapply(res,function(x) x$y.norm)) #Corrected microenvironment betas
temp3<-do.call("rbind",lapply(res,function(x) x$y.orig)) #Original data
temp4 <- do.call("rbind",lapply(res,function(x) x$groups)) #Clustered populations per each CpG
temp5 <- do.call("rbind",lapply(res,function(x) x$model.slopes)) #Slopes of the models
temp6 <- do.call("rbind",lapply(res,function(x) x$model.intercepts)) #Intercepts of the models
temp7 <- do.call("rbind",lapply(res,function(x) x$model.RSE)) #Mean of the residual's absolute value per population
temp8 <- do.call("rbind",lapply(res,function(x) x$model.SSx)) #Sum of squares of the x axis
temp9 <- do.call("rbind",lapply(res,function(x) x$model.SSy)) #Sum of squares of the y axis
temp10 <- do.call("rbind",lapply(res,function(x) x$model.df)) #Mean of the residual's absolute value per population

#Adding column names to temp4
colnames(temp4) <- betaNames

#Storing the software's output as a list of dataframes
output_list <- list(
  betas.original=temp3, #Original beta values
  betas.tumor=temp1, #Corrected tumor beta values
  betas.microenvironment=temp2, #Corrected microenvironment beta values
  cpg.populations=temp4, #Methylation patterns (populations) of eachCpG
  reg.slopes=temp5, #Slopes of the populations
  reg.intercepts=temp6, #Intercepts of the populations
  reg.mean_abs_residuals=temp7, #Mean of the residual's absolute value per population
  reg.RSE=temp7, #Residual standard error
  reg.SSx=temp8, #Sum of squares of the y axis
  reg.SSy=temp9, #Sum of squares of the x axis
  reg.df=temp10 #Degrees of freedom
)

# =====================
# CREATING OUTPUT FILES
# =====================



##Creating a prefix for the output files
out_prefix <- "output_training"

#Defining a function to export the dataframes as csv files
df_to_csv <- function(df, filename) {
  write.csv(df, filename)
}

#Creating output files per each dataframe of the output_list list
lapply(names(output_list), function(n) {
  df_to_csv(output_list[[n]],filename=paste("./",out_prefix,"_",n,".csv",sep=""))
})

#cat("\n=================\n")
#cat ("PROCESS FINISHED")
#cat("\n=================\n")

##END