  one_minus_purity_vec <- c()  

  for (sample in colnames(unadj_validation)) {
  
  ##set/create own home directory below:
  setwd("~/Methylation/adjustBetas")
  
  #Sourcing the predicting_purity function
  source("./predicting_purity.R", encoding = "UTF-8")
  
  #Reading slopes and intervcepts values obtained from training data
  slopes <- read.csv("./output_training_reg.slopes.csv", sep=",", row.names=1)
  intercepts <- read.csv("./output_training_reg.intercepts.csv", sep=",", row.names=1)
  
  result_vector <- c()
  
  #Test the function for all the samples in the validation dataset
  for (cpg in 1:length(unadj_validation[,sample])) {
  
      result_vector <- c(result_vector, predicting_purity(unadj_validation[cpg,sample], slopes[rownames(unadj_validation)[cpg],], slopes[rownames(unadj_validation)[cpg],]))
  }
  
  one_minus_purity_vec <- c(one_minus_purity_vec, mean(result_vector,na.rm = TRUE))
  
}

plot(purity_validation, 1-one_minus_purity_vec, xlim=c(0,1), ylim=c(0,1))
