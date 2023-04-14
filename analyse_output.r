#Setting the working directory
setwd("~/Methylation/adjustBetas")

#out_ls1 <- readRDS("purity_est_out.RData")
out_ls1 <- readRDS("out_with_reg.RData")

if(!requireNamespace("caTools", quietly = TRUE)) {
  install.packages("caTools") }

library(caTools)


##data set reduced to top 5000 CpGs in set of 630 TCGA BRCA tumors
load("workspace_tcgaBrca_top5000.RData")

#Setting the seed for spliting the data
set.seed(1)

#Splitting the data
train_samples <- sample.split(colnames(betaUnadjusted), SplitRatio=0.8)

#Creating training and validation datasets for the purity_vector
purity_training <- purityVector[train_samples]
purity_validation <- purityVector[!train_samples]

vec_of_vec_nums <- c()

#How many intervals per sample
for (sample in names(purity_validation)) {

    vec_of_vec_nums <- c(vec_of_vec_nums, length((out_ls1[[sample]][['interval(s)']])))

}

#How many intervals per sample
for (sample in names(purity_validation[vec_of_vec_nums==1])) {

    cat("\n____________________________________________\n")
    print(sample)
    cat("\nInterval: \n")
    print(out_ls1[[sample]][["interval(s)"]])
    cat("\nEstimates (maxs) \n")
    print(out_ls1[[sample]][["1-Pur_estimates"]])
    cat("\n\nActual 1-purity", 1-purity_validation[sample], "")
    cat("\n____________________________________________\n")

}


for (i in 1:max(vec_of_vec_nums)){
cat(i,"-> ",length(vec_of_vec_nums[vec_of_vec_nums==i]), "\n")
}

included_in_interval <- 0
close_to_interval <- 0

for (sample in names(purity_validation[vec_of_vec_nums==1])) {
  
  print(out_ls1[[sample]][["interval(s)"]][[1]][1])
  print(out_ls1[[sample]][["interval(s)"]][[1]][2])
  print(1-purity_validation[sample])
  if ((1-purity_validation[sample]) >= out_ls1[[sample]][["interval(s)"]][[1]][1] & (1-purity_validation[sample]) <= out_ls1[[sample]][["interval(s)"]][[1]][2] ) {

    included_in_interval <- included_in_interval + 1
    
  } 
}

cat("\nIncluded_in_interval:", included_in_interval)
