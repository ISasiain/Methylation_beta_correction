#Setting the working directory
setwd("~/Methylation/adjustBetas")

vec_of_vec_nums <- c()
included_in_interval <- 0

a <- load("out_with_reg.RData")

#How many intervals per sample
for (sample in names(purity_validation)) {

    vec_of_vec_nums <- c(vec_of_vec_nums, length((out_ls[[sample]][['interval(s)']])))

}

#How many intervals per sample
for (sample in names(purity_validation[vec_of_vec_nums==1])) {

    cat("\n____________________________________________\n")
    print(sample)
    cat("\nInterval: \n")
    print(out_ls[[sample]][["interval(s)"]])
    cat("\nEstimates (maxs) \n")
    print(out_ls[[sample]][["1-Pur_estimates"]])
    cat("\n\nActual 1-purity", 1-purity_validation[sample], "")
    cat("\n____________________________________________\n")

}

print(included_in_interval)

for (i in 1:max(vec_of_vec_nums)){
cat(i,"-> ",length(vec_of_vec_nums[vec_of_vec_nums==i]), "\n")
}

