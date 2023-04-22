#Setting the working directory
setwd("~/Methylation/adjustBetas")

#Importing the required modules
library(ggplot2)

out_ls <- readRDS("out_without_reg_test.RData")
#out_ls <- readRDS("out_without_reg.RData")
out_ls_corr <- readRDS("out_with_reg_test.RData")

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


# for (i in 1:max(vec_of_vec_nums)){
# cat(i,"-> ",length(vec_of_vec_nums[vec_of_vec_nums==i]), "\n")
# }
# 
# included_in_interval <- 0
# close_to_interval <- 0
# 
# for (sample in names(purity_validation[vec_of_vec_nums==1])) {
#   
#   print(out_ls[[sample]][["interval(s)"]][[1]][1])
#   print(out_ls[[sample]][["interval(s)"]][[1]][2])
#   print(1-purity_validation[sample])
#   if ((1-purity_validation[sample]) >= out_ls[[sample]][["interval(s)"]][[1]][1] & (1-purity_validation[sample]) <= out_ls[[sample]][["interval(s)"]][[1]][2] ) {
# 
#     included_in_interval <- included_in_interval + 1
#     
#   } 
# }
# 
# cat("\nIncluded_in_interval:", included_in_interval)

#Mean-width modified distance score

MWMDS <- function(output_list, purity_vector) {
  
  columns <- c("Interval's width", "Distance", "Distance*Width")
  results <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(results) <- columns
  
  quality_values <- lapply(names(output_list), function(sample) {
    
    # Determining the interval's width
    width <- output_list[[sample]][["interval(s)"]][[1]][2] - output_list[[sample]][["interval(s)"]][[1]][1]
    
    # Determining the distance
    if ((length(intersect(seq(output_list[[sample]][["interval(s)"]][[1]][1], output_list[[sample]][["interval(s)"]][[1]][2], by=0.001),purity_vector[sample])) == 1)) {
      distance <- 0
    } else {
      distance <- round(min(abs(c(output_list[[sample]][["interval(s)"]][[1]][1] - (1-purity_vector[sample]), output_list[[sample]][["interval(s)"]][[1]][2] - (1-purity_vector[sample])))),3)
    }
    
    # Adding the calculated values to a dataframe
    data.frame("Interval's width" = width, 
               "Distance" = distance, 
               "Distance*Width" = round(distance * width,3))
  })
  
  # Append the results for all samples to the final data frame
  for (i in seq_along(quality_values)) {
    results <- rbind(results, quality_values[[i]])
  }
  
  row.names(results) <- names(output_list)
  return(results)
}

# ==============================
#      PLOTTING THE RESULTS
# ==============================

par(mfrow=c(3,2))

out_df <- MWMDS(out_ls, purity_validation)
out_df_corr <- MWMDS(out_ls_corr, purity_validation)

# Histogram of Distance
ggplot(out_df, aes(x=Distance)) +
  geom_histogram(binwidth=0.01, color="black", fill="lightblue") +
  xlim(0,0.8) +
  ggtitle("Not corrected") +
  xlab("Distances") +
  ylab("Frequency") +
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggplot(out_df_corr, aes(x=Distance)) +
  geom_histogram(binwidth=0.01, color="black", fill="lightblue") +
  xlim(0,0.8) +
  ggtitle("Corrected with reg") +
  xlab("Distances") +
  ylab("Frequency") +  
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

# Histogram of Interval's width
ggplot(out_df, aes(x=Interval.s.width)) +
  geom_histogram(binwidth=0.0025, color="black", fill="lightgreen") +
  xlim(0,0.16) +
  ggtitle("Not corrected") +
  xlab("Interval's width") +
  ylab("Frequency") +
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggplot(out_df_corr, aes(x=Interval.s.width)) +
  geom_histogram(binwidth=0.0025, color="black", fill="lightgreen") +
  xlim(0,0.16) +
  ggtitle("Corrected with reg") +
  xlab("Interval's width") +
  ylab("Frequency")  +
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

# # Histogram of Distance.Width
# ggplot(out_df, aes(x=Distance.Width)) +
#   geom_histogram(binwidth=0.0025, color="black", fill="pink") +
#   xlim(0,0.07) +
#   ggtitle("Not corrected") +
#   xlab("Distance * Width") +
#   ylab("Frequency")
# 
# ggplot(out_df_corr, aes(x=Distance.Width)) +
#   geom_histogram(binwidth=0.0025, color="black", fill="pink") +
#   xlim(0,0.07) +
#   ggtitle("Corrected with reg") +
#   xlab("Distance * Width") +
#   ylab("Frequency")

summary(out_df)
summary(out_df_corr)

#Scatterplot with errors

#Corrected values
upper_lim_corr <- c()
estimate_corr <- c()
lower_lim_corr <- c()

#UNcorrected values
upper_lim_uncorr <- c()
estimate_uncorr <- c()
lower_lim_uncorr <- c()

#Actual 1-P value
actual_1_minus_P <- 1-purity_validation
index <- sort(actual_1_minus_P, index.return=TRUE)$ix

par(mfrow=c(1,2))

for (sample in names(out_ls_corr)) {
  
  lower_lim_corr <- c(lower_lim_corr, out_ls_corr[[sample]][["interval(s)"]][[1]][1])
  estimate_corr <- c(estimate_corr ,mean(out_ls_corr[[sample]][["1-Pur_estimates"]]))
  upper_lim_corr <- c(upper_lim_corr, out_ls_corr[[sample]][["interval(s)"]][[1]][2])
}

corr_plot_df <- data.frame(
  actual=actual_1_minus_P[order(index)],
  upper=upper_lim_corr[order(index)],
  est=estimate_corr[order(index)],
  lower=lower_lim_corr[order(index)]
)

ggplot(data=corr_plot_df, aes(x=actual, y=est)) +
  xlim(0,1) + ylim(0,1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=lower,ymax=upper), col="#336633", size=0.8) +
  geom_abline(slope = 1, intercept = 0, size=1.5, col="grey") + 
  ggtitle("Corrected 1-Purity") +
  xlab("Actual 1-Purity") + 
  ylab("Estimated 1-Purity") +
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

for (sample in names(out_ls)) {
  
  upper_lim_uncorr <- c(upper_lim_uncorr, out_ls[[sample]][["interval(s)"]][[1]][1])
  estimate_uncorr <- c(estimate_uncorr ,mean(out_ls[[sample]][["1-Pur_estimates"]]))
  lower_lim_uncorr <- c(lower_lim_uncorr, out_ls[[sample]][["interval(s)"]][[1]][2])
}

uncorr_plot_df <- data.frame(
  actual=actual_1_minus_P[order(index)],
  upper=upper_lim_uncorr[order(index)],
  est=estimate_uncorr[order(index)],
  lower=lower_lim_uncorr[order(index)]
)

ggplot(data=uncorr_plot_df, aes(x=actual, y=est)) +
  xlim(0,1) + ylim(0,1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=lower,ymax=upper), col="#336633", size=0.8) +
  geom_abline(slope = 1, intercept = 0, size=1.5, col="grey") + 
  ggtitle("Uncorrected 1-Purity") +
  xlab("Actual 1-Purity") + 
  ylab("Estimated 1-Purity") +
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

#Plotting the evoliution of the distances in function of the actual and estimated purity

par(mfrow=c(1,2))

# Plot of Distance vs Actual 1-Purity (not corrected)
ggplot(out_df[order(index),], aes(x=actual_1_minus_P[order(index)], y=Distance)) +
  geom_point(color="blue") +
  xlab("Actual 1-Purity") +
  ylab("Distance") + 
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

# Plot of Distance vs Actual 1-Purity (corrected)
ggplot(out_df_corr[order(index),], aes(x=actual_1_minus_P[order(index)], y=Distance)) +
  geom_point(color="blue") +
  xlab("Actual 1-Purity") +
  ylab("Distance") + 
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

# Plot of Distance vs Estimated 1-Purity (not corrected)
ggplot(out_df, aes(x=uncorr_plot_df$est, y=out_df$Distance[order(index)])) +
  geom_point(color="red") +
  xlab("Estimated 1-Purity") +
  ylab("Distance") + 
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

# Plot of Distance vs Estimated 1-Purity (corrected)
ggplot(out_df_corr, aes(y=out_df_corr$Distance[order(index)], x=corr_plot_df$est)) +
  geom_point(color="red") +
  xlab("Estimated 1-Purity") +
  ylab("Distance") + 
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())


