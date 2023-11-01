#!/usr/bin/Rscript

# Script Name: analyze_output.r

# Description:
# This script performs an analysis of sample purity predictions from different datasets.
# It takes R objects containing predicted purity values for various folds and cancer types as input and generates boxplots
# along with quality metrics to assess the accuracy of these predictions.

# Required R Packages:
# - ggplot2: Used for creating data visualizations.
# - stringr: Employed for string manipulation.
# - tibble: Useful for data manipulation.
# - tidyr: Essential for data tidying.
# - dplyr: Needed for data manipulation.
# - optparse: Utilized for parsing command-line arguments.

# Loading Required Packages:
# The script first ensures a specific CRAN mirror is set, and it installs the required packages if not already installed.
# Then, it loads the 'ggplot2', 'stringr', 'tibble', 'tidyr', 'dplyr', and 'optparse' packages for data analysis and argument parsing.

# Configuration of Command Line Arguments:
# This script uses the 'optparse' package to configure and parse command-line arguments.
# The available options include:
# - "--path_to_estimated_1-P": Path to an R object containing estimated 1-Purity values.
# - "--path_to_cpgs_used_per_sample": Path to an R object with cpgs used for estimating purity per sample.
# - "--path_to_original_betas": Path to an R object with original betas used for estimating purity per sample.
# - "--path_to_actual_1-P": Path to an R object containing actual 1-Purity values.
# - "--output_prefix": Prefix for the output plots (default: "result_analysis").
# - "--plot_dis_vs_ploidy": Set to TRUE to plot the Distance to estimate VS ploidy (default: FALSE).
# - "--path_to_ploidy": Path to an R object containing actual ploidy values if "--plot_dis_vs_ploidy" is set to TRUE.

# Parsing Command Line Arguments:
# The script uses the 'optparse' package to parse the command-line arguments, including paths to various R objects and configuration options.

# Processing Command Line Arguments:
# The provided arguments are processed, and relevant data is loaded, including estimated purities, actual purities, cpgs used, and original betas.

# Quality Metric Calculation:
# The script calculates quality metrics, such as interval width, distance to interval, distance to estimate, and distance multiplied by width, for each prediction.

# Quality Metrics Determination:
# Quality metrics for each prediction and fold are determined and stored in a data frame.

# Plotting Results:
# The script generates various data visualizations to assess the quality of predictions. These include:
# - Histogram of Distance to Interval
# - Histogram of Distance to Estimate
# - Histogram of Interval Width
# - Scatterplot of Actual Purity vs. Estimated Purity
# - Scatterplot of Distance to Estimate vs. Actual Purity
# - Scatterplot of Distance to Estimate vs. Estimated Purity
# - Scatterplot of Distance to Estimate vs. Ploidy (if enabled)
# - Density Plot of Error Percentage (Distance to Estimate)
# - Cumulative Density Plot of Error Percentage (Distance to Estimate)
# - Histogram of the Percentage of Times Each CpG is Used to Estimate Purity
# - Histogram of the Number of CpGs Used to Estimate Purity per Sample
# - Heatmap of CpGs Used in Each Sample
# - Scatterplot of the Number of CpGs vs. Actual Purity of Samples

# Example Usage:
# Users can run the script from the command line with various options to analyze and visualize sample purity predictions.

# Author: Iñaki Sasiain Casado
# Affiliation: Johan Staaf lab, LU



# ==========================================
#      IMPORTING THE REQUIRED PACKAGES 
# ==========================================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

if(!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2") }

library(ggplot2)

if(!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2") }

library(reshape2)

if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

library(optparse)

# =============================================================
#      PARSING COMMAND LINE ARGUMENTS AND LOADING THE DATA
# =============================================================

argument_list <- list(

  make_option(c("-e", "--path_to_estimated_1-P"), type="character",  
              help="The user must especify the path of the R object containing the estimated 1-Pur values",
              metavar = "[path]"),

  make_option(c("-c", "--path_to_cpgs_used_per_sample"), type="character",  
              help="The user must especify the path of the R object containing the cpgs used to estimate purity per sample",
              metavar = "[path]"),

  make_option(c("-b", "--path_to_original_betas"), type="character",
              help="The user must especify the path of the R object containing the original betas used to estimate purity per sample",
              metavar="[path]"),

  make_option(c("-a", "--path_to_actual_1-P"), type="character",
              help="The user must especify the path of the R object containing the actual 1-Pur values",
              metavar = "[path]"),

  make_option(c("-o", "--output_prefix"), type="character", default="result_analysis",
              help="The user can especify a prefox for the output plots. Default [%default]",
              metavar = "[prefix]"),

  make_option(c("-P", "--plot_dis_vs_ploidy"), type="logical", default=FALSE,
              help="This argument must be set to TRUE in order to plot the Distance to estimate VS ploidy. Default [%default]",
              metavar = "[TRUE/FALSE]"),

  make_option(c("-p", "--path_to_ploidy"), type="character",
              help="The user must especify the path of the R object containing the actual ploidy values if the -P flag has been set to TRUE.",
              metavar = "[path]")

)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program produces plots to analyse the performance of the purity prediction algorithm."))


#Loading the estimated and actual purities and the cpgs used
out_ls <- readRDS(arguments$"path_to_estimated_1-P")
purity_validation <- readRDS(arguments$"path_to_actual_1-P")
cpgs_ls <- readRDS(arguments$"path_to_cpgs_used_per_sample")
original_betas <- readRDS(arguments$"path_to_original_betas")

print(purity_validation)

#Adding this to adapt the sample names
#names(purity_validation) <- paste(names(purity_validation), "-01A", sep="")


# ===========================================================
#      CREATING FUNCTION TO DETERMINE QUALITY PARAMETERS
# ===========================================================


det_qual <- function(output_list, purity_vector) {
  
  columns <- c("Interval's width", "Dis_to_int", "Dis_to_est", "Distance*Width")
  results <- data.frame(matrix(nrow = 0, ncol = length(columns)))

  colnames(results) <- columns
  
  quality_values <- lapply(names(output_list), function(sample) {

    # Determining the interval's width
    width <- output_list[[sample]][["interval(s)"]][[1]][2] - output_list[[sample]][["interval(s)"]][[1]][1]
    
    # Determining the distance to the interval
    if ((length(intersect(seq(output_list[[sample]][["interval(s)"]][[1]][1], output_list[[sample]][["interval(s)"]][[1]][2], by=0.001),purity_vector[sample])) == 1)) {
      dis_to_interval <- 0
    } else {
      dis_to_interval <- round(min(abs(c(output_list[[sample]][["interval(s)"]][[1]][1] - (1-purity_validation[sample]), output_list[[sample]][["interval(s)"]][[1]][2] - (1-purity_validation[sample])))),3)
    }

    # Determining the distance to the estimate (the mean is calculated in case there were more than one maximums. THis must be removed after softening the coverage plot)
    dis_to_estimate <- mean(round(abs(1-purity_vector[sample] - output_list[[sample]][["1-Pur_estimates"]]),3))
    
    # Adding the calculated values to a dataframe
    data.frame("Interval's width" = width, 
               "Dis_to_int" = dis_to_interval,
               "Dis_to_est" = dis_to_estimate,
               "Distance*Width" = round(dis_to_interval * width,3))
  })
  
  # Append all the samples' to the final data frame
  for (i in seq_along(quality_values)) {
    results <- rbind(results, quality_values[[i]])
  }
  
  row.names(results) <- names(output_list)
  return(results)
}

out_df <- det_qual(out_ls, purity_validation)

# ============================================================
#     PROCESS DATA TO ANALYSE CPG USAGE IN THE PREDICTION
# ============================================================

# Storing the ids of all the samples and cpgs
all_the_cpgs <- row.names(original_betas)
all_the_samples <- names(cpgs_ls)

# Creating a boolean matrix with cpgs as columns and samples as rows
matrix_to_heatmap <- matrix(FALSE, ncol = length(all_the_cpgs), nrow = length(all_the_samples))
colnames(matrix_to_heatmap) <- all_the_cpgs
rownames(matrix_to_heatmap) <- all_the_samples

# Setting the positions of the matrix to TRUE if the cpg is included for each sample
for (sample in all_the_samples) {

  included_cpgs <- all_the_cpgs %in% cpgs_ls[[sample]]
  matrix_to_heatmap[sample, included_cpgs] <- TRUE

}

# Determine how many times is each cpgs used to predict the purity of each sample
cpg_counts_df <- data.frame(
  TRUES = colSums(matrix_to_heatmap),
  FALSES = colSums(!matrix_to_heatmap),
  per_TRUES = colSums(matrix_to_heatmap)/length(all_the_samples),
  per_FALSES = colSums(!matrix_to_heatmap)/length(all_the_samples)
)

# Sorting the cpg dataframe based on the number of samples in which each cpg is included
cpg_counts_df <- cpg_counts_df[order(cpg_counts_df$TRUES, decreasing=TRUE),]

# Creating and sorting a vector with the number of cpgs used to predict purity per sample
num_cpgs <- sapply(cpgs_ls, length)
num_cpgs <- num_cpgs[order(num_cpgs)]

#Converting te vector into a dataframe to be plotted using ggplot
cpgs_per_sample_df <- data.frame()
cpgs_per_sample_df <- cbind(num_cpgs)
cpgs_per_sample_df <- cbind(cpgs_per_sample_df, purity_validation[order(num_cpgs)])

colnames(cpgs_per_sample_df) <- c("num_of_cpgs", "actual_1_minus_P")



#Sorting matrix_to_heatmap based on the num_cpgs and cpg_counts data
matrix_to_heatmap <- matrix_to_heatmap[names(num_cpgs),rownames(cpg_counts_df)]

# Transforming the format of the matrix to be plotted using ggplot
df_to_heatmap <- melt(matrix_to_heatmap)

# =======================================
#      CALCULATING SUMMARY STATISTICS
# =======================================

print(summary(out_df)) 

# ==============================
#      PLOTTING THE RESULTS
# ==============================

## HISTOGRAM OF DISTANCE TO INTERVAL
ggplot(out_df, aes(x=Dis_to_int)) +
  geom_histogram(binwidth=0.01, color="black", fill="lightblue") +
  xlim(0,0.8) +
  xlab("Distance to interval") +
  ylab("Frequency") +
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())
ggsave(paste(arguments$output_prefix, "dis_to_int.barplot.png", sep="."))

## HISTOGRAM OF DISTANCE TO ESTIMATE
ggplot(out_df, aes(x=Dis_to_est)) +
  geom_histogram(binwidth=0.01, color="black", fill="lightblue") +
  xlim(0,0.8) +
  xlab("Distance to estimate") +
  ylab("Frequency") +
  theme_classic() +
  theme(plot.title = element_text(size = 22),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())
ggsave(paste(arguments$output_prefix, "dis_to_int.barplot.png", sep="."))

## HISTOGRAM OF WIDTH
ggplot(out_df, aes(x=Interval.s.width)) +
  geom_histogram(binwidth=0.0025, color="black", fill="lightgreen") +
  xlim(0,0.16) +
  xlab("Interval's width") +
  ylab("Frequency") +
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())
ggsave(paste(arguments$output_prefix, "int_width.barplot.png",sep="."))


## SCATTERPLOT WITH INTERVAL. Actual purity

#Estimated values
upper_lim <- c()
estimate <- c()
lower_lim <- c()

#Actual 1-P value
actual_1_minus_P <- 1-purity_validation[names(out_ls)]
index <- sort(actual_1_minus_P, index.return=TRUE)$ix

if (arguments$plot_dis_vs_ploidy) {

ploidy <- readRDS(arguments$path_to_ploidy)


#Actual NC 
actual_NC <- ((2 * (1-purity_validation[names(out_ls)])))/(purity_validation[names(out_ls)] * ploidy[names(out_ls)] + 2 * (1-purity_validation[names(out_ls)]))

#Getting the values of all the samples in vectors
for (sample in names(out_ls)) {
  
  lower_lim <- c(lower_lim, out_ls[[sample]][["interval(s)"]][[1]][1])
  estimate <- c(estimate ,mean(out_ls[[sample]][["1-Pur_estimates"]]))
  upper_lim <- c(upper_lim, out_ls[[sample]][["interval(s)"]][[1]][2])
}

#Creating a dataframe with the sorted vectors
plot_df <- data.frame(
  actual=actual_NC[order(index)],
  upper=upper_lim[order(index)],
  est= estimate[order(index)], 
  lower=lower_lim[order(index)]
)

print(plot_df)

ggplot(data=plot_df, aes(x=actual, y=est)) +
  xlim(0,1) + ylim(0,1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=lower,ymax=upper), col="#336633", size=0.8) +
  geom_abline(slope = 1, intercept = 0, size=1.5, col="grey") + 
  xlab("Actual proportion of non-cancer DNA") + 
  ylab("Estimated proportion of non-cancer DNA") +
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())
ggsave(paste(arguments$output_prefix, "NC_Act_vs_Est.scatterplot.png",sep="."))

}


#Getting the values of all the samples in vectors
for (sample in names(out_ls)) {
  
  lower_lim <- c(lower_lim, out_ls[[sample]][["interval(s)"]][[1]][1])
  estimate <- c(estimate ,mean(out_ls[[sample]][["1-Pur_estimates"]]))
  upper_lim <- c(upper_lim, out_ls[[sample]][["interval(s)"]][[1]][2])
}

#Creating a dataframe with the sorted vectors
plot_df <- data.frame(
  actual=actual_1_minus_P[order(index)],
  upper=upper_lim[order(index)],
  est=estimate[order(index)],
  lower=lower_lim[order(index)]
)


## SCATTERPLOT WITH INTERVAL. 

ggplot(data=plot_df, aes(x=actual, y=est)) +
  xlim(0,1) + ylim(0,1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=lower,ymax=upper), col="#336633", size=0.8) +
  geom_abline(slope = 1, intercept = 0, size=1.5, col="grey") + 
  xlab("Actual 1-Purity") + 
  ylab("Estimated 1-Purity") +
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())
ggsave(paste(arguments$output_prefix, "Act_vs_Est.scatterplot.png",sep="."))




## EVOLUTION OF DISTANCE TO ESTIMATE IN FUNCTION OF ACTUAL AND ESTIMATED PURITY

# Plot of Distance to interval vs Actual 1-Purity
ggplot(out_df[order(index),], aes(x=actual_1_minus_P[order(index)], y=Dis_to_est)) +
  geom_point(color="blue") +
  xlab("Actual 1-Purity") +
  ylab("Dis_to_estimate") + 
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(paste(arguments$output_prefix, "Dis_vs_ac1-P.scatterplot.png", sep="."))

# Plot of Distance to interval vs Estimated 1-Purity
ggplot(out_df, aes(x=plot_df$est, y=out_df$Dis_to_est[order(index)])) +
  geom_point(color="red") +
  xlab("Estimated 1-Purity") +
  ylab("Distance to estimate") + 
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())
ggsave(paste(arguments$output_prefix, "Dis_vs_es1-P.scatterplot.png",sep="."))


## DISTANCE TO INTERVAL VS PLOIDY
if (arguments$plot_dis_vs_ploidy) {

    ploidy <- readRDS(arguments$path_to_ploidy)

    ggplot(out_df, aes(x=ploidy[row.names(out_df)], y=out_df$Dis_to_est[order(index)])) +
    geom_point(color="red") +
    xlab("Ploidy") +
    ylab("Disance to estimate") + 
    theme_classic() +
    theme(axis.title = element_text(size = 22),
          axis.text = element_text(size = 20),
          panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
          panel.grid.minor = element_blank())
  ggsave(paste(arguments$output_prefix, "Dis_vs_ploidy.scatterplot.png",sep="."))

}

## ERROR PERCENTAGE (DISTANCE TO ESTIMATE) DISTRIBUTION

ggplot(out_df, aes(x=Dis_to_est)) +
  geom_density(alpha=0.5, fill="red") +
  xlab("Deviation") +
  xlim(0,0.2) +
  ylab("% of the samples") + 
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(paste(arguments$output_prefix, "error.densityplot.png",sep="."))

## CUMULAIVE ERROR PERCENTAGE (DISTANCE TO ESTIMATE) DISTRIBUTION

ggplot(out_df, aes(x=Dis_to_est)) +
  stat_ecdf(color="red", size=2) +
  xlab("Deviation") +
  xlim(0,0.2) +
  ylab("% of the samples") + 
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(paste(arguments$output_prefix, "error.cum_densityplot.png",sep="."))

# PLOT THE PERCENTAGE OF TIMES IN WHICH EACH CPG IS USED TO ESTIMATE PURITY

ggplot(cpg_counts_df, aes(x=per_TRUES)) +
  geom_histogram(binwidth=0.0025, color="black", fill="red") +
  xlab("% of samples in which each CpGs are included") +
  ylab("CpGs") +
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(paste(arguments$output_prefix, "prop_cpgs_included.png",sep="."))


## HISTOGRAM OF THE CPGS USED TO ESTIMATE SAMPLE PURITY

ggplot(data.frame(cpgs_per_sample_df), aes(x=num_cpgs)) +
  geom_histogram(color="black", fill="red") +
  xlab("Number of CpGs used to estimate purity") +
  ylab("Number of samples") +
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(paste(arguments$output_prefix, "num_cpg_per_sample.png",sep="."))

## PLOTTING HEATMAP OF THE CPGS USED IN EACH SAMPLE

# Define custom color scale
custom_colors <- c("green", "red")
custom_breaks <- c(TRUE, FALSE)
custom_labels <- c("Included", "Not included")

# Plotting the heat map
ggplot(df_to_heatmap, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = custom_colors,
                    breaks = custom_breaks,
                    labels = custom_labels) +
  labs(x = "CpGs", y = "Samples", title = "CpGs included per sample")

ggsave(paste(arguments$output_prefix, "heatmap_cpg_per_sample.png",sep="."), width = 10, height = 30, limitsize=FALSE)


##  CpG USAGE VS ACTUAL PURITY OF SAMPLES

ggplot(data=data.frame(cpgs_per_sample_df), aes(x= actual_1_minus_P, y= num_of_cpgs )) +
  geom_point() +
  xlab("Actual 1-P") +
  ylab("Number of CpGs used to estimate purity") +
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(paste(arguments$output_prefix, "num_cpgs_actual_1-P.scatterplot.png",sep="."))