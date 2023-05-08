#!/usr/bin/env Rscript

# ==========================================
#      IMPORTING THE REQUIRED PACKAGES 
# ==========================================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

if(!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2") }

library(ggplot2)

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

  make_option(c("-a", "--path_to_actual_1-P"), type="character",
              help="The user must especify the path of the R object containing the actual 1-Pur values",
              metavar = "[path]"),

  make_option(c("-o", "--output_prefix"), type="character", default="result_analysis",
              help="The user can especify a prefox for the output plots. Default [%default]",
              metavar = "[prefix]")

)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program produces plots to analyse the performance of the purity prediction algorithm."))


#Loading the estimated and actual purities
out_ls <- readRDS(arguments$"path_to_estimated_1-P")
purity_validation <- readRDS(arguments$"path_to_actual_1-P")

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
  
  # Append the results for all samples to the final data frame
  for (i in seq_along(quality_values)) {
    results <- rbind(results, quality_values[[i]])
  }
  
  row.names(results) <- names(output_list)
  return(results)
}

out_df <- det_qual(out_ls, purity_validation)

# =======================================
#      CALCULATING SUMMARY STATISTICS
# =======================================

print(summary(out_df)) 

# ==============================
#      PLOTTING THE RESULTS
# ==============================

## HISTOGRAM OF DISTANCE
ggplot(out_df, aes(x=Dis_to_int)) +
  geom_histogram(binwidth=0.01, color="black", fill="lightblue") +
  xlim(0,0.8) +
  ggtitle("Not corrected") +
  xlab("Distance to interval") +
  ylab("Frequency") +
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())
ggsave(paste(arguments$output_prefix, "dis_to_int.barplot.png", sep="."))

## HISTOGRAM OF WIDTH
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
ggsave(paste(arguments$output_prefix, "int_width.barplot.png",sep="."))


## SCATTERPLOT WITH INTERVAL

#Estimated values
upper_lim <- c()
estimate <- c()
lower_lim <- c()

#Actual 1-P value
actual_1_minus_P <- 1-purity_validation
index <- sort(actual_1_minus_P, index.return=TRUE)$ix

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

ggplot(data=plot_df, aes(x=actual, y=est)) +
  xlim(0,1) + ylim(0,1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=lower,ymax=upper), col="#336633", size=0.8) +
  geom_abline(slope = 1, intercept = 0, size=1.5, col="grey") + 
  ggtitle("Prediction of 1-Purity") +
  xlab("Actual 1-Purity") + 
  ylab("Estimated 1-Purity") +
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())
ggsave(paste(arguments$output_prefix, "Act_vs_Est.scatterplot.png",sep="."))

## EVOLUTION OF DISTANCE TO INTERVAL IN FUNCTION OF ACTUAL AND ESTIMATED PURITY

# Plot of Distance to interval vs Actual 1-Purity
ggplot(out_df[order(index),], aes(x=actual_1_minus_P[order(index)], y=Dis_to_int)) +
  geom_point(color="blue") +
  xlab("Actual 1-Purity") +
  ylab("Dis_to_interval") + 
  ggtitle("Distance to interval vs Actual 1-Purity") +
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(paste(arguments$output_prefix, "Dis_vs_ac1-P.scatterplot.png",sep="."))

# Plot of Distance to interval vs Estimated 1-Purity
ggplot(out_df, aes(x=plot_df$est, y=out_df$Dis_to_int[order(index)])) +
  geom_point(color="red") +
  ggtitle("Distance to interval vs Estimated 1-Purity") +
  xlab("Estimated 1-Purity") +
  ylab("Disance to interval") + 
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())
ggsave(paste(arguments$output_prefix, "Dis_vs_es1-P.scatterplot.png",sep="."))


## ERROR PERCENTAGE (DISTANCE TO ESTIMATE) DISTRIBUTION

ggplot(out_df, aes(x=out_df$Dis_to_est*100), fill="lighred") +
  geom_density(alpha=0.5) +
  ggtitle("Distribution of the deviation") +
  xlab("% Deviation") +
  xlim(0,20) +
  ylab("% of the samples") + 
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(paste(arguments$output_prefix, "error.densityplot.png",sep="."))

## CUMULAIVE ERROR PERCENTAGE (DISTANCE TO ESTIMATE) DISTRIBUTION

ggplot(out_df, aes(x=out_df$Dis_to_est*100)) +
  stat_ecdf(color="lightred", size=2) +
  ggtitle("Cumulative distribution of the deviation") +
  xlab("% Deviation") +
  xlim(0,20) +
  ylab("% of the samples") + 
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(paste(arguments$output_prefix, "error.cum_densityplot.png",sep="."))