#Opening slopes file
slopes_pops <- read.csv("output_training_reg.mean_abs_residuals.csv", sep=",",row.names = 1)

#Potting the data in a histogram
par(mfrow=c(1,1))
hist(na.omit(as.vector((as.matrix(slopes_pops)))),
     breaks=50,
     xlab="Slopes of regressions",
     main="Slopes of regressions (Beta VS 1-Purity)")

#This function creates a dataframe to remove the populations whose slopes 
#are around 0, and therefore, are not informative.

slope_filter <- function(min_th, max_th, slopes_df) {

  slopes_df <= min_th | slopes_df >= max_th

}

hist(na.omit(as.vector((as.matrix(slopes_pops))))[na.omit(as.vector((as.matrix(slope_filter(-0.5,0.5,slopes_pops)))))],
     breaks=50,
     xlab="Slopes of regressions",
     main="Slopes of regressions (Beta VS 1-Purity)")