## Function to predict purity values based on the individual betas

predicting_purity <- function(beta,slopes,intercepts,RSE,SS_x,SS_y,degrees_of_freedom,slope_threshold=0.4, RSE_threshold=0.2, alpha=0.85) {

    #Identifying the regressions whose residual standard error is over the threshold
    over_residual_th <- RSE > RSE_threshold
    
    # Identifying the regressions whose slope is in the (0-slope_threshold, 0+slope_threshold)
    # interval to NA. This regressions will be ignored in the estimation of each sample's 
    # purity value as they are not informative.
    uninformative <- slopes > (0-slope_threshold) & slopes < (0+slope_threshold) 

    # Setting the parameters of the uninformative and overt the mean residual 
    # threshold regressions to NA
    slopes[over_residual_th | uninformative] <- NA
    intercepts[over_residual_th | uninformative] <- NA
    RSE[over_residual_th | uninformative] <- NA
    SS_x[over_residual_th | uninformative] <- NA
    SS_y[over_residual_th | uninformative] <- NA
    degrees_of_freedom[over_residual_th | uninformative] <- NA

    # Removing NA values from all the previously defined vectors. This values 
    # would be on the vectors if there were less than three methylation patterns 
    # identified in the CpG analysed or if there were uninformative regressions.
    slopes <- slopes[!is.na(slopes)]
    intercepts <- intercepts[!is.na(intercepts)]
    RSE <- RSE[!is.na(RSE)]
    SS_x <- SS_x[!is.na(SS_x)]
    SS_y <- SS_y[!is.na(SS_y)]
    degrees_of_freedom <- degrees_of_freedom[!is.na(degrees_of_freedom)]


    # Check that all the regressions are not uninformative and execute the following code 
    # to estimate the 1-Purity value. If they weren't NA would be assigned to the 1-Purity value,
    # so that CpG would not be used to estimate the final 1-Purity of the sample

    if (length(slopes[!is.na(slopes)])!=0) {

        #Creating a vector to store the regressions in which the beta value analysed is included
        matches <- rep(NA, length(slopes))
        
        #Iterate through the regression lines predicted per methylation pattern identified
        for (reg in 1:length(slopes)) {

            #Check if the beta value analysed is in the regression line
            if ((beta - intercepts[reg]) / slopes[reg] >= 0 && (beta - intercepts[reg]) / slopes[reg] <= 1) {
                #Is the beta value is in the regression matches[reg] will be set to TRUE
                matches[reg] <- TRUE
            #Is the beta value is not in the regression matches[reg] will be set to FALSE
            } else {
                matches[reg] <- FALSE
            }
        }
        
        # If the length of the matches vector when matches==TRUE is one a unique CpG methylation pattern (regression)
        # can be assigned to each beta value. If this value is different to one no signoficant association to a pattern can 
        # be achieved
        if (length(matches[matches])==1) {

            #Store the identified population into a variable
            identified_pop <- match(TRUE, matches)

            #Calculate the parameters of the inverse regression. 1-Purity VS beta instead of beta VS 1-Purity
            inv_slope <- 1/slopes[identified_pop]
            inv_intercept <- -intercepts[identified_pop]/slopes[identified_pop]
            inv_df <- degrees_of_freedom[identified_pop]
            inv_RSE <- RSE[identified_pop] * sqrt(SS_x[identified_pop]/SS_x[identified_pop])

            #Predicting the 1-purity value from the inverse regression
            predicted_1_minus_p <- inv_slope*beta + inv_intercept

            #Getting the prediction interval
            one_minus_purity <- c(
                predicted_1_minus_p - qt((1-alpha)/2, inv_df) * inv_RSE,
                predicted_1_minus_p + qt((1-alpha)/2, inv_df) * inv_RSE
            )

        } else {
            #If the CpG can not be assigned to a single population (methylation pattern) the 1-purity value will not be assigned.
            one_minus_purity <- c(NA,NA)
            }
        
        #If all the regressions are uninformative NA value will be assigned to one_minus_purity.
        } else {
            one_minus_purity <- c(NA,NA)
        }

    if (! is.na(one_minus_purity[1] && is.na(one_minus_purity[2]))) {

        #Setting the interval's endpoints to 0 or 1 if one of the limits lower or higher than those values.
        if (one_minus_purity[1]<0) {one_minus_purity[1]=0}
        if (one_minus_purity[1]>1) {one_minus_purity[1]=1}

        if (one_minus_purity[2]<0) {one_minus_purity[2]=0}
        if (one_minus_purity[2]>1) {one_minus_purity[2]=1} 

    }   

    #Returning the 1-purity value and the identified population
    return(one_minus_purity)
}

setwd("~/Methylation/adjustBetas")

interval_df = data.frame()

sample <- colnames(read.csv("output_training_betas.original.csv", sep=",", row.names=1))[1]

for (cpg in rownames(read.csv("output_training_betas.original.csv", sep=",", row.names=1))) {

my_beta <- read.csv("output_training_betas.original.csv", sep=",", row.names=1)[cpg,sample]
my_slopes <- read.csv("output_training_reg.slopes.csv", sep=",", row.names=1)[cpg,]
my_intercepts <- read.csv("output_training_reg.intercepts.csv", sep=",", row.names=1)[cpg,]
my_RSE <- read.csv("output_training_reg.RSE.csv", sep=",", row.names=1)[cpg,]
my_SSx <- read.csv("output_training_reg.SSx.csv", sep=",", row.names=1)[cpg,]
my_SSy <- read.csv("output_training_reg.SSy.csv", sep=",", row.names=1)[cpg,]
my_df <- read.csv("output_training_reg.df.csv", sep=",", row.names=1)[cpg,]

interval_df <- rbind(interval_df, predicting_purity(my_beta, my_slopes, my_intercepts, my_RSE, my_SSx, my_SSy, my_df))

}
