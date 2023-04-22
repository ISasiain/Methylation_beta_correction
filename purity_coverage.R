#!/usr/bin/Rscript

purity_value_per_sample <- function(pred_purity_confidence,interval_threshold=4 ,min_endpoint=0, max_endpoint=1) {
  
  #Creating a list to store the output
  output_list = list()

  #Removing NA values from the predicted purity dataframe and formatting the values
  #to three decimal positions
  pred_purity_confidence <- format(round(na.omit(pred_purity_confidence),3),nsmall=3)
  

  # =============================================
  # DETERMINE THE COVERAGE OF THE PURITY INTERVAL
  # =============================================

  #Creating a vector to store the sections with the maximum coverage
  coverage_per_section <- c()

  # Divide the (min_endpoint,max_endpoint) interval in sections of length 0.001.
  #Generating a list with the each section of the interval as keys (the keys 
  #are written as characters)
  sections <- seq(min_endpoint,max_endpoint,by=0.001)
  coverage_per_section <- setNames(rep(0,length(sections)), as.character(sections))
  

    #Check if the sections of the 0-1 purity interval are covered in each interval and determine the coverage (how 
    #many times is each section includedin the intervals) per section

  for (int in 1:nrow(pred_purity_confidence)) {

    pos <- as.character(seq(pred_purity_confidence[int,1], pred_purity_confidence[int,2], by=0.001))
    included_section <- names(coverage_per_section) %in% pos
    coverage_per_section[included_section] <- coverage_per_section[included_section] + 1

  }

  #par(mfrow=c(1,2))
  
  #reg <- lm(unname(coverage_per_section)~as.numeric(names(coverage_per_section)))
  #plot(y=coverage_per_section, x=as.numeric(names(coverage_per_section)), type="l", pwd=10, xlab="1 - Purity", ylab = "Coverage")
  #abline(reg, col="lightgreen", lwd=2)

  #Correcting the overrepresentation of purity values between 0.8 and 1. Fitting linear regression and using the resiuduals
  coverage_per_section <- setNames(residuals(lm(unname(coverage_per_section)~as.numeric(names(coverage_per_section)))),names(coverage_per_section))
  
  #plot(y=coverage_per_section, x=as.numeric(names(coverage_per_section)), type="l", pwd=10, xlab="1 - Purity", ylab = "Adapted coverage")
  #abline(h=0, col="lightgreen", lwd=2)
  
  # ======================================================
  # DETERMINE THE MAXIMUM COVERAGE ESTIMATES AND INTERVALS
  # ======================================================

  ## GETTING THE ESTIMATES

  #Appending the max value(s), the 1-Purity estimate(s), to the output_list 
  output_list[["1-Pur_estimates"]] <- sections[coverage_per_section >= max(coverage_per_section)]

  ## GETTING THE INTERVALS

  # Get the intervals with the maximum coverage. The minimum coverage threshold is the maximum value minus the interval_
  #threshold percentage selected (default value is 10%)
  selected_values <- sections[coverage_per_section >= max(coverage_per_section)-round(max(coverage_per_section)*interval_threshold/100)]

  #Defining variables to store parameters of the intervals
  start_val <- selected_values[1] # Start point of the interval
  ref_val <- NULL # Value to use as the refernce for the ref value of the loop
  end_val <- NULL # End point of the interval

  interval_list <- list() # A list to store the intervals detected

  # Iterate through the values over the coverage threshold (except the first element, as 
  # it has already been assigned to start_val)
  for (val in selected_values[-1]) {
    
    #If ref_val is not defined assign it to the start_value
    if (is.null(ref_val)) {
      ref_val <- start_val
    } 

    # Compare the current value with ref_val to check if they are contiguous (the differnece between them is 0.001)
    if (format(val - ref_val,2)==0.001) {
      #If they were contiguous the current value would be assigned to ref_val
      ref_val <- val
    } else {
      #If they are not contiguous the reference value will be selected as the endpoint of the interval
      end_val <- ref_val
      #The interval is added to the interval list
      interval_list <- c(interval_list, list(c(start_val, end_val)))
        

      #Reestarting the loop after identifying an interval
      start_val <- val # Assign the current value to start_val
      end_val <- NULL # Assign NULL to end_val
      ref_val <- NULL # Assign NULL to ref_val
      
    } 
  }
  
  #Add the last selected interval to the interval list at the end of the for loop
  if (!is.null(ref_val)) { 
    end_val <- ref_val #Assign the ref_val to end_val if ref_val is defined
  } else {end_val <- start_val}

  #Add the interval to the interval list
  interval_list <- c(interval_list, list(c(start_val, end_val)))
  

  # ==============================================
  # REMOVING THE INTERVALS WITHOUT A MAXIMUM VALUE
  # ==============================================
  
# Checking if any of the maxs identified are inside the detected intervals. If they are the interval will be removed
interval_list <- Filter(function(interval) {
    length(intersect(seq(interval[1], interval[2], by = 0.001), format(output_list[["1-Pur_estimates"]], 3))) != 0
}, interval_list)
  
  #Adding the intervals detected to the output list
  output_list[["interval(s)"]] <- interval_list

  #The maximum coverage interval list will be returned
  return(output_list)
}