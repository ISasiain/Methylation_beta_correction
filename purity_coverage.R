purity_value_per_sample <- function(pred_purity_confidence,interval_threshold=5 ,min_endpoint=0, max_endpoint=1) {
  
  #Creating a list to store the output
  output_list = list()

  #Removing NA values from the predicted purity dataframe and formatting the values
  #to three decimal positions
  pred_purity_confidence <- round(na.omit(pred_purity_confidence),3)
  

  # =============================================
  # DETERMINE THE COVERAGE OF THE PURITY INTERVAL
  # =============================================

  #Creating a vector to store the sections with the maximum coverage
  coverage_per_section <- c()

  # Divide the (min_endpoint,max_endpoint) interval in sections of length 0.001.
  #Generating a list with the each section of the interval as keys (the keys 
  #are written as characters)
  for (section in seq(0,max_endpoint-min_endpoint,by=0.001)) {
    coverage_per_section[as.character(section)] <- 0
  }
  
  #Iterate through the CpG's predicted beta value intervals per sample
  for (beta_est in 1:nrow(pred_purity_confidence)) {
    
    #Check if the sections of the 0-1 purity interval are covered in each interval and determine the coverage (how 
    #many times is each section includedin the intervals) per section
    for ( pos in as.character(seq(pred_purity_confidence[beta_est,1],pred_purity_confidence[beta_est,2],by=0.001))) {
      coverage_per_section[pos] = coverage_per_section[pos] + 1}
  }

  coverage_per_section <- setNames(residuals(lm(unname(coverage_per_section)~as.numeric(names(coverage_per_section)))),names(coverage_per_section))

  # =======================================================================
  # DETERMINE THE MAXIMUM COVERAGE INTERVALS WITHIN THE PURITY 0-1 INTERVAL
  # =======================================================================

  # Get the intervals with the maximum coverage. The minimum coverage threshold is the maximum value minus the interval_
  #threshold percentage selected (default value is 10%)
  selected_values <- as.numeric(names(coverage_per_section[coverage_per_section>=max(coverage_per_section)-round(max(coverage_per_section)*interval_threshold/100)]))

  #Creating a vector with the maximum coverage sections of the 0-1 interval. ONLY THE ACTUAL MAXS 
  output_list[["1-Pur_estimates"]] <- as.numeric(names(coverage_per_section[coverage_per_section>=max(coverage_per_section)]))

  #Getting the maximum coverage intervals fro the data

  #Defining variables
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
  end_val <- ref_val #Assign the ref_val to end_val

  #Add the interval to the interval list
  interval_list <- c(interval_list, list(c(start_val, end_val)))

  output_list[["interval(s)"]] <- interval_list

  #The maximum coverage interval list will be returned
  return(output_list)
  }