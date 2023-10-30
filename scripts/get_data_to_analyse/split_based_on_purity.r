# Script Name: split_based_on_purity.r

# Description:
# This script is designed to split beta values and purity data into eight intervals based on specified purity intervals. It loads beta values and purity data, then creates separate data subsets for each purity interval.
# Important: The purity_validation.r and betas_validation.r scripts must be located a directory before the workin directory.

# Loading Data:
# - The script reads beta values from an RData file ("betas_validation.RData") and purity data from another RData file ("purity_validation.RData").

# Creating Purity Intervals:
# - The script defines a list, "list_to_split," containing eight purity intervals with specified purity range values.
# - Each interval is named from "from_1_to_0.7" to "from_0.2_to_0."

# Iterating through Purity Intervals:
# - The script iterates through the "list_to_split," creating data subsets for each purity interval.
# - For each interval, it filters and assigns the beta values and purity values corresponding to that interval.
# - The filtered beta values and purity values are saved as separate RData files with specific names, such as "betas_from_1_to_0.7.RData" and "pur_from_1_to_0.7.RData."

# Example Usage:
# - Users can run the script to split beta values and purity data into specific purity intervals as defined in the "list_to_split."

# Author: Iñaki Sasiain Casado
# Affiliation: Johan Staaf lab


#Loading the data
betas <- readRDS("../betas_validation.RData");
purities <- readRDS("../purity_validation.RData");

#Creating a list with the 8 purity intrevals to split the data into (The names are 1-Purity, not purities)
list_to_split <- list(
    from_1_to_0.7 = c(0, 0.3),
    from_0.7_to_0.6 = c(0.3, 0.4),
    from_0.6_to_0.5 = c(0.4, 0.5),
    from_0.5_to_0.4 = c(0.5, 0.6),
    from_0.4_to_0.35 = c(0.6, 0.65),
    from_0.35_to_0.3 = c(0.65, 0.7),
    from_0.3_to_0.2 = c(0.7, 0.8),
    from_0.2_to_0 = c(0.8, 1)
)

#Iterate through the list to create a vector with the sample names included in each interval
for (interval in names(list_to_split)) {
    
    assign(paste("pur_",interval,sep=""),purities[purities >= list_to_split[[interval]][1] & purities < list_to_split[[interval]][2]])
    assign(paste("betas_",interval,sep=""), betas[,names(purities[purities >= list_to_split[[interval]][1] & purities < list_to_split[[interval]][2]])])

    saveRDS(get(paste("betas_",interval,sep="")), file=paste("betas_",interval,".RData",sep=""))
    saveRDS(get(paste("pur_",interval,sep="")), file=paste("pur_",interval,".RData",sep=""))
}