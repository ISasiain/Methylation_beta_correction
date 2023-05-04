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