setwd("/home/isc/Methylation/adjustBetas/12_plots_for_methods/sample_TCGA-EW-A1P7-01A")

#Sourcing functions
source("../../scripts/calculate_purity/predicting_purity.r")
source("../../scripts/calculate_purity/purity_coverage.r")

#Loading data
betas <- readRDS("bet_TCGA-EW-A1P7-01A.RData")
purities <- readRDS("pur_TCGA-EW-A1P7-01A.RData")
my_slopes <- readRDS(list.files("./regressions", pattern="*reg.slopes.RData", full.names=TRUE))
my_intercepts <- readRDS(list.files("./regressions", pattern="*reg.intercepts.RData", full.names=TRUE))
my_RSE <- readRDS(list.files("./regressions", pattern="*reg.RSE.RData", full.names=TRUE))
my_df <- readRDS(list.files("./regressions", pattern="*reg.df.RData", full.names=TRUE))

#Running the functions (some lines of purity_coverage.r were uncommented to get and save the plots)

# Defining an empty matrix with the cpg ids as rownames to add the all the 1-Purity predicted intervals for all 
# the CpGs of a sample
interval_mat <- matrix(ncol=2, nrow=length(betas))
rownames(interval_mat) <- names(betas)

# Predicting all the 1-Purity intervals for each CpG of each sample and append them to the empty interval_mat
for (cpg in names(betas)) {

    # The following if statement will be used to take into account only cpgs included into the
    # refernce regression dataset
    if (cpg %in% rownames(my_slopes)) {

      interval_mat[cpg,] <- predicting_purity(beta=betas[cpg],
                                              slopes=my_slopes[cpg, ],
                                              intercepts=my_intercepts[cpg, ],
                                              RSE=my_RSE[cpg, ],
                                              degrees_of_freedom=my_df[cpg, ],
                                              slope_threshold=0.2,
                                              alpha=0.7)
    
    }
  }

purity_value_per_sample(
    interval_mat,
    5
)