#!/usr/bin/Rscript

# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

options(repos = "https://cran.r-project.org/")

if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))


# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

argument_list <- list(

  make_option(c("-R", "--ref_regressions"), type="character",  
              help="Path of the directory containing the RObjects with the parameters of the refernce regressions",
              metavar = "[file path]"),

  make_option(c("-p", "--est_purity"), type="character",  
              help="Path to the tsv file with the predicted sample purity values of the samples whose betas have to be corrected. The file must be the tsv text file generated as an output of run_all_validation.r.",
              metavar = "[file path]"),

  make_option(c("-b", "--betas_to_correct"), type="character",  
              help="Path to the file with the beta values to be corrected whose sample purity has been estimated. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.",
              metavar = "[file path]"),

  make_option(c("-o", "--output"), type="character", default="./",
              help="The path to the location where the output files will be saved must be entered here. The output is an R object. Default [%default]",
              metavar = "[file path]"),

  make_option(c("-n", "--output_name"), type="character", default="output",
              help="The prefix to be used to name the output files. Default [%default]",
              metavar = "[file path]")

)

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program corrects methylation beta values providing parameters of the regressions used for the correction."))


# ============
# LOADING DATA
# ============

#Reading the R objects containing the regression data as dataframes
my_slopes <- readRDS(list.files(arguments$ref_regressions, pattern="*reg.slopes.RData", full.names=TRUE))
my_intercepts <- readRDS(list.files(arguments$ref_regressions, pattern="*reg.intercepts.RData", full.names=TRUE))


# Loading betas to correct
to_correct_betas <- readRDS(arguments$betas_to_correct)

#Adapt name of betas
colnames(to_correct_betas) <- lapply(colnames(betas), 
                         function (name) {strsplit(name, "-01")[[1]][1]})

# Loading predicted 1 - purities
predicted_1mPurities <- read.table(arguments$est_purity, 
                                header=TRUE,
                                sep="\t")

# Removing samples with more than one estimates (if any)
predicted_1mPurities <- predicted_1mPurities[which(predicted_purities[,"num_of_est"]==1),]

# Transforming the predicted_1mPurities dataframe into a vector
predicted_1mPurities_vec <- predicted_1mPurities[,"estimate"]
names(predicted_1mPurities_vec) <- predicted_1mPurities[,"sample"]



# ================
# FILTERING CPGS
# ================

cat("\nChecking cpgs...\n\n")

if (length(!(rownames(to_correct_betas) %in% rownames(my_slopes))) != 0) {

    # Printing warning message
    cat("\n",  length(!(rownames(to_correct_betas) %in% rownames(my_slopes))), "CpGs are not included into the refernce regressions, so they can not be corrected.\n\n")

    # Filtering not included CpGs
    to_correct_betas <- to_correct_betas[rownames(to_correct_betas) %in% rownames(my_slopes),]
}

# Remove CpGs from the regressions that are not included into the data to correct to speed up the process.
my_slopes <- my_slopes[rownames(my_slopes) %in% rownames(to_correct_betas),]
my_intercepts <- my_intercepts[rownames(my_intercepts) %in% rownames(to_correct_betas),]



# ==============================================
# CORRECTING BETAS BASED ON REFERNCE REGRESSIONS
# ==============================================

#Generating function to identify the refernce regression to which each CpG of each sample belongs
identify_regression <- function(beta, estimated_1mPurity, vec_slopes, vec_intercepts) {

    # Generating a vector of distances
    vec_of_dis <- sapply(c(1:length(vec_slopes)), 
                         FUN = function(pop) {

                                    #Ignoring NA values that appear when less than 3 populations are detected
                                    #in the refernce rrgressions
                                    if (!is.na(vec_slopes[pop])) {

                                        vec_slopes[pop] * estimated_1mPurity + vec_intercepts[pop] - beta

                                    }
                                }
                            )

    #Determining the population (vector index) with the lowest absolute distance. If the distances are equal the first
    #population with be chosen by default
    pop_identified <- which.min(abs(vec_of_dis))

    #Generating and returning output named vector
    return(
        c("Slope"=vec_slopes[pop],
          "Itercepts"=vec_intercepts[pop],
          "Distance"=vec_of_dis[pop])
        )
    }


#Generate function to correct betas based on the identified regression. The parametres specified must be 
#from the beta VS 1-P regressions.
correcting_betas <- function(slope, intercept, distance, to_correct) {

    if (to_correct=="Tumor") {
        #The tumor beta value will be obtained using the intercept and the calculated distance.
        #The maximum possible value will allways be kept below or equal to 1
        tum_beta <- if(intercept + distance <= 1) {
                        intercept + distance
                    } else {
                        1
                    }

        return(tum_beta)

    } else if (to_correct=="Micoenvironment") {
        #The microenvironment beta value will be obtained using the intercept and slope when 1-P=1 and the calculated distance.
        #The minimum possible value will allways be kept below or equal to 1
        env_beta <- if (intercept + slope + distance >= 0) {
                        intercept + slope + distance
                    } else {
                        0
                    }
        
        return(env_beta)

    }
}

 # Creating dataframes to store the corrected values
 corrected_tumor <- data.frame()
 corrected_microenvironment <- data.frame()



 for (cpg in rownames(to_correct_betas)) {

    corrected_tumor[cpg,rownames(to_correct_betas)] <- sapply(names(predicted_1mPurities), 

        function(sample) {

        identified_reg <- identify_regression(
            beta=to_correct_betas[cpg, sample],
            estimated_1mPurity=predicted_1mPurities[sample],
            vec_slopes=my_slopes[cpg,],
            vec_intercepts=my_intercepts[cpg,]
            )

        corrected_tum <- correcting_betas(
            slope=identified_reg["Slope"],
            intercept=identified_reg["Intercept"],
            distance=identified_reg["Distance"],
            to_correct="Tumor"
            )

        return(corrected_tum)

        }
    )

    corrected_microenvironment[cpg,rownames(to_correct_betas)] <- sapply(names(predicted_1mPurities), 

        function(sample) {

        identified_reg <- identify_regression(
            beta=to_correct_betas[cpg, sample],
            estimated_1mPurity=predicted_1mPurities[sample],
            vec_slopes=my_slopes[cpg,],
            vec_intercepts=my_intercepts[cpg,]
            )

        corrected_env <- env_beta(
            slope=identified_reg["Slope"],
            intercept=identified_reg["Intercept"],
            distance=identified_reg["Distance"],
            to_correct="Microenvironment"
            )

        return(corrected_env)

        }
    )
 }



# =======================
# GENERATING OUTPUT FILES
# =======================

# Generating RObject files
saveRDS(corrected_tumor, file=paste(arguments$output_location, arguments$output_filename, ".tumor.samples_to_correct.RData", sep=""))
saveRDS(corrected_microenvironment, file=paste(arguments$output_location, arguments$output_filename, ".microenvironment.samples_to_correct.RData", sep=""))

# Generating tsv files
write.table(corrected_tumor, 
            file=paste(arguments$output_location, arguments$output_filename, ".tumor.samples_to_correct.tsv", sep=""),
            col.names=NA, 
            sep="\t")
write.table(corrected_microenvironment, 
            file=paste(arguments$output_location, arguments$output_filename, ".microenvironment.samples_to_correct.tsv", sep=""),
            col.names=NA, 
            sep="\t")