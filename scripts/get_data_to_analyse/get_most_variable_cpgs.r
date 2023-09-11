# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

#Loading optparse to parse command line arguments
if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))


# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

#Creating argument list
argument_list <- list(

  make_option(c("-r", "--reference_dataset"), type="character",  
              help="The path to the RObject containing the reference data (betas of each CpG) must be entered here.",
              metavar = "[PATH_TO_REFERNCE]"),

  make_option(c("-a", "--data_to_analyse"), type="character",  
              help="The path to the RObject containing the data to analyse (betas of each CpG) must be entered here.",
              metavar = "[PATH_TO_DATA_TO_ANALYSE]"),

  make_option(c("-n", "--number_of_CpGs"), type="numeric",  
              help="The nomber of the most variable CpGs to keep must be entered here.",
              metavar = "[CpG_NUM]"),

  make_option(c("-l", "--include_CpG_list"), type="logical",  
              help="This argument must be set to TRUE to generate a file containing the included CpG names as an output. The default argument is [%default]",
              default=FALSE,
              metavar = "[TRUE/FALSE]"),

  make_option(c("-p", "--prefix"), type="character",  
              help="The prefix of the output files can be selected here. The default argument is [%default]",
              default="output",
              metavar = "[PREFIX]")
)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This programm gets purity and original beta Robjects from a R object containing those or more elements. It also allow the user to cerate training anjd validation dataset"))


# ===========
# LODING DATA
# ===========

# Reading the reference betas
reference <- readRDS(arguments$reference_dataset)

# Reading the data to analyse
to_analyse <- readRDS(arguments$data_to_analyse)


# =================================
#    GETTING MOST VARIABLE CPGs
# =================================

# Detemine the variance of all the rows (CpGs)
 cpgs_variance <- apply(
                        reference,
                        MARGIN=1,
                        FUN=var
                   )

# Create a vector to sort the rows. Get only the number of rows containing the CpGs that want to be included
sorting_vec <- order(cpgs_variance, decreasing=TRUE)[1:arguments$number_of_CpGs]

# Sorting the betas dataframe and getting only the CpGs to include
cpgs <- rownames(reference[sorting_vec,])

# ===========================
#     GENERATING OUTPUT
# ===========================

to_analyse <- to_analyse[cpgs,]

saveRDS(to_analyse, file=paste(arguments$prefix, "_most_variable_CpGs.RData", sep=""))

if (arguments$include_CpG_list) {

    saveRDS(cpgs, file=paste(arguments$prefix, "_CpG_vector.RData", sep="")) 

}