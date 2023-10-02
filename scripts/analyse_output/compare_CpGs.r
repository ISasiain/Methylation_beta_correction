#!/usr/bin/Rscript

# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

#Loading optparse to parse command line arguments
if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))

#Loading ggplot2 for plotting
if(!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2") }

suppressPackageStartupMessages(library(ggplot2))

#Loading ggplot2 for plotting
if(!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram") }

suppressPackageStartupMessages(library(VennDiagram))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

#Creating argument list
argument_list <- list(

  make_option(c("-c", "--cpgs_to_compare"), type="character",  
              help="The path of all the R objects containing vectors with CpGs to analyse have to be specified here separated by a comma",
              metavar = "[COMMA_SEPARATED_PATHS]"),

  make_option(c("-a", "--context_annotation"), type="character",  
              help="The path to the file containing the annotated CpG context must be entered here.",
              metavar = "[PATH]"),

  make_option(c("-p", "--prefix"), type="character",  
              help="A prefix for the output plots can be entered here. The default prefix is [%output]",
              default="output",
              metavar = "[PREFIX]")

)

# Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This programm gets purity and original beta Robjects from a R object containing those or more elements. It also allow the user to cerate training anjd validation dataset"))


# ============
# LOADING DATA
# ============

# Loading CpG context csv
cpg_context <- read.csv(arguments$context_annotation, row.names=1)

# Loading cpg vectors

#Splitting argument by the comma
path_vec <- strsplit(arguments$cpgs_to_compare, ",")[[1]]

#Creating a list to append cpgs to analyse
cpgs_ls <- list()

#Filling the list
for(file in path_vec) {

    #Generate file name
    name <- strsplit(tail(strsplit(file, "\\/")[[1]], n=1), "_")[[1]][1]

    cpgs_ls[[name]] <- readRDS(file)

}



# ====================
# CHECKING COMMON CPGS
# ====================

#Creating and saving Venn plot
venn.diagram(cpgs_ls,
             filename=paste(arguments$prefix, ".VennDiagram.png", sep=""))




# ============================
# CHECKING CONTEXT OF THE CpGs
# ============================

# Create dataframe to check the cpgs included
context_df <- as.data.frame(matrix(NA, nrow=length(c(names(cpgs_ls)), "ALL"), ncol=length(colnames(cpg_context)[-1])))

#The names of the contexts to eb analysed have bee hardcoded.
colnames(context_df) <- c("promoter","proximal","distal","cgi","shore","ocean", "nonAtac", "atac")
rownames(context_df) <- c(names(cpgs_ls), "ALL")

cat("\n\nGenerating dataframe...\n")

#Filling the dataframe with the number of CpGs in each context

for (type in rownames(context_df)) {

  #Determine the cpgs whose context has to be checked
  if (type == "ALL") {

    #Using sapply to get a single vector
    cpgs <- rownames(cpg_context)

  } else {

    #Using sapply to get a single vector
    cpgs <- cpgs_ls[[x]]

  }

  # Adding context counts to dataframe
  for (context in colnames(context_df)) {
    context_df[type, context] <- sum(cpg_context[cpgs, context]==1, na.rm=TRUE)
  }

}

cat("\n\nPlotting...\n")

for (type in rownames(context_df)) {

  # Create pie chart for promoter, proximal, and distal
  gg_prom_prox_dis <- ggplot(data=context_df, aes(x="", y=as.vector(context_df[type, c("promoter","proximal","distal")]), fill=c("promoter","proximal","distal"))) +
         geom_bar(width=1, stat="identity") +
         labs(title=paste(type, ": promoter, proximal and distal", sep="")) +
         theme_classic() +
         coord_polar("y", start=0)

  ggsave(filename=paste(type, "prom_prox_dis.CpG_context_piechart.png", sep=""))
  
  # Create pie chart for cgi, shore, and ocean
  gg_cgi_shore_ocean <- ggplot(data=context_df, aes(x="", y=as.vector(context_df[type, c("cgi","shore","ocean")]), fill=c("cgi","shore","ocean"))) +
         geom_bar(width=1, stat="identity") +
         labs(title=paste(type, ": cgi, shore and ocean", sep="")) +
         theme_classic() +
         coord_polar("y", start=0)

  ggsave(filename=paste(type, "cgi_shore_ocean.CpG_context_piechart.png", sep=""))

  # Create pie chart for chromatin accessibility
  gg_atac <- ggplot(data=context_df, aes(x="", y=as.vector(context_df[type, c("atac","nonAtac")]), fill=c("atac","nonAtac"))) +
         geom_bar(width=1, stat="identity") +
         labs(title=paste(type, ": chromatin accessibility", sep="")) +
         theme_classic() +
         coord_polar("y", start=0)

  ggsave(filename=paste(type, "atac.CpG_context_piechart.png", sep=""))

}