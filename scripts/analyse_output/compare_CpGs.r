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

#Determining single classes and combinations to check.

#Creating a list of all teh combinations to analyse
combinations <- list()

for (included_classes in 1:length(names(cpgs_ls))) {

  for (combi in 1:ncol(combn(names(cpgs_ls), included_classes))) {

    combinations <- append(combinations, list(combn(names(cpgs_ls), included_classes)[,combi]))

  }

}


#Adding names to list
names(combinations) <- sapply(1:length(combinations), function(i) {toString(combinations[[i]])})


# Create dataframe to check the cpgs included
context_df <- as.data.frame(matrix(NA, nrow=length(combinations), ncol=length(colnames(cpg_context)[-1])))


colnames(context_df) <- colnames(cpg_context)[-1]
rownames(context_df) <- names(combinations)

cat("\n\nGenerating dataframe...\n")
#Filling the dataframe with the number of CpGs in each context

for (combi in names(combinations)) {

  #Determine the cpgs whose context has ho be checked
  if (length(combinations[[combi]]) == 1) {

    #Using sapply to get a single vector
    cpgs <- sapply(combinations[[combi]], function(x) {cpgs_ls[[x]]})

  } else {

    cpgs <- Reduce(intersect, lapply(combinations[[combi]], function(x) {cpgs_ls[[x]]}))

  }

  # Adding context counts to dataframe
  for (context in colnames(cpg_context)[-1]) {
    context_df[combi, context] <- sum(cpg_context[cpgs, context]==1, na.rm=TRUE)
  }

}

cat("\n\nPlotting...\n")

# Plotting results
for (context in colnames(context_df)) {

  gg <- ggplot(data=context_df, aes(y=context_df[,context], x=as.factor(rownames(context_df)))) +
         geom_bar(stat="identity", fill="lightgreen") +
         xlab("CpGs used in...") + 
         ylab("Number of CpGs") +
         theme_classic() +
         theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(filename=paste(context, ".CpG_context_barplot.png", sep=""))
}
