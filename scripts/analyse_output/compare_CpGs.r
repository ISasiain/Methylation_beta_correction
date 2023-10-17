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
              help="The path to the RObject file containing the annotated CpG context must be entered here.",
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

cat("\n\nLoading the data...\n")

# Loading CpG context csv
cpg_context <- readRDS(arguments$context_annotation)

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

cat("\n\nGenerating Venn diagram...\n")

#Creating and saving Venn plot
venn.diagram(cpgs_ls,
             filename=paste(arguments$prefix, ".VennDiagram.png", sep=""))




# ============================
# CHECKING CONTEXT OF THE CpGs
# ============================

# Create dataframe to check the cpgs included
context_df <- as.data.frame(matrix(NA, nrow=length(c(names(cpgs_ls), "ALL")), ncol=length(colnames(cpg_context))))
rownames(context_df) <- c(names(cpgs_ls), "ALL")
colnames(context_df) <- colnames(cpg_context)


cat("\n\nGenerating dataframe...\n")

#Filling the dataframe with the number of CpGs in each context

for (type in rownames(context_df)) {

  #Determine the cpgs whose context has to be checked
  if (type == "ALL") {

    #Using sapply to get a single vector
    cpgs <- rownames(cpg_context)

  } else {

    #Using sapply to get a single vector
    cpgs <- cpgs_ls[[type]]

  }

  # Adding context counts to dataframe
  for (context in colnames(context_df)) {
    context_df[type, context] <- sum(cpg_context[cpgs, context]==1, na.rm=TRUE)
  }

}



cat("\n\nPlotting...\n")

for (type in rownames(context_df)) {

  # Create a dataframe per each type to be plotted 
  type_df <- as.data.frame(cbind(as.numeric(context_df[type,]), as.character(colnames(context_df))))
  colnames(type_df) <- c("count", "context")
  rownames(type_df) <- colnames(context_df)
  
  type_df$"count" <- as.numeric(type_df$"count")



  # Create pie chart for promoter, proximal, and distal
  prom_prox_dis <- type_df[c("promoter","proximal","distal"),]

  gg_prom_prox_dis <- ggplot(data=prom_prox_dis, aes(x="", y=count, fill=context)) +
         geom_bar(stat="identity") +
         labs(title=paste(type, ": promoter, proximal and distal", sep="")) +
         theme_void() +
         coord_polar("y", start=0) +
         geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

  ggsave(filename=paste(type, "_prom_prox_dis.CpG_context_piechart.png", sep=""))
  
  # Create pie chart for cgi, shore, and ocean
  cgi_shore_ocean <- type_df[c("cgi","shore","ocean"),]

  gg_cgi_shore_ocean <- ggplot(data=cgi_shore_ocean, aes(x="", y=count, fill=context)) +
         geom_bar(stat="identity") +
         labs(title=paste(type, ": cgi, shore and ocean", sep="")) +
         theme_void() +
         coord_polar("y", start=0) +
         geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

  ggsave(filename=paste(type, "_cgi_shore_ocean.CpG_context_piechart.png", sep=""))


  if (type!="ALL") {
    # Create pie chart for chromatin accessibility
    atac <- type_df[c(paste("atac", type, sep=""),paste("nonAtac", type, sep="")),]

    gg_atac <- ggplot(data=atac, aes(x="", y=count, fill=context)) +
          geom_bar(stat="identity") +
          labs(title=paste(type, ": chromatin accessibility", sep="")) +
          theme_void() +
          coord_polar("y", start=0) +
          geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

    ggsave(filename=paste(type, "_atac.CpG_context_piechart.png", sep=""))

  } else {

    #All BRCA
    atac <- type_df[c("atacBRCA","nonAtacBRCA"),]

    gg_atac <- ggplot(data=atac, aes(x="", y=count, fill=context)) +
          geom_bar(stat="identity") +
          labs(title=paste(type, ": chromatin accessibility", sep="")) +
          theme_void() +
          coord_polar("y", start=0) +
          geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

    ggsave(filename=paste("ALL_BRCA_atac.CpG_context_piechart.png", sep=""))


    #All LUSC
    atac <- type_df[c("atacLUSC","nonAtacLUSC"),]

    gg_atac <- ggplot(data=atac, aes(x="", y=count, fill=context)) +
          geom_bar(stat="identity") +
          labs(title=paste(type, ": chromatin accessibility", sep="")) +
          theme_void() +
          coord_polar("y", start=0) +
          geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

    ggsave(filename=paste("ALL_LUSC_atac.CpG_context_piechart.png", sep=""))



    #All LUAC
    atac <- type_df[c("atacLUAC","nonAtacLUAC"),]

    gg_atac <- ggplot(data=atac, aes(x="", y=count, fill=context)) +
          geom_bar(stat="identity") +
          labs(title=paste(type, ": chromatin accessibility", sep="")) +
          theme_void() +
          coord_polar("y", start=0) +
          geom_text(aes(label = count), position = position_stack(vjust = 0.5))  # Add labels

    ggsave(filename=paste("ALL_LUAC_atac.CpG_context_piechart.png", sep=""))

  }
}

cat("\n=================\n")
cat ("PROCESS FINISHED")
cat("\n=================\n")