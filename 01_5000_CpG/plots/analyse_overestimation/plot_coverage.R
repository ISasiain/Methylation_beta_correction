# ==================================
# INSTALLING AND/OR LOADING PACKAGES
# ==================================

if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", quietly = TRUE) }

suppressPackageStartupMessages(library(optparse))

if(!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", quietly = TRUE) }

suppressPackageStartupMessages(library(ggplot2))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

argument_list <- list(

  make_option(c("-i", "--input"), type="character",   
              help="Name (path) of the input file",
              metavar= "[path to data]"),

  make_option(c("-t", "--title"), type="character",
              help="The title of the plot must be entered here",
              metavar= "[Plot's title]"),
  
  make_option(c("-o", "--output"), type="character", 
              help="Name of the output plot (without extension)",
              metavar= "[name of plot]")  

)

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program plots the coverage of different samples together"))


# ==================
# READING INPUT FILE
# ==================

file_to_plot <- readRDS(arguments$input)
plot_created <- FALSE
my_col <- 0
all_cols <- rainbow(length(file_to_plot))

for (sample in names(file_to_plot)) {

  my_col <- my_col + 1

  if (!plot_created)  {
    
    my_plot <- ggplot(data=data.frame(file_to_plot[[sample]][["Coverage_plot"]]),
                      aes(x=x, y=y, color = all_cols[my_col])) +
                      geom_line(color = all_cols[my_col],linewidth=0.5, alpha=0.5) +
                      ggtitle(arguments$title) +
                      xlab("1-Purity") +
                      ylab("Coverage") +
                      theme_classic() +
                      theme(plot.title = element_text(size = 20),
                        axis.title = element_text(size = 16),
                        axis.text = element_text(size = 14),
                        panel.grid.major = element_line(colour = "lightgrey", linetype = "dotted"),
                        panel.grid.minor = element_blank())
                      

    
    plot_created <- TRUE
    
  } else {
    
    my_plot <- my_plot + geom_line(data=data.frame(file_to_plot[[sample]][["Coverage_plot"]]),
                                   aes(x=x, y=y), 
                                   color = all_cols[my_col],
                                   linewidth=0.5, alpha=0.5)
    
  }
  
}


ggsave(paste(arguments$output, ".png", sep=""))