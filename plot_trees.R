#!/usr/bin/env Rscript

############################################################
# General purpose tree plotting script for consistency.    #
# This script can be slow as it was designed to be pretty  #
# bullet-proof when loading/installing necessary packages. #
#                                                          #
# Rendering may need to be adjusted for speed.             #
#                                                          #
# By J. Healey                                             #
############################################################

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")

# ggtree is installed by:
#biocLite("ggtree")

# Biostrings is installed by:
##biocLite("Biostrings")

## biocLite("BiocUpgrade") ## you may need this if you get version
                           ## errors from Bioconductor

# Standard install if missing
list.of.packages <- c("ggplot2", "argparse","ape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Bioconductor install if missing
list.of.bioc.packages <- c("Biostrings", "ggtree")
new.bioc.packages <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]
if(length(new.bioc.packages)) biocLite(new.bioc.packages)

suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
suppressMessages(library("ape"))
suppressMessages(library("Biostrings"))
suppressMessages(library("ggtree"))
suppressMessages(library("tools"))
options(warn=-1)

# Parse commandline arguments

parser <- ArgumentParser()

parser$add_argument('-t',
                    '--tree',
                    action='store',
                    required=TRUE,
                    help="Tree file for plotting with ggtree in Newick format.")

parser$add_argument('-v',
                    '--verbose',
                    action='count',
                    default=0,
                    help='Verbose behaviour, printing up to 3 levels of different output. Default off.')

parser$add_argument('-o',
                     '--outfile',
                     action='store',
                     default='None',
                     help='Filename and path to save the image to. Defaults to the same as the infile with .png appended')

args <-parser$parse_args()

# Reassign to standard variables to pass through functions
tree <- args$tree
outfile <- args$outfile
verbose <- args$verbose

# Function to synthesise an output path/name
getOutfile <- function(tree){
  outfile <- paste(paste(dirname(tree),
                         file_path_sans_ext(basename(tree)),
                         sep="/"),
                  "png",
                   sep = ".")
  return(outfile)
}

# If no output file specified, redefine outfile to infile + ".png"
if (outfile == 'None'){
  outfile <- getOutfile(tree)
}

if (verbose > 0){
  cat("\n", "Reading in tree:", "\n")
  cat("=================", "\n")
  cat(tree, "\n")
  
}

# Import tree data
tree.obj <- read.tree(args$tree)

if (verbose > 1){
  cat("\n", "Your data structure:", "\n")
  cat("=====================", "\n","\n")
  print.phylo(tree.obj, printlen = length(tree.obj$tip.label))
  cat("\n")
  write.tree(tree.obj)
  
  }

# Plot and customise tree space
tree.img <- ggtree(tree.obj, size=0.75) + xlim(NA, 1.4)

# Add tip labels (right aligned)
tree.img <- tree.img + geom_tiplab(align=T,
                                   linesize = 0.3,
                                   offset = 0.1,
                                   size = 8)

# Add bootstraps/node labels
tree.img <- tree.img + geom_text2(size=5,
                                  aes(subset = !isTip, label=label),
                                  nudge_x = 0.02)

# Add tip shapes
tree.img <- tree.img + geom_tippoint(color='firebrick')

# Add scalebar top left
tree.img <- tree.img + geom_treescale(fontsize=5,
                                      linesize = 2,
                                      offset = -1)

# Make transparent
tree.img <- tree.img + theme(plot.background = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank(),
#                             panel.grid = element_blank(),
                             legend.background = element_blank(),
                             plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))



tree.img

if (verbose > 0){
  cat("\n", "Saving tree image to:", "\n")
  cat("======================", "\n")
  cat(outfile, "\n")
}

ggsave(outfile,
       plot = tree.img,
       height = 8,
       width = 12,
       units = "in",
       dpi = 1200, 
       bg = "transparent")


