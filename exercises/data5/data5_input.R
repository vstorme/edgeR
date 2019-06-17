## ------------------------------------------------------------------------
sessionInfo()

#' 
#' this RNAseq experiment involves a control-treatment case
#' there are 3 independent biological samples for each treatment group
#' 
#' start by setting the working directory
#' 
#' case 2.6 from manual: 1 factor with 2 factor levels and batch effect
#' ---------------------------------------------------------------------
#' deltahrcC challenged and mock-inoculated samples
#' samples were collected in 3 batches
#' 
#' set the working dir
#' 
## ------------------------------------------------------------------------
setwd("C:/MYCOURSES/RNAseq_v2/EXERCISES/data5")

#' 
#' when working from Athena type:
#' > setwd("//Client/D$/MYCOURSES/RNAseq_v2/EXERCISES/data5")
#' 
#' load the necessary libraries
#' 
## ------------------------------------------------------------------------

source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
biocLite("NBPSeq")

library(knitr)
library(limma)
library(edgeR)
library(locfit)
library(statmod)
library(Rcpp)
library(pheatmap)
library(gplots)
library(qvalue)
library(NBPSeq)

#' 
#' reading the data
#' -------------------------------------
#' 
#' load the data from the NBPSeq package
#' 
## ------------------------------------------------------------------------
data(arab)
head(arab)

#' batch effect is caused by the time that each replicate was constructed
#' use "substring" to extract a substring in a character vector
#' usage substring(x,start,stop)
#' use relevel to re-order so that the level specified by ref is first and the others are moved down
#' 
## ------------------------------------------------------------------------
Treat <- factor(substring(colnames(arab),1,4))
Treat <- relevel(Treat, ref="mock")
Time <- factor(substring(colnames(arab),5,5))
