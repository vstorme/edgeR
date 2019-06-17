#'   
#' ---
#' title: "case study 1: two group comparison"
#' author: "Veronique Storme"
#' date: "14/11/2018"
#' output: html_document
#' ---
#' 
#' 
#' bioconductor must have been installed as well as the edgeR package
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("BiocUpgrade")
#' biocLite("edgeR")
#' biocLite("limma")
#' 
#' The locfit package is also needed
#' biocLite("locfit")
#' 
#' other packages are needed for plotting and clustering, eg ggplo2
#' install.packages("ggplot2")
#' 
#' load the necessary libraries
#' 
library(knitr)
library(limma)
library(edgeR)
library(locfit)
library(statmod)
library(Rcpp)
library(pheatmap)
library(gplots)
library(tidyr)

#' 
#' getting the version and correct citations
#' 
sessionInfo()
citation("base")
citation("edgeR")

#' 
#' this RNAseq experiment involves a control-treatment case
#' there are 3 independent biological samples for each treatment group
#' 
#' ### reading the counts from separate files
#' 
#' start by setting the working directory
#' 
## ------------------------------------------------------------------------
setwd("//client/C$/MYCOURSES/RNAseq_v2/DEMO/data1")

setwd("C:/MYCOURSES/RNAseq_v2/DEMO/data1")


#' 
#' counts are stored in 6 separate plain text files
#' 
## ------------------------------------------------------------------------
dir()

#' 
#' the file target_data1.txt gives the filename, the group and a brief description for each sample
#' 
## ------------------------------------------------------------------------
targets <- readTargets("target_data1.txt")
targets

#' 
#' read the tables of counts, calculate the sizes of the libraries and produce a DGElist object
#' 
## ------------------------------------------------------------------------
d <- readDGE(targets)
class(d)
names(d)
d$samples
head(d$counts)
dim(d)

#' 
#' ### data exploration
#' 
#' (asinh is kind of log transformation but handles the zeroes better)
#' 
## ------------------------------------------------------------------------
hist(d$counts)
hist(asinh(d$counts))

#' 
## ------------------------------------------------------------------------
barplot(d$samples$lib.size*1e-6, names=1:6, ylab="Library size (millions)",cex.names=0.5)

#' 
#' get the counts in cpm (=counts divided by the libsize and multiplied by a million)
#' cpm(x, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25, ...)
#' prior.count: average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE
#' ?cpm
#' 
## ------------------------------------------------------------------------
d$samples
d$counts[1:5,]
cpm(d)[1:5,]
?cpm

#' 
#' ### filtering
#' 
#' filter low expression tags (little power to detect DE for filtered tags)
#' 
#' a gene is required to have a count of 5-10 in a library to be considered expressed in that library
#' the expression must exceed 5-10 read counts in at 1 condition (if there are 3 replicates for each condition, then in at least 3 samples)
#' libsize divided by 1 million gives an idea about the relationship between raw count and cpm
#' thus 5 read counts correspond to 5 divided by the libsize and multiplied by a million
#' 
## ------------------------------------------------------------------------
d$counts[1:5,]
prop <- d$samples$lib.size/1000000
prop # nr of reads corresponding to 1 cpm
cutoff = 5/prop
cutoff

#' 
#' this means that 5 read counts correspond to 0.5 to 1.12 cpm
#' filter on 1 cpm
#' 
## ------------------------------------------------------------------------
keep <- rowSums(cpm(d)>1) >=3
df <- d[keep,]
dim(df)

#' 
#' cpm(d) >= 3 returns a dataframe with TRUE and FALSE
#' 
## ------------------------------------------------------------------------
cpm(d)[1:10,]
eg = cpm(d) > 1
eg[1:10,]
keep[1:10]

#' 
#' thus the command "rowSums(cpm(d)>1) >=3" keeps those rows where the expression is true in at least 3 samples 
#' 
#' reset library sizes (must be done manually, using the colSums function)
#' 
## ------------------------------------------------------------------------
df$samples$lib.size <- colSums(df$counts)
df$samples

#' 
## ------------------------------------------------------------------------
hist(asinh(df$counts))

#' 
#' ### Normalization
#' 
#' apply TMM normalisation (TMM is the default method)
#' 
#' usage:
#' calcNormFactors(object, method=c("TMM","RLE","upperquartile","none"), refColumn = NULL,
#'       logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
#' refColumn: column (sample) to use as reference for method="TMM".
#' If refColumn is unspecified, the library whose upper quartile is closest to the mean upper quartile is used.
#' 
## ------------------------------------------------------------------------
df <- calcNormFactors(df,refColumn=1)
df$samples

#' 
#' Note on cpm:
#' When TMM normalisation has been performed, then value for norm.factors, is then also used in cpm calulation
#' 
#' the effective library size = product of the original library
#' size (lib.size) and the scaling factor (norm.factors)
#' 
## ------------------------------------------------------------------------
lib.size.eff <- df$sample$lib.size*df$sample$norm.factors
lib.size.eff
# compare with
df$samples$lib.size

# cpm.sample3 <- (df$counts[,3]/lib.size.eff[3])*1000000
# cpm.sample3[1:5]
# cpm(df, normalized.lib.sizes=TRUE)[1:5,3]

#' 
#' The performance of the TMM normalization procedure can be examined using mean-difference (MD) plots.
#' The expression profiles of individual samples can be explored more closely with mean-difference (MD) plots. An MD plot visualizes the library size-adjusted log-fold change between two libraries (the difference) against the average log-expression across those libraries (the mean). 
#' 
#' The following command produces an MD plot that compares sample 1 to an artificial reference library constructed from the average of all the other samples:
#' 
## ------------------------------------------------------------------------
plotMD(df, column=1)
abline(h=0, col="red", lty=2, lwd=2)

#' 
#' Ideally, the bulk of genes should be centred at a log-fold change of zero. This
#' indicates that any composition bias between libraries has been successfully removed.
#' This quality check should be repeated by constructing a MD plot for each sample
#' 
## ------------------------------------------------------------------------
for (i in 1:6)
  {
  jpeg(file=paste("MDplot_sample", i, ".jpg", sep=""),width = 10, height = 10, 
       units = "cm", res = 200) 
  plotMD(df, column = i)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
  }

#' 
#' comparing densities
#' 
## ------------------------------------------------------------------------

# reshape the count data from wide to long
head(df$counts)

countdata = as.data.frame(df$counts)

longdata = gather(countdata, key=sample, value=counts, 1:6, factor_key=TRUE)
longdata[16940:16950,]

library(ggplot2)
ggplot(longdata,aes(x=asinh(counts), fill=sample)) + geom_density(alpha=0.25)

#' 
#' 
#' ### multidimensional scaling plot
#' 
#' The data can be explored by generating multi-dimensional scaling (MDS) plots.
#' This visualizes the differences between the expression profiles of different samples in two dimensions.
#' 
#' The distance between each pair of samples in the MDS plot is calculated as the leading fold change, defined as the root-mean-square of the largest 500 log2-fold changes between that pair of samples (default).
#' 
## ------------------------------------------------------------------------

plotMDS(df,top = 500, gene.selection="pairwise",labels = d$samples$description, cex=0.5)

#' 
#' ### dispersion estimation
#' 
#' > define the design of the experiment
#' 
## ------------------------------------------------------------------------

df$samples
levels(df$samples$group)
df$samples$group <- relevel(df$samples$group, ref="control")

design = model.matrix(~ group, data=df$samples)

design
rownames(design) <- df$samples$group
design
#' 
#' > non-robust estimate of the prior degrees of freedom
#' 
#' common, trended and tagwise dispersions are estimated
#' 
#' When robust=FALSE, there is only one prior.df estimated, when robust=TRUE, each gene gets a prior.df value 
#' Ref: Phipson et al. 2013 "empirical Bayes in the presence of exceptional cases with application to microarray data"
#' 
## ------------------------------------------------------------------------
df <- estimateDisp(df,design=design, robust=FALSE)
df$prior.df

#' 
#' > robust estimate of the prior degrees of freedom
#' 
#' robust=TRUE : highly variable genes are less squeezed to the trend, resulting in a lower estimate of prior.df (=weight given to the trend)
#' 
## ------------------------------------------------------------------------
dfr <- estimateDisp(df,design=design, robust=TRUE)
dfr$prior.df[1:5]

#' 
#' The following plot shows the biological coefficient of variation for each gene against its abundance with the estimated common dispersion and trended dispersions
#' 
## ------------------------------------------------------------------------
plotBCV(dfr,cex=0.4)

#' 
#' The value of the trended NB dispersions should range between 0.005 to 0.05 for laboratory-controlled biological systems like mice or cell lines, though larger values will be observed for patient- derived data (> 0.1) or single-cell data (> 1). 
#' Note that tagwise and common estimates are also shown here but will not be used further.
#' 
#' outliers are marked by low prior.df values
#' 
## ------------------------------------------------------------------------
quantile(dfr$prior.df,c(0,0.005, 0.01,0.05,0.5,0.95,0.99,0.995, 1))

#' 
#' 0.5% of the genes have prior.df values < 1.45
#' 
#' outliers with a prior.df value < 1
#' 
## ------------------------------------------------------------------------
outliers = dfr$counts[dfr$prior.df < 1,]
dim(outliers)

round(cpm(outliers), digits=0)[1:10,]

#' 
#' ### fitting a quasi-likelihood model
#' 
## ------------------------------------------------------------------------
fit = glmQLFit(dfr, design, robust=TRUE)

#' 
#' ### Testing for DE genes
#' 
#' The model fitted is:
#' 
#' log(E(Y/efflib)) = b0 + b1X (X=1 when group=treated and 0 when group=control)
#' 
#' verify with
#' 
## ------------------------------------------------------------------------
design

#' 
#' H0: log(E(Y/efflib|treated) = log(E(Y/efflib|control) 
#' 
#' or log(E(Y/efflib|treated) - log(E(Y/efflib|control) = 0
#' 
#' log(E(Y/efflib)|treated) = b0 + b1
#' 
#' log(E(Y/efflib)|control) = b0
#' 
#' Thus log(E(Y/efflib|treated) - log(E(Y/efflib|control) = b1 = 0
#' 
#' Thus we need to test whether b1=0
#' 
## ------------------------------------------------------------------------

DE = glmQLFTest(fit, coef=2)
head(DE$table)
head(fit$coefficients)

#' 
#' Note: logFC is here log2FC
#' 
#' topTags() produces a table of the top differentially expressed tags (default n=10)
#' default adjustment for multiple testing with the BH method (FDR column)
#' sorted by (unadjusted) pvalue (default) or by absolute log2-fold change or none
#' 
## ------------------------------------------------------------------------

DE.fdr = topTags(DE,n=Inf,sort.by="none")
names(DE.fdr)
head(DE.fdr$table)

#' 
#' export as csv file
#' 
## ------------------------------------------------------------------------
write.csv(DE.fdr$table,file="de_fdr_logFC.csv")

#' 
#' > use of decideTestsDGE function
#' 
#' total nr of DE genes at FDR < 0.05 is :
#' 
## ------------------------------------------------------------------------
de <- decideTestsDGE(DE, p=0.05, adjust="BH")
de[1:5,]
summary(de)

#' 
#' entries for -1, 0, and 1  are for down-regulated, nonDE and up-regulated tags, respectively
#' 
#' Note that with the previous non-robust nonQL method, 15 genes were found to be DE and upregulated
#' 
#' The QL method with the robust option has a better control of the FDR
#' 
#' ### using a different parameterization: the no intercept model
#' 
## ------------------------------------------------------------------------

group = c(rep("control",3),rep("treated",3))
design = model.matrix(~ 0 + group)
rownames(design) <- group
design

#' 
#' re-estimate dispersions
#' 
## ------------------------------------------------------------------------
dfr <- estimateDisp(df,design=design, robust=TRUE)

#' 
#' re-fit the model
#' 
## ------------------------------------------------------------------------
fit = glmQLFit(dfr, design, robust=TRUE)
head(fit$coefficients)

#' 
#' Now we have:
#' 
#' log(E(Y/efflib)) = b1X1 + b2X2 
#' (X1=1 when group=control, X2=1 when group=treated)
#' for geneID00001 : log(E(Y/efflib)) = -12.06X1 -12.10X2
#' log(E(Y/efflib|control)) = -12.06
#' log(E(Y/efflib|treatment)) = -12.10
#' 
#' testing for DE
#' 
#' now we need a different contrast
#' 
## ------------------------------------------------------------------------
con <- makeContrasts(grouptreated - groupcontrol, levels=design)
DE = glmQLFTest(fit, contrast=con)

de <- decideTestsDGE(DE, p=0.05, adjust="BH")
summary(de)

#' 
#' > EXTRA: Understanding the intercept model in more detail: 
#' 
## ------------------------------------------------------------------------
head(fit$coefficients)

#' 
#' We model the count RATE
#' 
#' average predicted count rates are given by :
#' 
#' log(E(Y/efflib)) = b0 + b1X (X=1 when group=treated)
#' for geneID00001 : log(E(Y/efflib)) = -12.1 -0.04X
#' log(E(Y/efflib|control)) = -12.1
#' log(E(Y/efflib|treatment)) = -12.1 -0.04
#' 
#' From this model, fold changes can be deduced:
#' 
#' logFC (treatment vs control)
#' = log(E(Y|treatment)) - log(E(Y|control))
#' = -0.04
#' thus FC = exp(-0.04) = 0.96
#' and log2FC = log2(0.96) = -0.059
#' 
#' However, fit$fitted.values does not return average predicted count rates but predicted counts for the given libsize :
#' 
#' log(E(Y/efflib)) = b0 + b1X
#' => log(E(Y)) - log(effective libsize) = -12.1 -0.04X
#' => log(E(Y)) = -12.1 -0.04X + log(effective libsize)
#' => E(Y) = exp(-12.1 -0.04X + log(effective libsize))
#' 
#' Verify with code for the first gene
#' 
## ------------------------------------------------------------------------
Ecounts.1 = exp(fit$coefficients[1,1] + log(lib.size.eff[1]))
Ecounts.2 = exp(fit$coefficients[1,1] + log(lib.size.eff[2]))
Ecounts.3 = exp(fit$coefficients[1,1] + log(lib.size.eff[3]))
Ecounts.4 = exp(fit$coefficients[1,1] + fit$coefficients[1,2] + log(lib.size.eff[4]))
Ecounts.5 = exp(fit$coefficients[1,1] + fit$coefficients[1,2] + log(lib.size.eff[5]))
Ecounts.6 = exp(fit$coefficients[1,1] + fit$coefficients[1,2] + log(lib.size.eff[6]))

c(Ecounts.1, Ecounts.2, Ecounts.3, Ecounts.4, Ecounts.5, Ecounts.6)

#' 
#' compare with
#' 
## ------------------------------------------------------------------------
fit$fitted.values[1,]

#' 
#' and the predicted count rates for a libsize of 1 million
#' 
## ------------------------------------------------------------------------
Erate.control = exp(fit$coefficients[1,1])*1000000
Erate.treated = exp(fit$coefficients[1,1] + fit$coefficients[1,2])*1000000

c(Erate.control, Erate.treated)

#' 
#' This is approximately the same as the average cpm for each treatment
#' 
## ------------------------------------------------------------------------
cpm(dfr, normalized.lib.sizes=TRUE, log=FALSE)[1,]
mean(cpm(dfr)[1,1:3])
mean(cpm(dfr)[1,4:6])

#' 
#' calculating fold changes treatment/control
#' 
## ------------------------------------------------------------------------
FC = exp(fit$coefficients[,2])
FC[1:20]

