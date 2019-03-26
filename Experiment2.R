### DESeq2 - December 2018 


#performs read counting and differential gene expression analysis
#plotting/visualisation of the data
#scatter plots
#PCA of samples
#MA-plot
#Volcano plot
#reporting significant results in excel format etc.

.libPaths()

source("https://bioconductor.org/biocLite.R")

library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("backports")
library("DESeq2")
library("ggplot2")
library("org.Hs.eg.db")

#BiocInstaller::biocLite(c("GenomicFeatures", "AnnotationDbi"))

#Loading the count data into R 

CountsTable_Experiment2 <- as.matrix(read.table("Experiment2_count_data.txt", header = TRUE))

#View the top of the "CountsTable_PPP2R3B" matrix object using the head command in R:

head(CountsTable_Experiment2)

#Here we are creating a new R matrix object called DESeq2Experiment2 using the as.matrix() function on our CountsTable_PPP2R3B object:

DESeq2Experiment2 <- as.matrix(CountsTable_Experiment2)

head(DESeq2Experiment2)

#-------------------------------------------------
######### Building our DESeq2 Experiemnt #########
#-------------------------------------------------

#lebel our data columns with untreated or treated status
condition <- factor(c(rep("untreated", 3), rep("treated", 3), rep("untreated", 3), rep("treated", 3)))

#likewise label our two cell lines
cellline <- factor(c(rep("SKMEL2", 6), rep("SKMEL30", 6)))

#now we create a dataframe with our 'condition' and cell line information
coldata <- data.frame(row.names=colnames(DESeq2Experiment2), condition, cellline)

#Next we will define our DESeq2 experiment design to the "condition" and "celline" of the count data of the cell lines as defined above. 
#This is acheived with the DESeq2 function DESeqDataSetFromMatrix()

DESeq2Experiment2 <- DESeqDataSetFromMatrix(countData=DESeq2Experiment2, colData=coldata, design = ~ cellline + condition)

#here we need to check our base level
DESeq2Experiment2$condition

#and relevel to ensure the untreated is base level
DESeq2Experiment2$condition <- relevel(DESeq2Experiment2$condition, "untreated")

#and juct checking the level again
DESeq2Experiment2$condition


#-------------------------------------------------
########### QC and Data Visualization ############
#-------------------------------------------------


rld <- rlog(DESeq2Experiment2, blind = FALSE)

################### SCATTER PLORS ################

# Lets have a look at the scatter plots using the pairs function. 
	# We will only show 10 representative graphs of randomly selected samples

# Define a function to draw a scatter plot for a pair of variables (samples) with density colours
plotFun <- function(x,y) { 

  dns <- densCols(x,y); 
  points(x,y, 
		col=dns, 
		pch=".", 
		panel.first=grid());  
		
  # abline(a=0, b=1, col="brown")

  }

# Plot the scatter plot for a few pairs of samples selected at random
set.seed(123) # forces the random number generator to produce fixed results

pairs(log2(counts(DESeq2Experiment2)[,sample(ncol(counts(DESeq2Experiment2)),5)] + 1), 
      panel=plotFun, 
	  lower.panel = NULL)


set.seed(123) # forces the random number generator to produce fixed results

pairs(log2(counts(DESeq2Experiment2)[,sample(ncol(counts(DESeq2Experiment2)),5)] + 1), 
      panel=plotFun, 
	  lower.panel = NULL)


####################### PCA ######################

#here we are using the plotPCA function of DESeq2 to visialise the principle components

plotPCA(rld, intgroup=c("condition", "cellline"))

#-------------------------------------------------
#### Differential gene expression analysis ######
#-------------------------------------------------

# DESeq2 analysis 
DESeq2Experiment2 <- DESeq(DESeq2Experiment2)

# Build a results table
	# Note that by default, the FDR cutoff of results() is 0.1
res2 <- results(DESeq2Experiment2)

# Get metadata of results table (res)
mcols(res2, use.names = TRUE)

# Get more summary of res
summary(res2)

# Filter results (res) by setting a lower cutoff
	# FRD (padj) - set it to < 0.05 to see how many significant results you have
res.05 <- results(DESeq2Experiment2,
					alpha = 0.05)

print("padj < 0.05")

table(res.05$padj < 0.05)



##################### PLOT-MA ####################

plotMA(res2, 
		ylim=c(-5,5))

#-------------------------------------------------
#################### Reporting ###################
#-------------------------------------------------

# Annotating and exporting results
	# Annotate with gene IDs (gene symbol and EntrezID); #org.Hs.eg.db
res2$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res2),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res2$entrez <- mapIds(org.Hs.eg.db,
                   keys=row.names(res2),
                   column="ENTREZID",
                   keytype="ENSEMBL",
                   multiVals="first")

#Ordering results by adjusted pvalues
resOrdered <- res2[order(res2$padj), ]
head(resOrdered)

##STEP13## Writing results to files 

# Exporting expression results
resOrderedDF <- as.data.frame(resOrdered)

output.CSV <- "Differential_expression_results_experiemnt2.csv"

write.csv(resOrderedDF, file = output.CSV)

################### PLOT-VOLCANO #################

alpha <- 0.05 # Threshold on the p-value

# par(mfrow=c(1,2))

# Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
res2$sig <- -log10(res2$padj)

sum(is.infinite(res$sig))

res2[is.infinite(res2$sig),"sig"] <- 350
# View(res[is.na(res$pvalue),])

# Select genes with a defined p-value (DESeq2 assigns NA to some genes)
genes.to.plot <- !is.na(res2$pvalue)
# sum(genes.to.plot)

range(res2[genes.to.plot,
		"log2FoldChange"])

# Volcano plot of adjusted p-values

cols <- densCols(res2$log2FoldChange, res2$sig)
cols[res2$pvalue ==0] <- "purple"
res2$pch <- 19
res2$pch[res$pvalue ==0] <- 6
plot(res2$log2FoldChange, 
     res2$sig, 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     pch=res2$pch, cex=0.4)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

## Plot the names of a reasonable number of genes, by selecting those begin not only significant but also having a strong effect size
gn.selected <- abs(res2$log2FoldChange) > 1 & res2$padj < alpha 
text(res2$log2FoldChange[gn.selected],
     -log10(res2$padj)[gn.selected],
     lab=res2$symbol[gn.selected ], cex=0.4)
	 

#-------------------------------------------------
#################### Finished ####################
#-------------------------------------------------
