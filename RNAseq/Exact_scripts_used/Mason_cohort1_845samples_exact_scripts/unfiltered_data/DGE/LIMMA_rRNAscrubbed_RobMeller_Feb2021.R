######################################################################################################
## Multiple Linear Regression of COVID expression in rRNA scrubbed data
## use Limma: dataformat-genes on the rows and samples on columns
## RMeller 09-Dec-2020

rm(list=ls());gc()
# load packages

library(limma)
library(edgeR)
library(lattice)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Glimma)
library(ggfortify)
library(broom) 
library(heatmap3)
library(tsne)
library(EnhancedVolcano)
setwd("D:/COVID/COV-IRT_Human_LIMMA"); getwd()
dir.create("PLOTS-ALL"); dir.create("RESULTS")



#############################################################################
################# 1. Load data and check a few factors ######################
#############################################################################
# Skip if you have already ran Section One

# Data files should be in this container
DATA="D:/COVID/NASA_DATA/Ribo_Cleaned_NASA Human_Reads"
#Read in the COVID gene expression data set 
Data <- read.csv(file.path(DATA,"STAR_Unnormalized_Counts_rRNA_removed.csv"), header = TRUE, check.names = FALSE)
str(Data); dim(Data)

# Make the first column rownames and transpose the table.
# This makes the genes on the columns and IDs on the rows. 
Data0 <- data.frame(t(Data[,-1]))
colnames (Data0) <- Data[,1]
rownames(Data0) <- str_remove(rownames(Data0), "_all-reads")
dim(Data0); Data0[1:20,1:20]; str(Data0)
#now add the rownames to a column called ID
Data0$ID <- rownames(Data0)
dim(Data0)

#need to remove rRNA reads
#rRNA <- unique(read.table(file.path(DATA,"rRNA_geneIDs.txt"), sep ="\t", check.names = FALSE))
# then match ENSG IDs ... using "-match"
#geneIDs <- colnames(Data0)
#datarows = match(rRNA[,1],geneIDs)
#clean_Data <- Data0[,-datarows]
#dim(clean_Data)

# just incase one is curious....
#ribodata <- Data0[, datarows,]
#write.table(ribodata, file = "RESULTS/rRNA_reads.txt", sep = "\t"); dim(ribodata)

## Load clinical trait data
traitData = read.csv(file.path(DATA,"MASON_Metadata_Combined.txt"), sep ="\t", header = TRUE);
dim(traitData); head(traitData)

# make the names shorter for nicer graphs
colnames(traitData)[1] <- ("ID")
colnames(traitData)[24] <- ("H.Sapiens_Kraken_Reads")
colnames(traitData)[25] <- ("SARS_CoV2_Kraken_Reads")
colnames(traitData)[28] <- ("H.Sapiens_Reads")
colnames(traitData)[27] <- ("SARS_CoV2_Reads")
rownames(traitData) <- traitData[,1]

# Now match the Samples and expression Data
Combined_data0 <- inner_join(traitData,Data0, by = "ID")
dim(Combined_data0)
# 732 Samples

# Now remove data that are duplicates
Combined_data <- Combined_data0[ Combined_data0$DeDupList != "FALSE", , drop=FALSE]
dim(Combined_data)
# 670 samples

# Save it just in case
write.table(Combined_data, file = "RESULTS/Combined_data.txt", sep = "\t")

# How many transcriptome aligned reads/ sample?
plotdata <- sort(as.numeric(rowSums(Combined_data[,32:60166]))*1e-6)
tiff(file="PLOTS-ALL/Counts_barplot.tiff", unit= "in", res = 300, width = 6, height = 5)
barplot(plotdata, ylim=c(0,50), xaxt = "n", xlab= "Samples", ylab = "Million counts", main="Distribution of Sample Reads")
abline(h = 5, col = "red"); dev.off()

plotdata2 <- sort(as.numeric(Combined_data$"H.Sapiens_Reads")*1e-6)
tiff(file="PLOTS-ALL/Counts_AlignedtoGenome_barplot.tiff", unit= "in", res = 300, width = 6, height = 5)
barplot(plotdata2, ylim=c(0,500), xaxt = "n", xlab= "Samples", ylab = "Million counts", main="Distribution of Sample Reads")
abline(h = 5, col = "red"); dev.off()

tiff(file="PLOTS-ALL/Counts_AlignedtoGenome_vs_TranscriptomeCounts.tiff", unit= "in", res = 300, width = 6, height = 5)
plot (log10(Combined_data$"H.Sapiens_Reads"), log10(rowSums(Combined_data[,32:60166])), pch=20)
abline(h = 5, col = "red"); dev.off()

# a few statistsics
cor.test(Combined_data$"H.Sapiens_Reads", rowSums(Combined_data[,32:60166]), method="pearson")
cor.test(Combined_data$"H.Sapiens_Reads", rowSums(Combined_data[,32:60166]), method="spearman")


## test with PCA as well, this works best on normalized data such as cpm data
# raw DATA
data0_PCA <- Data0[match(Combined_data$ID,rownames(Data0)),]
data0_PCA$ID <- NULL
data0_PCA <- sapply(data0_PCA[,], as.numeric)
Counts_as_cpm_Data0 <- cpm(t(data0_PCA), lib.size= Combined_data$H.Sapiens_Reads, log=TRUE, prior.count=0.25)
PCA1<- prcomp(Counts_as_cpm_Data0)
tiff(file="PLOTS-ALL/Counts_PCA.tiff", unit= "in", res = 300, width = 5, height = 5)
autoplot(PCA1); dev.off()

## Clean data to remove samples with low numbers of Transcriptome_reads
A <- rowSums(Combined_data[,32:60166])
keep <- A > 2000000
Data_2Mill <- Combined_data[keep,]
dim(Data_2Mill)

cpm_Data_2Mill <- cpm(t(Data_2Mill[,32:60166]), lib.size= Data_2Mill$H.Sapiens_Reads, log=TRUE, prior.count=0.25)
PCA3<- prcomp(cpm_Data_2Mill)
tiff(file="PLOTS-ALL/Counts_PCA_clean_2Mill.tiff", unit= "in", res = 300, width = 5, height = 5)
autoplot(PCA3, colour ="red"); dev.off()




## Save Datsheet for analysis
write.table(Data_2Mill, file = "RESULTS/Clean_Com_Data.txt", sep = "\t")


## how do alignments determine number of counted genes???
Data1 <- data.frame(t(Data_2Mill[,32:60166]))
colnames(Data1)<- Data_2Mill$ID
gene_counts <- apply(Data1, 2, function (x) sum(x>10))

tiff(file="PLOTS-ALL/Genescounted_vs_TranscriptomeCounts.tiff", unit= "in", res = 300, width = 6, height = 5)
plot (log10(rowSums(Data_2Mill[,32:60166])), (gene_counts), pch=20)
dev.off()
tiff(file="PLOTS-ALL/Genescounted_vs_GenomeCounts.tiff", unit= "in", res = 300, width = 6, height = 5)
plot (log10(Data_2Mill$"H.Sapiens_Reads"),(gene_counts), pch=20, col ="red")
dev.off()


## Save data files in RESULTS folder

# Save the Trait information
write.table(Data_2Mill[,1:32], "RESULTS/traitdata-LIMMA.txt", sep = "\t")

# Create the counts matrix 
#  For LIMMA data needs adapting so that genes are on the rows and sample on columns
# remove SARS data to new file (last 12 genes with names ENSSAS...
Data2 <- Data1[1:60123, ]
SARS <- Data1[60124:60135, ]
dim(Data2); write.table(as.matrix(Data2), "RESULTS/Counts-LIMMA.txt", sep = "\t")
dim(SARS); write.table(as.matrix(SARS), "RESULTS/Counts-LIMMA_SARS.txt", sep = "\t")


# now clean up
rm(list=ls()); gc()

###################################################################################
## 2. Now lets start to Limma#########################################################
###################################################################################
## 2.1. Limma requires a count table and trait data.

counts <- read.table("RESULTS/counts-LIMMA.txt", header=TRUE, sep = "\t")
traitData <- read.table(file = "RESULTS/traitData-LIMMA.txt", header=TRUE, sep = "\t") 

## 2.2 Now create an annotation file 

DATA="D:/COVID/NASA_DATA/ANALYSIS_DATA"
annot0 = read.table(file.path(DATA,"HumanAnnotation.txt"), sep = "\t", header = TRUE)
head(annot0)
# split up the annotation chromosome column
temp <- annot0 %>% separate(Chromosome, c("chr", "pos"), ":")
annot0 <- temp %>% separate(pos, c("start", "stop"), "-")
dim(annot0); head (annot0)

# now match to the Gene names in counts
samples = rownames(counts)
rowdata <- match(samples, annot0$gene_id)
annot<- annot0[rowdata,]
dim(annot); dim (counts)



# save the annotation information
write.table(annot, file = "RESULTS/annot-LIMMA.txt", sep = "\t")


## 2.3 Now create the design matrix and save it
head(traitData)
design <- model.matrix(~ 0 + Pred_sex + SequencingBatch + PCR_SEQ, data = traitData)
design <- model.matrix(~ 0 + SequencingBatch + PCR_SEQ, data = traitData)
write.table(design, file = "RESULTS/design-LIMMA.txt", sep = "\t")


## 2.4 Now create a contrasts matrix and save it
head(design)
contr.matrix <- makeContrasts(HIGHvsNONE = PCR_SEQHigh_Pos - PCR_SEQNone_Neg, 
	HIGHvsNONE_Pos = PCR_SEQHigh_Pos - PCR_SEQNone_Pos,
	MEDvsNONE = PCR_SEQMed_Pos - PCR_SEQNone_Neg, 
	LOWvsNONE = PCR_SEQLow_Pos - PCR_SEQNone_Neg,
	ViralvsNONE = PCR_SEQViral_Neg - PCR_SEQNone_Neg, 
	NONE_PosvsNeg = PCR_SEQNone_Pos - PCR_SEQNone_Neg, 
	MEDvsNONE_Pos = PCR_SEQMed_Pos - PCR_SEQNone_Pos, 
	LOWvsNONE_Pos = PCR_SEQLow_Pos - PCR_SEQNone_Pos,levels = colnames(design))
contr.matrix
write.table(contr.matrix, file = "RESULTS/contrast-LIMMA.txt", sep = "\t")


## 2.5 Now assemble the DGE object and assess
dge <- DGEList(counts=counts, gene=annot)
# filter genes by expression 566 samples, if >10 rum = 5000
A <- rowSums(dge$counts); isexpr <- A > 5000
dge <- dge[isexpr, , keep.lib.size = TRUE]


#apply TMM normalization
dge <- calcNormFactors(dge)

# create list of expressed genes
exp_gene_list <- dge$genes$GeneSymbol
write.table(exp_gene_list, file = "RESULTS/Exp_Genes_Filt.txt", sep = "\t")


## 2.5.1 analyse date to see if we have any strong covariates....
#tiff(file="PLOTS-ALL/MDS_10_pairwise.tiff", unit= "in", res = 300, width = 8, height = 5)
#plotMDS(dge, labels=Load, top=10, cex = 1, col=ifelse(Load=="m","blue","red"), gene.selection="pairwise", prior.count = 3
#dev.off()
#tiff(file="PLOTS-ALL/MDS_4_common.tiff", unit= "in", res = 300, width = 8, height = 5)
#plotMDS(dge, labels=Load, top=2, cex = 1, col=ifelse(Load=="m","blue","red"), gene.selection="pairwise", prior.count = 3)
#dev.off()
#So correction for sex is appropriate?  Used to show sex effect now gone.  


## 2.6 Now create our voom model for testing DGE(and get a plot) voon because library size > 3 fold diff
v <- voom(dge, design, plot = TRUE)
# then fit to the linear model using the design  
fit <-lmFit(v, design)
# fit to contrast matix, to identify contrast variables
cfit <- contrasts.fit(fit, contrasts=contr.matrix)
# finally apply Bayesian correction
efit <- eBayes(cfit)
# Mason paper used 1.5 fold changes 
tfit <- treat(efit, lfc=(log2(1.5)))
dt <- decideTests(tfit)
summary(dt)

# lets print this to a file
write.table((summary(dt)), file = "RESULTS/Down_Up_DGE_LIMMA.txt", sep = "\t")

# 2.7  MA plot of results
tiff(file="PLOTS-ALL/MD_HIGH-vs-Neg.tiff", unit= "in", res = 300, width = 8, height = 5)
plotMD(efit, column=1, status=dt[,1], main=colnames(efit)[1],legend="bottomright")
abline(h=0,col="darkgrey"); dev.off()

tiff(file="PLOTS-ALL/MD_HIGH-vs-NONE_Pos.tiff", unit= "in", res = 300, width = 8, height = 5)
plotMD(efit, column=2, status=dt[,2], main=colnames(efit)[2],legend="bottomright")
abline(h=0,col="darkgrey"); dev.off()

tiff(file="PLOTS-ALL/MD_MED-vs-Neg.tiff", unit= "in", res = 300, width = 8, height = 5)
plotMD(efit, column=3, status=dt[,3], main=colnames(efit)[3],legend="bottomright")
abline(h=0,col="darkgrey"); dev.off()

tiff(file="PLOTS-ALL/MD_LOW-vs-Neg.tiff", unit= "in", res = 300, width = 8, height = 5)
plotMD(efit, column=4, status=dt[,4], main=colnames(efit)[4],legend="bottomright")
abline(h=0,col="darkgrey"); dev.off()

tiff(file="PLOTS-ALL/MD_Virus-vs-Neg.tiff", unit= "in", res = 300, width = 8, height = 5)
plotMD(efit, column=5, status=dt[,5], main=colnames(efit)[5],legend="bottomright")
abline(h=0,col="darkgrey"); dev.off()

tiff(file="PLOTS-ALL/MD_NONE_Pos-vs-Neg.tiff", unit= "in", res = 300, width = 8, height = 5)
plotMD(efit, column=6, status=dt[,6], main=colnames(efit)[6],legend="bottomright")
abline(h=0,col="darkgrey"); dev.off()

tiff(file="PLOTS-ALL/MD_Med-vs-NONE_Pos.tiff", unit= "in", res = 300, width = 8, height = 5)
plotMD(efit, column=7, status=dt[,7], main=colnames(efit)[7],legend="bottomright")
abline(h=0,col="darkgrey"); dev.off()

tiff(file="PLOTS-ALL/MD_Low-vs-NONE_Pos.tiff", unit= "in", res = 300, width = 8, height = 5)
plotMD(efit, column=8, status=dt[,8], main=colnames(efit)[8],legend="bottomright")
abline(h=0,col="darkgrey"); dev.off()


# 2.5 Now compare between conditions (pairwise)
contr.matrix
de.common1 <- which(dt[,1]!=0 & dt[,2]!=0)
head(efit$genes$GeneSymbol[de.common1])
#compare viral vs High covid
de.common2 <- which(dt[,1]!=0 & dt[,3]!=0)
head(efit$genes$GeneSymbol[de.common2])
#high pos or none pos vs none_neg
de.common3 <- which(dt[,1]!=0 & dt[,4]!=0)
head(efit$genes$GeneSymbol[de.common3])

# Lists are boring, so lets look at the data in a VennDiagram
tiff(file="PLOTS-ALL/VennGenes_HIGHvsNeg_&_Pos.tiff", unit= "in", res = 300, width = 3, height = 3)
vennDiagram(dt[,1:2], circle.col=c("blue", "salmon"), mar=rep(1,4), names=c("High vs Neg", "High vs Pos"), cex=c(1.0, .5, 0.25), main="COVID(High)vs NONE")
dev.off()

tiff(file="PLOTS-ALL/VennGenes_COVID-high-med-low.tiff", unit= "in", res = 300, width = 3, height = 3)
vennDiagram(dt[,c(1,3,4)], circle.col=c("blue", "yellow", "pink"), names=c("High", "Med", "Low"), cex=c(1.0, .5, 0.25), main="COVID Viral Load")
dev.off()

tiff(file="PLOTS-ALL/VennGenes_COVIDvsVIRAL.tiff", unit= "in", res = 300, width = 3, height = 3)
vennDiagram(dt[,c(1,5)], circle.col=c("blue", "red"), names=c("COVID", "Other Viral"), cex=c(1.0, .5, 0.25), main="COVID vs virus")
dev.off()

tiff(file="PLOTS-ALL/VennGenes_Seq vs PCR.tiff", unit= "in", res = 300, width = 3, height = 3)
vennDiagram(dt[,c(1,6)], circle.col=c("blue", "green"),names=c("High", "None-SeqPos"), cex=c(1.0, .5, 0.25), main="Seq vs PCR")
dev.off()
tiff(file="PLOTS-ALL/VennGenes_Med vs None_POS.tiff", unit= "in", res = 300, width = 3, height = 3)
vennDiagram(dt[,c(3,6)], circle.col=c("salmon", "green"),names=c("Med", "None-SeqPos"), cex=c(1.0, .5, 0.25), main="Seq vs PCR")
dev.off()

tiff(file="PLOTS-ALL/VennGenes_Low vs None_POS.tiff", unit= "in", res = 300, width = 3, height = 3)
vennDiagram(dt[,c(4,6)], circle.col=c("pink", "green"),names=c("Low", "None-SeqPos"), cex=c(1.0, .5, 0.25), main="Seq vs PCR")
dev.off()

summary(dt)
# print table of results (p<0.05)
Table1 <- topTable(efit, coef = 1, n =Inf, sort = "p", p = 0.050)[,-1]
head(Table1); write.table(Table1, file = "RESULTS/SIG_Genes_Table1.txt", sep = "\t"); dim(Table1)
Table2 <- topTable(efit, coef = 2, n =Inf, sort = "p", p = 0.050)[,-1]
head(Table2); write.table(Table2, file = "RESULTS/SIG_Genes_Table2.txt", sep = "\t"); dim(Table2)
Table3 <- topTable(efit, coef = 3, n =Inf, sort = "p", p = 0.050)[,-1]
head(Table3); write.table(Table3, file = "RESULTS/SIG_Genes_Table3.txt", sep = "\t")
Table4 <- topTable(efit, coef = 4, n =Inf, sort = "p", p = 0.050)[,-1]
head(Table4); write.table(Table4, file = "RESULTS/SIG_Genes_Table4.txt", sep = "\t")
Table5 <- topTable(efit, coef = 5, n =Inf, sort = "p", p = 0.050)[,-1]
head(Table5); write.table(Table5, file = "RESULTS/SIG_Genes_Table5.txt", sep = "\t")
Table6 <- topTable(efit, coef = 6, n =Inf, sort = "p", p = 0.050)[,-1]
head(Table6); write.table(Table6, file = "RESULTS/SIG_Genes_Table6.txt", sep = "\t")
Table7 <- topTable(efit, coef = 7, n =Inf, sort = "p", p = 0.050)[,-1]
head(Table7); write.table(Table7, file = "RESULTS/SIG_Genes_Table7.txt", sep = "\t")
Table8 <- topTable(efit, coef = 8, n =Inf, sort = "p", p = 0.050)[,-1]
head(Table8); write.table(Table8, file = "RESULTS/SIG_Genes_Table8.txt", sep = "\t")

dir.create("RESULTS/AB")
## print table of results for Afshin
Table1A <- topTable(efit, coef = 1, n =Inf, sort = "p", p = 1.0)[,-1]
head(Table1A); write.table(Table1A, file = "RESULTS/AB/SIG_Genes_Table1.txt", sep = "\t")
Table2A <- topTable(efit, coef = 2, n =Inf, sort = "p", p = 1.0)[,-1]
head(Table2A); write.table(Table2A, file = "RESULTS/AB/SIG_Genes_Table2.txt", sep = "\t")
Table3A <- topTable(efit, coef = 3, n =Inf, sort = "p", p = 1.0)[,-1]
head(Table3A); write.table(Table3A, file = "RESULTS/AB/SIG_Genes_Table3.txt", sep = "\t")
Table4A <- topTable(efit, coef = 4, n =Inf, sort = "p", p = 1.0)[,-1]
head(Table4A); write.table(Table4A, file = "RESULTS/AB/SIG_Genes_Table4.txt", sep = "\t")
Table5A <- topTable(efit, coef = 5, n =Inf, sort = "p", p = 1.0)[,-1]
head(Table5A); write.table(Table5A, file = "RESULTS/AB/SIG_Genes_Table5.txt", sep = "\t")
Table6A <- topTable(efit, coef = 6, n =Inf, sort = "p", p = 1.0)[,-1]
head(Table6A); write.table(Table6A, file = "RESULTS/AB/SIG_Genes_Table6.txt", sep = "\t")
Table7A <- topTable(efit, coef = 7, n =Inf, sort = "p", p = 1.0)[,-1]
head(Table7A); write.table(Table7A, file = "RESULTS/AB/SIG_Genes_Table7.txt", sep = "\t")
Table8A <- topTable(efit, coef = 8, n =Inf, sort = "p", p = 1.0)[,-1]
head(Table8A); write.table(Table8A, file = "RESULTS/AB/SIG_Genes_Table8.txt", sep = "\t")


## 2.6. Volcanoplots..../ \......

# if starting fresh
#Table1 <- read.table("SIG_Genes_Table1.txt", sep = "\t", header = TRUE)
# Create dataframe of log2FC and adj pvalues....
str(Table1)

res1 <- data.frame((Table1[,c(7,11)]));
rownames(res1) <- make.names(Table1$GeneSymbol, unique=TRUE) 
colnames(res1)<- c("log2FC","adj.P.Value")
# now draw the figure.....
tiff(file="PLOTS-ALL/Volcano-High COVID vs Neg.tiff", unit= "in", res = 300, width = 5, height = 5)
EnhancedVolcano(res1,lab = rownames(res1),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "HIGH COVID vs Neg"); dev.off()

# and repeat...
tiff(file="PLOTS-ALL/Volcano-High COVID vs NONE_POS.tiff", unit= "in", res = 300, width = 5, height = 5)
res2 <- data.frame(Table2[,c(7,11)]);
rownames(res2) <- make.names(Table2$GeneSymbol, unique=TRUE) 
colnames(res2)<- c("log2FC","adj.P.Value")
# now draw the figure.....
EnhancedVolcano(res2,lab = rownames(res2),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "HIGH COVID vs NONE_Pos");dev.off()

tiff(file="PLOTS-ALL/Volcano-Med COVID vs Neg.tiff", unit= "in", res = 300, width = 5, height = 5)
res3 <- data.frame(Table3[,c(7,11)]);
rownames(res3) <- make.names(Table3$GeneSymbol, unique=TRUE) 
colnames(res3)<- c("log2FC","adj.P.Value")
# now draw the figure.....
EnhancedVolcano(res3,lab = rownames(res3),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "MED COVID vs Neg");dev.off()

tiff(file="PLOTS-ALL/Volcano-Low COVID vs Neg.tiff", unit= "in", res = 300, width = 5, height = 5)
res4 <- data.frame(Table4[,c(7,11)]);
rownames(res4) <- make.names(Table4$GeneSymbol, unique=TRUE) 
colnames(res4)<- c("log2FC","adj.P.Value")
# now draw the figure.....
EnhancedVolcano(res4,lab = rownames(res4),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "LOW COVID vs Neg");dev.off()

tiff(file="PLOTS-ALL/Volcano-Virus vs Neg.tiff", unit= "in", res = 300, width = 5, height = 5)
res5 <- data.frame(Table5[,c(7,11)]);
rownames(res5) <- make.names(Table5$GeneSymbol, unique=TRUE) 
colnames(res5)<- c("log2FC","adj.P.Value")
# now draw the figure.....
EnhancedVolcano(res5,lab = rownames(res5),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "Virus vs Neg");dev.off()

tiff(file="PLOTS-ALL/Volcano- NONE_Pos vs Neg.tiff", unit= "in", res = 300, width = 5, height = 5)
res6 <- data.frame(Table6[,c(7,11)]);
rownames(res6) <- make.names(Table6$GeneSymbol, unique=TRUE) 
colnames(res6)<- c("log2FC","adj.P.Value")
# now draw the figure.....
EnhancedVolcano(res6,lab = rownames(res6),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "NONE_Pos vs Neg");dev.off()

tiff(file="PLOTS-ALL/Volcano-Med-v-NONE_POS.tiff", unit= "in", res = 300, width = 5, height = 5)
res7 <- data.frame(Table7[,c(7,11)]);
rownames(res7) <- make.names(Table7$GeneSymbol, unique=TRUE) 
colnames(res7)<- c("log2FC","adj.P.Value")
# now draw the figure.....
EnhancedVolcano(res7,lab = rownames(res7),
    x = 'log2FC',
    y = 'adj.P.Value',     
    xlim = c(-5, 8), title = "MED vs NONE_Pos");dev.off()

tiff(file="PLOTS-ALL/Volcano-Low-vs-NONE_POS.tiff", unit= "in", res = 300, width = 5, height = 5)
res8 <- data.frame(Table8[,c(7,11)]);
rownames(res8) <- make.names(Table8$GeneSymbol, unique=TRUE) 
colnames(res8)<- c("log2FC","adj.P.Value")
# now draw the figure.....
EnhancedVolcano(res8,lab = rownames(res8),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "LOW vs NONE_Pos");dev.off()


###########################################################################

###  3. Cluster analysis of results
## Use data file saved at end of section 1, convert to cpm and then filter by DEGs.  

Combined_data <- read.table(file = "RESULTS/Clean_Com_Data.txt", sep = "\t", header = TRUE)
dim(Combined_data); # !! 1-31 are trait data; 32-60166 are count data

# convert counts to CPM and rebuild table
Counts_as_cpm <- cpm(t(Combined_data[ ,32:60166]), lib.size= Combined_data$H.Sapiens_Reads, log=TRUE, prior.count=0.25)
Combined_data_cpm <- (as.data.frame(cbind(Combined_data[ ,1:31], t(Counts_as_cpm))))
write.table(Combined_data_cpm, file = "RESULTS/Clean_Com_Data_cpm.txt", sep = "\t")

# quick check of all data
dim(Combined_data_cpm)
tiff(file="PLOTS-ALL/PCA_All_Genes.tiff", unit= "in", res = 300, width = 8, height = 5)
prcomp_genes <- prcomp(Combined_data_cpm[,32:60166], scale = FALSE)
autoplot(prcomp_genes, data=Combined_data_cpm, colour='PCR_SEQ')
dev.off()

# Now filter by DEGs for High
#highDEG<-rownames(Table1[Table1$"adj.P.Val"  <= 0.05 & Table1$logFC  <= -1.48 | Table1$logFC  >= 1.48,, drop=FALSE])
# Use dt table
res_table <- data.frame(dt)
head (res_table)

# filter by col where value does not = 0
highDEG_sig<-rownames(res_table[res_table[,1] != 0,, drop=FALSE])
highDEG_sig

high_counts <- match(highDEG_sig, colnames(Combined_data_cpm))
Counts <- Combined_data_cpm[ ,high_counts]
dim(Counts)
Counts2 <- sapply(Counts[,1:325],as.numeric)
comb_DEGHigh_counts_trait <- cbind(Combined_data_cpm[ ,1:31], Counts2)

##now filter off high, low and none
comb_DEGHigh_counts_trait <- comb_DEGHigh_counts_trait[ comb_DEGHigh_counts_trait$PCR_SEQ == "High_Pos"|comb_DEGHigh_counts_trait$PCR_SEQ == "None_Neg"|comb_DEGHigh_counts_trait$PCR_SEQ == "None_Pos", , drop=FALSE]
dim(comb_DEGHigh_counts_trait)


## 3.1. SUPERVISED CLUSTERING_PCA
tiff(file="PLOTS-ALL/PCA_HighvsNeg.tiff", unit= "in", res = 300, width = 8, height = 5)
prcomp_genes <- prcomp(comb_DEGHigh_counts_trait[,32:356], scale = FALSE)
autoplot(prcomp_genes, data=comb_DEGHigh_counts_trait, colour="PCR_SEQ")+ scale_color_manual(values=c("red", "gray", "blue")); 
dev.off()


## 3.2. UNSUPERVISED CLUSTERING TSNE
# try tsne may need to split the data a bit too many factors
library(tsne)
set.seed=42
tsne <- tsne(comb_DEGHigh_counts_trait[,32:325], initial_config = NULL, k = 2, 
	initial_dims = 30, perplexity = 100, max_iter = 500, min_cost = 0, epoch_callback = NULL, whiten = TRUE,epoch=100)
cols=as.factor((comb_DEGHigh_counts_trait[, 30]))
tiff(file="PLOTS-ALL/TSNE_HighvsNeg.tiff", unit= "in", res = 300, width = 5, height = 5)
plot(tsne[,1], tsne[,2], pch =19, col=cols)
legend("bottomleft", legend = unique(cols), col=unique(cols), pch=19);dev.off()


## 3.3 HIERARCHICAL CLUSTERING
# create color vector
map_treatment_to_color <- function(comb_DEGHigh_counts_trait){
	colorsVector <- ifelse (comb_DEGHigh_counts_trait$PCR_SEQ == "High_Pos", "blue", 
	ifelse (comb_DEGHigh_counts_trait$PCR_SEQ == "None_Neg", "red", "yellow"))
	return(colorsVector)
	}
sample_colors = map_treatment_to_color(comb_DEGHigh_counts_trait)

tiff(file="PLOTS-ALL/HC_HighvsNeg.tiff", unit= "in", res = 300, width = 5, height = 5)
heatmap3(t(scale(comb_DEGHigh_counts_trait[ ,32:325])), ColSideColors=sample_colors)
dev.off()

## Now repeat for LOW
names(res_table)
LowDEG<-rownames(res_table[res_table[,4] != 0,, drop=FALSE])

low_counts <- match(LowDEG, colnames(Combined_data_cpm))
CountsB <- Combined_data_cpm[ ,low_counts]
dim(CountsB)
Counts2B <- sapply(CountsB[,1:57],as.numeric)
Comb_DEGLOW_counts_trait <- cbind(Combined_data_cpm[ ,1:31], Counts2B)
dim(Comb_DEGLOW_counts_trait)

##now filter off high, low and none
Comb_DEGLOW_counts_trait <- Comb_DEGLOW_counts_trait[ Comb_DEGLOW_counts_trait$PCR_SEQ == "Low_Pos"|Comb_DEGLOW_counts_trait$PCR_SEQ == "None_Neg"|Comb_DEGLOW_counts_trait$PCR_SEQ == "None_Pos", , drop=FALSE]
dim(Comb_DEGLOW_counts_trait)


## 3.4. SUPERVISED CLUSTERING_PCA  -LOW
prcomp_genes <- prcomp(Comb_DEGLOW_counts_trait[,32:88], scale = FALSE)
#autoplot(prcomp_genes, data=Comb_DEGLOW_counts_trait, colour='PCR_SEQ', frame.type = 't') + scale_color_manual(values=c("red", "gray", "blue"))
tiff(file="PLOTS-ALL/PCA_LowvsNeg.tiff", unit= "in", res = 300, width = 8, height = 5)
autoplot(prcomp_genes, data=Comb_DEGLOW_counts_trait, colour='PCR_SEQ') + scale_color_manual(values=c("red", "gray", "blue")); dev.off()

## 3.5. UNSUPERVISED CLUSTERING TSNE  -LOW
library(tsne)
set.seed=42
tsne <- tsne(Comb_DEGLOW_counts_trait[,32:88], initial_config = NULL, k = 2, initial_dims = 30, perplexity = 100, max_iter = 500, min_cost = 0, epoch_callback = NULL, whiten = TRUE,epoch=100)
cols=as.factor((Comb_DEGLOW_counts_trait[, 30]))
tiff(file="PLOTS-ALL/TSNE_LowvsNeg.tiff", unit= "in", res = 300, width = 5, height = 5)
plot(tsne[,1], tsne[,2], pch =19, col=cols)
legend("bottomleft", legend = unique(cols), col=unique(cols), pch=19); dev.off()



## 3.6 HIERARCHICAL CLUSTERING -LOW
map_treatment_to_color <- function(Comb_DEGLOW_counts_trait){
	colorsVector <- ifelse (Comb_DEGLOW_counts_trait$PCR_SEQ == "Low_Pos", "blue", 
	ifelse (Comb_DEGLOW_counts_trait$PCR_SEQ == "None_Neg", "red", "yellow"))
	return(colorsVector)
	}
sample_colors = map_treatment_to_color(comb_DEGHigh_counts_trait)

tiff(file="PLOTS-ALL/HC_LowvsNeg.tiff", unit= "in", res = 300, width = 5, height = 5)
heatmap3(t(scale(Comb_DEGLOW_counts_trait[ ,32:88])), ColSideColors=sample_colors); dev.off()



###################################################################################
############  4. Pipe Data to Gene Set Analysis    ################################
###################################################################################
# load Libraries
setwd
library(topGO)
library(GO.db)
library(org.Hs.eg.db)
library(Rgraphviz)
rm(list=ls())
##....if starting here
counts <- read.table("RESULTS/counts-LIMMA.txt", header=TRUE, sep = "\t"); head(counts)
annot <- read.table("RESULTS/annot-LIMMA.txt", header=TRUE, sep = "\t");head(annot)
design <- read.table("RESULTS/design-LIMMA.txt", header=TRUE, sep = "\t");head(design)
contr.matrix <- read.table("RESULTS/counts-LIMMA.txt", header=TRUE, sep = "\t");head(annot)
expressed_genes <-read.table("RESULTS/Exp_Genes_All.txt", header=TRUE, sep = "\t");head(expressed_genes)
Table1 <- read.table("RESULTS/SIG_Genes_Table1.txt", header=TRUE, sep = "\t"); head(Table1)
Table2 <- read.table("RESULTS/SIG_Genes_Table2.txt", header=TRUE, sep = "\t"); head(Table2)
Table3 <- read.table("RESULTS/SIG_Genes_Table3.txt", header=TRUE, sep = "\t"); head(Table3)
Table4 <- read.table("RESULTS/SIG_Genes_Table4.txt", header=TRUE, sep = "\t"); head(Table4)
Table5 <- read.table("RESULTS/SIG_Genes_Table5.txt", header=TRUE, sep = "\t"); head(Table5)
Table6 <- read.table("RESULTS/SIG_Genes_Table6.txt", header=TRUE, sep = "\t"); head(Table6)

## 3.1 Analysis using topGO
# 3.1.1 Table 1 analysis High vs None
all_genes <- as.vector(annot$GeneSymbol)
table1_genes <- as.vector(Table1$GeneSymbol)


# then make a factor that is 1 if the probeset is "interesting" and 0 otherwise
geneList <- factor(as.integer(all_genes %in% table1_genes))
# name the factor with the probeset names
names (geneList) <- all_genes
str(geneList)

GOdata_BP1 <-new ("topGOdata", 
    ontology = "BP", 
    allGenes = geneList, 
    nodeSize = 5, 
    annot = annFUN.org, 
    mapping = "org.Hs.eg.db", 
    ID = "symbol")

# Run Fisher test
resultFisher_BP <- runTest(GOdata_BP1, algorithm = "classic", statistic = "fisher")
resultKS_BP <- runTest(GOdata_BP1, algorithm = "classic", statistic = "ks")
resultKS.elim_BP <- runTest(GOdata_BP1, algorithm = "elim", statistic = "ks")
allRes_BP1 <- GenTable(GOdata_BP1, classicFisher = resultFisher_BP,
                   classicKS = resultKS_BP, elimKS = resultKS_BP,
                   orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 30)

# Collate results...
allRes_BP1
write.table(allRes_BP1, file = "RESULTS/GO_results_BP_Table1.txt", sep = "\t")

# Draw node image
tiff(file = "PLOTS-ALL/GO_BP_Nodes_Table1.tiff", res=500, unit = "in", width = 12, height = 9);
showSigOfNodes(GOdata_BP1, score(resultKS_BP), firstSigNodes = 15, useInfo = 'all')
dev.off()

## 3.1.2 Table 2 analysis HighPos vs NonePos
all_genes <- as.vector(annot$GeneSymbol)
table2_genes <- as.vector(Table2$GeneSymbol)

# then make a factor that is 1 if the probeset is "interesting" and 0 otherwise
geneList <- factor(as.integer(all_genes %in% table2_genes))
# name the factor with the probeset names
names (geneList) <- all_genes
str(geneList)

GOdata_BP2 <-new ("topGOdata", 
    ontology = "BP", 
    allGenes = geneList, 
    nodeSize = 5, 
    annot = annFUN.org, 
    mapping = "org.Hs.eg.db", 
    ID = "symbol")

# run Fisher test
resultFisher_BP <- runTest(GOdata_BP2, algorithm = "classic", statistic = "fisher")
resultKS_BP <- runTest(GOdata_BP2, algorithm = "classic", statistic = "ks")
resultKS.elim_BP <- runTest(GOdata_BP2, algorithm = "elim", statistic = "ks")
allRes_BP2 <- GenTable(GOdata_BP2, classicFisher = resultFisher_BP,
                   classicKS = resultKS_BP, elimKS = resultKS_BP,
                   orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 30)

# collate results...
allRes_BP2
write.table(allRes_BP2, file = "RESULTS/GO_results_BP_Table2.txt", sep = "\t")

## draw node image

tiff(file = "PLOTS-ALL/GO_BP_Nodes_Table2.tiff", res=500, unit = "in", width = 12, height = 9);
showSigOfNodes(GOdata_BP2, score(resultKS_BP), firstSigNodes = 15, useInfo = 'all')
dev.off()

# 3.1.3 Table 3 analysis Med vs None
all_genes <- as.vector(annot$GeneSymbol)
table3_genes <- as.vector(Table3$GeneSymbol)

# then make a factor that is 1 if the probeset is "interesting" and 0 otherwise
geneList <- factor(as.integer(all_genes %in% table3_genes))
# name the factor with the probeset names
names (geneList) <- all_genes
str(geneList)

GOdata_BP3 <-new ("topGOdata", 
    ontology = "BP", 
    allGenes = geneList, 
    nodeSize = 5, 
    annot = annFUN.org, 
    mapping = "org.Hs.eg.db", 
    ID = "symbol")

# run Fisher test
resultFisher_BP <- runTest(GOdata_BP3, algorithm = "classic", statistic = "fisher")
resultKS_BP <- runTest(GOdata_BP3, algorithm = "classic", statistic = "ks")
resultKS.elim_BP <- runTest(GOdata_BP3, algorithm = "elim", statistic = "ks")
allRes_BP3 <- GenTable(GOdata_BP3, classicFisher = resultFisher_BP,
                   classicKS = resultKS_BP, elimKS = resultKS_BP,
                   orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 30)

# collate results...
allRes_BP3
write.table(allRes_BP3, file = "RESULTS/GO_results_BP_Table3.txt", sep = "\t")

## draw node image

tiff(file = "PLOTS-ALL/GO_BP_Nodes_Table3.tiff", res=500, unit = "in", width = 12, height = 9);
showSigOfNodes(GOdata_BP3, score(resultKS_BP), firstSigNodes = 15, useInfo = 'all')
dev.off()

# 3.1.4 Table 4 analysis Low vs NoneNeg
all_genes <- as.vector(annot$GeneSymbol)
table4_genes <- as.vector(Table4$GeneSymbol)

# then make a factor that is 1 if the probeset is "interesting" and 0 otherwise
geneList <- factor(as.integer(all_genes %in% table4_genes))
# name the factor with the probeset names
names (geneList) <- all_genes
str(geneList)

GOdata_BP4 <-new ("topGOdata", 
    ontology = "BP", 
    allGenes = geneList, 
    nodeSize = 5, 
    annot = annFUN.org, 
    mapping = "org.Hs.eg.db", 
    ID = "symbol")

# run Fisher test
resultFisher_BP <- runTest(GOdata_BP4, algorithm = "classic", statistic = "fisher")
resultKS_BP <- runTest(GOdata_BP4, algorithm = "classic", statistic = "ks")
resultKS.elim_BP <- runTest(GOdata_BP4, algorithm = "elim", statistic = "ks")
allRes_BP4 <- GenTable(GOdata_BP4, classicFisher = resultFisher_BP,
                   classicKS = resultKS_BP, elimKS = resultKS_BP,
                   orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 30)

# collate results...
allRes_BP4
write.table(allRes_BP4, file = "RESULTS/GO_results_BP_Table4.txt", sep = "\t")

## draw node image

tiff(file = "PLOTS-ALL/GO_BP_Nodes_Table4.tiff", res=500, unit = "in", width = 12, height = 9);
showSigOfNodes(GOdata_BP4, score(resultKS_BP), firstSigNodes = 15, useInfo = 'all')
dev.off()


# 3.1.5 Table 5 analysis Virus vs NoneNeg
all_genes <- as.vector(annot$GeneSymbol)
table5_genes <- as.vector(Table5$GeneSymbol)

# then make a factor that is 1 if the probeset is "interesting" and 0 otherwise
geneList <- factor(as.integer(all_genes %in% table5_genes))
# name the factor with the probeset names
names (geneList) <- all_genes
str(geneList)

GOdata_BP5 <-new ("topGOdata", 
    ontology = "BP", 
    allGenes = geneList, 
    nodeSize = 5, 
    annot = annFUN.org, 
    mapping = "org.Hs.eg.db", 
    ID = "symbol")

# run Fisher test
resultFisher_BP <- runTest(GOdata_BP5, algorithm = "classic", statistic = "fisher")
resultKS_BP <- runTest(GOdata_BP5, algorithm = "classic", statistic = "ks")
resultKS.elim_BP <- runTest(GOdata_BP5, algorithm = "elim", statistic = "ks")
allRes_BP5 <- GenTable(GOdata_BP5, classicFisher = resultFisher_BP,
                   classicKS = resultKS_BP, elimKS = resultKS_BP,
                   orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 30)

# collate results...
allRes_BP5
write.table(allRes_BP5, file = "RESULTS/GO_results_BP_Table5.txt", sep = "\t")

## draw node image

tiff(file = "PLOTS-ALL/GO_BP_Nodes_Table5.tiff", res=500, unit = "in", width = 12, height = 9);
showSigOfNodes(GOdata_BP5, score(resultKS_BP), firstSigNodes = 15, useInfo = 'all')
dev.off()


# 3.1.6 Table 6 analysis NonePos vs NoneNeg
all_genes <- as.vector(annot$GeneSymbol)
table6_genes <- as.vector(Table6$GeneSymbol)

# then make a factor that is 1 if the probeset is "interesting" and 0 otherwise
geneList <- factor(as.integer(all_genes %in% table6_genes))
# name the factor with the probeset names
names (geneList) <- all_genes
str(geneList)

GOdata_BP6 <-new ("topGOdata", 
    ontology = "BP", 
    allGenes = geneList, 
    nodeSize = 5, 
    annot = annFUN.org, 
    mapping = "org.Hs.eg.db", 
    ID = "symbol")

# run Fisher test
resultFisher_BP <- runTest(GOdata_BP6, algorithm = "classic", statistic = "fisher")
resultKS_BP <- runTest(GOdata_BP6, algorithm = "classic", statistic = "ks")
resultKS.elim_BP <- runTest(GOdata_BP6, algorithm = "elim", statistic = "ks")
allRes_BP6 <- GenTable(GOdata_BP6, classicFisher = resultFisher_BP,
                   classicKS = resultKS_BP, elimKS = resultKS_BP,
                   orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 30)

# collate results...
allRes_BP6
write.table(allRes_BP6, file = "RESULTS/GO_results_BP_Table6.txt", sep = "\t")

## draw node image

tiff(file = "PLOTS-ALL/GO_BP_Nodes_Table6.tiff", res=500, unit = "in", width = 12, height = 9);
showSigOfNodes(GOdata_BP6, score(resultKS_BP), firstSigNodes = 15, useInfo = 'all')
dev.off()

## 3.2. Compare GO terms for DATA
#LIMMA vennDiagram function requires data as 1 and 0 values, to use list of values use VennDiagram package.  
library(VennDiagram)

## 3.2.1 compare HIGHvsNONE vs LOWvsNONE  
#Compare using x$GO.ID or first term

venn.diagram(list(High = allRes_BP1[,1], Low = allRes_BP4[,1]), filename = "Plots-ALL/HighvsLow-GOterms.tiff", height = 3000, width = 3000, resolution =
    300, imagetype = "tiff", units = "px", main = " High vs Low GO Terms", 
main.cex = 3, cat.fontface="bold", scaled = TRUE, label.cex = 3, cat.pos = c(330, 30), cat.cex = 2.5, fill = c("red", "lightblue"), alpha = 0.5, cex = 3)

venn.diagram(list(NonePOS = allRes_BP6[,1], Low = allRes_BP4[,1]), filename = "Plots-ALL/NonePOSvsLow-GOterms.tiff", height = 3000, width = 3000, resolution =
    300, imagetype = "tiff", units = "px", main = " NonePOS vs Low GO Terms", 
main.cex = 3, cat.fontface="bold", scaled = TRUE, label.cex = 3, cat.pos = c(330, 30), cat.cex = 2.5, fill = c("yellow", "lightblue"), alpha = 0.5, cex = 3)

venn.diagram(list(High = allRes_BP1[,1], Virus = allRes_BP5[,1]), filename = "Plots-ALL/HighvsVirus-GOterms.tiff", height = 3000, width = 3000, resolution =
    300, imagetype = "tiff", units = "px", main = " High vs Low GO Terms", 
main.cex = 3, cat.fontface="bold", scaled = TRUE, label.cex = 3, cat.pos = c(330, 30), cat.cex = 2.5, fill = c("red", "pink"), alpha = 0.5, cex = 3)

venn.diagram(list(High = allRes_BP1[,1], Med = allRes_BP3[,1]), filename = "Plots-ALL/HighvsMed-GOterms.tiff", height = 3000, width = 3000, resolution =
    300, imagetype = "tiff", units = "px", main = " High vs Low GO Terms", 
main.cex = 3, cat.fontface="bold", scaled = TRUE, label.cex = 3, cat.pos = c(330, 30), cat.cex = 2.5, fill = c("red", "green"), alpha = 0.5, cex = 3)

venn.diagram(list(Low = allRes_BP4[,1], Med = allRes_BP3[,1]), filename = "Plots-ALL/LowvsMed-GOterms.tiff", height = 3000, width = 3000, resolution =
    300, imagetype = "tiff", units = "px", main = " Med vs Low GO Terms", 
main.cex = 3, cat.fontface="bold", scaled = TRUE, label.cex = 3, cat.pos = c(330, 30), cat.cex = 2.5, fill = c("lightblue", "green"), alpha = 0.5, cex = 3)

venn.diagram(list(NonePOS = allRes_BP6[,1], Med = allRes_BP3[,1]), filename = "Plots-ALL/NonePosvsMed-GOterms.tiff", height = 3000, width = 3000, resolution =
    300, imagetype = "tiff", units = "px", main = " Med vs NonePOS GO Terms", 
main.cex = 3, cat.fontface="bold", scaled = TRUE, label.cex = 3, cat.pos = c(330, 30), cat.cex = 2.5, fill = c("yellow", "green"), alpha = 0.5, cex = 3)


## now compare all vSARS viral loads
venn.diagram(list(High = allRes_BP1[,1], Med = allRes_BP3[,1], Low = allRes_BP4[,1]), filename = "Plots-ALL/HighvsMedvsLow-GOterms.tiff", height = 3000, width = 3000, resolution =
    300, imagetype = "tiff", units = "px", main = " High vs Low GO Terms", 
main.cex = 3, cat.fontface="bold", scaled = TRUE, label.cex = 3, cat.pos = c(330, 30, 180), cat.cex = 2.5, fill = c("red", "green", "lightblue"), alpha = 0.5, cex = 3)

#########################################################################################
#########  5. Gene Set enrichment analysis with FGSEA ###################################
#########################################################################################
setwd("D:/COVID/COV-IRT_Human_LIMMA"); getwd()
dir.create ("GSEA_Results")

####GSEA
#library(RGSEA)
library(fgsea)
set.seed(123)
GMT_C1 ="D:/SEQUENCING/MsigDB/msigdb_v7.2_GMTs/c1.all.v7.2.symbols.gmt"
GMT_C2 ="D:/SEQUENCING/MsigDB/msigdb_v7.2_GMTs/c2.all.v7.2.symbols.gmt"
GMT_C3 ="D:/SEQUENCING/MsigDB/msigdb_v7.2_GMTs/c3.all.v7.2.symbols.gmt"
GMT_C5="D:/SEQUENCING/MsigDB/msigdb_v7.2_GMTs/c5.all.v7.2.symbols.gmt"
pathways_C1 <-gmtPathways(GMT_C1)
pathways_C2 <-gmtPathways(GMT_C2)
pathways_C3 <-gmtPathways(GMT_C3)
pathways_C5 <-gmtPathways(GMT_C5)

##5.1 High vs None
Table1 <- read.table("RESULTS/AB/SIG_Genes_Table1.txt", header=TRUE, sep = "\t"); head(Table1)
gene_list = Table1$t
names(gene_list) = Table1$Gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
head(gene_list)
gene_list <- gene_list[!is.na(gene_list)]

head(Table1)

## ANALYSIS by Chr location
fgseaRes_1_C1 <-fgsea(pathways=pathways_C1, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_1_C1 <- fgseaRes_1_C1[order(padj),] 
sum(fgseaRes_1_C1[, padj < 0.05])
write.table(as.matrix(fgseaRes_1_C1), file = "GSEA_Results/GSEA_C1_1_HighvsNone.txt", sep = "\t")

## ANALYSIS by Curated Gene Sets
fgseaRes_1_C2 <- fgsea(pathways=pathways_C2, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_1_C2 <- fgseaRes_1_C2[order(padj),] 
sum(fgseaRes_1_C2[, padj < 0.05])
write.table(as.matrix(fgseaRes_1_C2), file = "GSEA_Results/GSEA_C2_1_HighvsNone.txt", sep = "\t")

## ANALYSIS by MIR and TFTs
fgseaRes_1_C3 <-fgsea(pathways=pathways_C3, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_1_C3 <- fgseaRes_1_C3[order(padj),] 
sum(fgseaRes_1_C3[, padj < 0.05])
write.table(as.matrix(fgseaRes_1_C3), file = "GSEA_Results/GSEA_C3_1_HighvsNone.txt", sep = "\t")

## ANALYSIS by Gene Ontology
fgseaRes_1_C5 <-fgsea(pathways=pathways_C5, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_1_C5 <- fgseaRes_1_C5[order(padj),] 
sum(fgseaRes_1_C5[, padj < 0.05])
write.table(as.matrix(fgseaRes_1_C5), file = "GSEA_Results/GSEA_C5_1_HighvsNone.txt", sep = "\t")


##5.2 Medium vs None
Table3 <- read.table("RESULTS/AB/SIG_Genes_Table3.txt", header=TRUE, sep = "\t"); head(Table3)
gene_list = Table3$t
names(gene_list) = Table3$Gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

## ANALYSIS by Chr location
fgseaRes_3_C1 <-fgsea(pathways=pathways_C1, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_3_C1 <- fgseaRes_3_C1[!is.na(fgseaRes_3_C1$padj)]
fgseaRes_3_C1 <- fgseaRes_3_C1[order(padj),] 
sum(fgseaRes_3_C1[, padj < 0.05])
write.table(as.matrix(fgseaRes_3_C1), file = "GSEA_Results/GSEA_C1_3_MedvsNone.txt", sep = "\t")

## ANALYSIS by Curated Gene Sets
fgseaRes_3_C2 <-fgsea(pathways=pathways_C2, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_3_C2 <- fgseaRes_3_C2[!is.na(fgseaRes_3_C2$padj)]
fgseaRes_3_C2 <- fgseaRes_3_C2[order(padj),] 
sum(fgseaRes_3_C2[, padj < 0.05])
write.table(as.matrix(fgseaRes_3_C2), file = "GSEA_Results/GSEA_C2_3_MedvsNone.txt", sep = "\t")

## ANALYSIS by MIR and TFTs
fgseaRes_3_C3 <-fgsea(pathways=pathways_C3, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_3_C3 <- fgseaRes_3_C3[!is.na(fgseaRes_3_C3$padj)]
fgseaRes_3_C3 <- fgseaRes_3_C3[order(padj),] 
sum(fgseaRes_3_C3[, padj < 0.05])
write.table(as.matrix(fgseaRes_3_C3), file = "GSEA_Results/GSEA_C3_3_MedvsNone.txt", sep = "\t")

## ANALYSIS by Gene Ontology
fgseaRes_3_C5 <-fgsea(pathways=pathways_C5, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_3_C5<- fgseaRes_3_C5[!is.na(fgseaRes_3_C5$padj)]
fgseaRes_3_C5[order(padj),] 
sum(fgseaRes_3_C5[, padj < 0.05])
write.table(as.matrix(fgseaRes_3_C5), file = "GSEA_Results/GSEA_C5_3_MedvsNone.txt", sep = "\t")


##5.3 LOW vs None
Table4 <- read.table("RESULTS/AB/SIG_Genes_Table4.txt", header=TRUE, sep = "\t"); head(Table4)
gene_list = Table4$t
names(gene_list) = Table4$Gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

## ANALYSIS by Chr location
fgseaRes_4_C1 <- fgsea(pathways=pathways_C1, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_4_C1 <- fgseaRes_4_C1[!is.na(fgseaRes_4_C1$padj)]
fgseaRes_4_C1 <- fgseaRes_4_C1[order(padj),] 
sum(fgseaRes_4_C1[, padj < 0.05])
write.table(as.matrix(fgseaRes_4_C1), file = "GSEA_Results/GSEA_C1_4_LowvsNone.txt", sep = "\t")

## ANALYSIS by Curated Gene Sets
fgseaRes_4_C2 <- fgsea(pathways=pathways_C2, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_4_C2 <- fgseaRes_4_C2[!is.na(fgseaRes_4_C2$padj)]
fgseaRes_4_C2 <- fgseaRes_4_C2[order(padj),] 
sum(fgseaRes_4_C2[, padj < 0.05])
write.table(as.matrix(fgseaRes_4_C2), file = "GSEA_Results/GSEA_C2_4_LowvsNone.txt", sep = "\t")

## ANALYSIS by MIR and TFTs
fgseaRes_4_C3 <- fgsea(pathways=pathways_C3, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_4_C3 <- fgseaRes_4_C3[!is.na(fgseaRes_4_C3$padj)]
fgseaRes_4_C3 <- fgseaRes_4_C3[order(padj),] 
sum(fgseaRes_4_C3[, padj < 0.05])
write.table(as.matrix(fgseaRes_4_C3), file = "GSEA_Results/GSEA_C3_4_LowvsNone.txt", sep = "\t")

## ANALYSIS by Gene Ontology
fgseaRes_4_C5 <- fgsea(pathways=pathways_C5, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_4_C5 <- fgseaRes_4_C5[!is.na(fgseaRes_4_C5$padj)]
fgseaRes_4_C5 <- fgseaRes_4_C5[order(padj),] 
sum(fgseaRes_4_C5[, padj < 0.05])
write.table(as.matrix(fgseaRes_4_C5), file = "GSEA_Results/GSEA_C5_4_LowvsNone.txt", sep = "\t")

## 5.4 NONE Neg vs Pos
Table6 <- read.table("RESULTS/AB/SIG_Genes_Table6.txt", header=TRUE, sep = "\t"); head(Table6)
gene_list = Table6$t
names(gene_list) = Table6$Gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]


## ANALYSIS by Chr location
fgseaRes_6_C1 <-fgsea(pathways=pathways_C1, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_6_C1 <- fgseaRes_6_C1[!is.na(fgseaRes_6_C1$padj)]
fgseaRes_6_C1 <- fgseaRes_6_C1[order(padj),] 
sum(fgseaRes_6_C1[, padj < 0.05])
write.table(as.matrix(fgseaRes_6_C1), file = "GSEA_Results/GSEA_C1_6_NegvsPos.txt", sep = "\t")

## ANALYSIS by Curated Gene Sets
fgseaRes_6_C2 <- fgseaMultilevel(pathways=pathways_C2, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_6_C2 <- fgseaRes_6_C2[!is.na(fgseaRes_6_C2$padj)]
fgseaRes_6_C2 <- fgseaRes_6_C2[order(padj),] 
sum(fgseaRes_6_C2[, padj < 0.05])
write.table(as.matrix(fgseaRes_6_C2), file = "GSEA_Results/GSEA_C2_6_NegvsPos.txt", sep = "\t")
# fg_Col_6_C2<-collapsePathways (fgseaRes_6_C2, pathways=pathways_C2, stats=gene_list, pval.threshold=0.05, nperm=100000, gseaParam=1)

## ANALYSIS by MIR and TFTs
fgseaRes_6_C3 <- fgsea(pathways=pathways_C3, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_6_C3 <- fgseaRes_6_C3[!is.na(fgseaRes_6_C3$padj)]
fgseaRes_6_C3 <- fgseaRes_6_C3[order(padj),] 
sum(fgseaRes_6_C3[, padj < 0.05])
write.table(as.matrix(fgseaRes_6_C3), file = "GSEA_Results/GSEA_C3_6_NegvsPos.txt", sep = "\t")

## ANALYSIS by Gene Ontology
fgseaRes_6_C5 <- fgsea(pathways=pathways_C5, stats=gene_list, eps = 0, minSize = 15, maxSize = 500, gseaParam=1)
fgseaRes_6_C5 <- fgseaRes_6_C5[!is.na(fgseaRes_6_C5$padj)]
fgseaRes_6_C5 <- fgseaRes_6_C5[order(padj),] 
sum(fgseaRes_6_C5[, padj < 0.05])
write.table(as.matrix(fgseaRes_6_C5), file = "GSEA_Results/GSEA_C5_6_NegvsPos.txt", sep = "\t")

fgseaRes_1_C3[which(fgseaRes_1_C3$pathway == "MIR2392"), ]
fgseaRes_3_C3[which(fgseaRes_3_C3$pathway == "MIR2392"), ]
fgseaRes_4_C3[which(fgseaRes_4_C3$pathway == "MIR2392"), ]
fgseaRes_6_C3[which(fgseaRes_6_C3$pathway == "MIR2392"), ]

##Filter expression matrix and the reassess volcano pots for the pathway. 
##########HERE############

getwd()
mir2392_targets <- read.table("mir2392_Targets.txt")

res1 <- data.frame((Table1[,c(7,11)]));
rownames(res1) <- make.names(Table1$GeneSymbol, unique=TRUE) 
colnames(res1)<- c("log2FC","adj.P.Value")
datarows = match(mir2392_targets[,1], rownames(res1))
new <- res1[datarows,]

# now draw the figure.....
tiff(file="PLOTS-ALL/Volcano-High COVID vs Neg.tiff", unit= "in", res = 300, width = 5, height = 5)
EnhancedVolcano(new,lab = rownames(new),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "HIGH COVID vs Neg"); dev.off()

res3 <- data.frame((Table3[,c(7,11)]));
rownames(res3) <- make.names(Table1$GeneSymbol, unique=TRUE) 
colnames(res3)<- c("log2FC","adj.P.Value")
datarows = match(mir2392_targets[,1], rownames(res3))
new <- res3[datarows,]


# now draw the figure.....
#tiff(file="PLOTS-ALL/Volcano-High COVID vs Neg.tiff", unit= "in", res = 300, width = 5, height = 5)
EnhancedVolcano(new,lab = rownames(new),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "MED COVID vs Neg")
#dev.off()


res4 <- data.frame((Table4[,c(7,11)]));
rownames(res4) <- make.names(Table1$GeneSymbol, unique=TRUE) 
colnames(res4)<- c("log2FC","adj.P.Value")
datarows = match(mir2392_targets[,1], rownames(res4))
new <- res4[datarows,]

# now draw the figure.....
#tiff(file="PLOTS-ALL/Volcano-High COVID vs Neg.tiff", unit= "in", res = 300, width = 5, height = 5)
EnhancedVolcano(new,lab = rownames(new),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "LOW COVID vs Neg")
#dev.off()

res6 <- data.frame((Table6[,c(7,11)]));
rownames(res6) <- make.names(Table1$GeneSymbol, unique=TRUE) 
colnames(res6)<- c("log2FC","adj.P.Value")
datarows = match(mir2392_targets[,1], rownames(res6))
new <- res6[datarows,]

# now draw the figure.....
#tiff(file="PLOTS-ALL/Volcano-High COVID vs Neg.tiff", unit= "in", res = 300, width = 5, height = 5)
EnhancedVolcano(new,lab = rownames(new),
    x = 'log2FC',
    y = 'adj.P.Value',
    xlim = c(-5, 8), title = "LOW COVID vs Neg")
#dev.off()



# probably want to review using excel then look for individual pathways. 




plotEnrichment(pathways_C5[["GO_RESPONSE_TO_VIRUS"]], gene_list) + labs(title="GO_RESPONSE_TO_VIRUS  ")

# search for a single
fgseaRes_MIR[which(fgseaRes_MIR$pathway == "MIR2392"), ]


#### Venn diagram of pathways vs NONE Neg/Pos

venn.diagram(list(High = fgseaRes_1_C2$pathway[fgseaRes_1_C2$padj < 0.05], Med = fgseaRes_3_C2$pathway[fgseaRes_3_C2$padj < 0.05], Low = fgseaRes_4_C2$pathway[fgseaRes_4_C2$padj < 0.05]), filename = "GSEA_RESULTS/GSEA_HighvsLow-C3.tiff", height = 3000, width = 3000, resolution =
    300, imagetype = "tiff", units = "px", main = " High vs Med GO Terms", 
main.cex = 3, cat.fontface="bold", scaled = TRUE, label.cex = 3, cat.pos = c(330, 30, 0), cat.cex = 2.5, fill = c("red", "lightblue", "yellow"), alpha = 0.5, cex = 3)





########################################################################################


session.info()

