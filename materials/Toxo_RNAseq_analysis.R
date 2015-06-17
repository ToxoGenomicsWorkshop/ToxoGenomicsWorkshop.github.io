##############################################################################################################################
#This script carries out the steps involved in analysis of RNAseq data.  
#depending on your data and interests, different parts of this script may not apply to you
##############################################################################################################################
#begin by downloading the packages you'll need for this analysis
source("http://bioconductor.org/biocLite.R")
biocLite(pkgs=c("Rsubread", "limma", "edgeR", "ShortRead", "ggvis", "ggplot", "reshape2", "dplyr")
?biocLite

#begin by loading the packages required for RNAseq data
library(Rsubread)
library(limma)
library(edgeR)
library(ShortRead)
options(digits=2)
#library(Biostrings)

#read in your study design file
targets <- read.delim("studyDesign.txt", row.names=NULL)
targets
groups <- factor(paste(targets$strain, targets$stage, sep="."))
batch <- factor(targets$rep)
#create some more human-readable labels for your samples using the info in this file
sampleLabels <- paste(targets$strain, targets$stage, targets$rep, sep=".")

#set-up your experimental design
design <- model.matrix(~0+groups)
#design2 <- model.matrix(~batch+groups)

colnames(design) <- c("intercept", levels(groups)
colnames(design2) <- levels(groups)
design

# ##############################################################################################################################
# #you can check read quality using shortRead package
# #but I usually find it is better to do this on the sequencer or Illumina's BaseSpace website
# ##############################################################################################################################
# myFastq <- targets$fastq
# #collecting statistics over the files
# qaSummary <- qa("B6-WT-untreat-rep2_S2_mergedLanes_read1.fastq", type="fastq")
# #create and view a report
# browseURL(report(qaSummary))

##############################################################################################################################
#build index from your reference genome (expect this to take about 20 min on 8G RAM for mouse genome)
#you must have already downloaded the fasta file for your genome of interest and have it in the working directory
#this only needs to be done once, then index can be reused for future alignments
##############################################################################################################################
buildindex(basename="TgME49",reference="ToxoDB-24_TgondiiME49_Genome.fasta", colorspace=TRUE)

##############################################################################################################################
#align your reads (in the fastq files) to your indexed reference genome that you created above
#expect this to take about 45min for a single fastq file containing 25 million reads
#the output from this is a .bam file for each of your original fastq files
##############################################################################################################################
reads <- targets$fastq[13:18]
align(index="TgME49", readfile1=reads, input_format="gzFASTQ",output_format="BAM", color2base=TRUE,
      output_file=targets$output[13:18], tieBreakHamming=TRUE, unique=TRUE, indels=5, nthreads=8)

##############################################################################################################################
#use the 'featureCounts' function to summarize read counts to genomic features (exons, genes, etc)
#will take about 1-2min per .bam file.
#for total transcriptome data summarized to mouse/human .gtf, expect about 50-60% of reads summarize to genes (rest is non-coding)
##############################################################################################################################
#read in text file with bam file names
myBAM <- targets$output
#summarize aligned reads to genomic features (i.e. exons)
fc <- featureCounts(files=myBAM, annot.ext="ToxoDB-24_TgondiiME49.gff", 
                    isGTFAnnotationFile=TRUE, useMetaFeatures=TRUE, strandSpecific=1, 
                    GTF.attrType="ID", GTF.featureType="gene", nthreads=8)
#use the 'DGEList' function from EdgeR to make a 'digital gene expression list' object
DGEList <- DGEList(counts=fc$counts, genes=fc$annotation)
save(DGEList, file="DGEList")



##############################################################################################################################
#beginning of Toxo13 Genomics Workshop
##############################################################################################################################

load("DGEList")
#retrieve all your gene/transcript identifiers from this DGEList object
myGeneIDs <- DGEList$genes$GeneID

##############################################################################################################################
#Normalize unfiltered data using 'voom' function in Limma package
#This will normalize based on the mean-variance relationship
#will also generate the log2 of counts per million based on the size of each library (also a form of normalization)
##############################################################################################################################
normData.unfiltered <- voom(DGEList, design2, plot=TRUE)
exprs.unfiltered <- normData.unfiltered$E
#note that because you're now working with Log2 CPM, much of your data will be negative number (log2 of number smaller than 1 is negative) 
head(exprs.unfiltered)

#if you need RPKM for your unfiltered, they can generated as follows
#Although RPKM are commonly used, not really necessary since you don't care to compare two different genes within a sample
rpkm.unfiltered <- rpkm(DGEList, DGEList$genes$Length)
rpkm.unfiltered <- log2(rpkm.unfiltered + 1)

##############################################################################################################################
#Filtering your dataset and normalize this
#Only keep in the analysis those genes which had >10 reads per million mapped reads in at least two libraries.
##############################################################################################################################
cpm.matrix.filtered <- rowSums(cpm(DGEList) > 10) >= 2
DGEList.filtered <- DGEList[cpm.matrix.filtered,]
dim(DGEList.filtered)

normData.filtered <- voom(DGEList.filtered, design2, plot=TRUE)
exprs.filtered <- normData.filtered$E
#note that because you're now working with Log2 CPM, much of your data will be negative number (log2 of number smaller than 1 is negative) 
head(exprs.filtered)

rpkm.filtered <- rpkm(DGEList.filtered, DGEList.filtered$genes$Length) #if you prefer, can use 'cpm' instead of 'rpkm' here
rpkm.filtered <- log2(rpkm.filtered + 1)

###############################################################################################
#explore your data using some standard approaches
###############################################################################################
#choose color scheme for graphs
cols.ALL <- topo.colors (n=18, alpha=1)
hist(exprs.filtered, xlab = "log2 expression", main = "normalized data - histograms", col=cols.ALL)
boxplot(exprs.filtered, ylab = "normalized log2 expression", main = "non-normalized data - boxplots", col=cols.ALL)

distance <- dist(t(exprs.filtered),method="maximum") # options for computing distance matrix are: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". 
clusters <- hclust(distance, method = "average") #options for clustering method: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
plot(clusters, label=sampleLabels)

###############################################################################################
#Principal component analysis of the filtered data matrix
###############################################################################################
pca.res <- prcomp(t(exprs.filtered), scale.=F, retx=T)
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
head(pca.res$rotation) #$rotation shows you how much each GENE influenced each PC (callled 'eigenvalues', or loadings)
head(pca.res$x) #$x shows you how much each SAMPLE influenced each PC (called 'scores')
plot(pca.res, las=1)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

#make some graphs to visualize your PCA result
library(ggplot2)
#lets first plot any two PCs against each other
#turn your scores for each gene into a data frame
data.frame <- as.data.frame(pca.res$x)
ggplot(data.frame, aes(x=PC1, y=PC2, colour=factor(groups))) +
  geom_point(size=5) +
  theme(legend.position="right")

#create a 'small multiples' chart to look at impact of each variable on each pricipal component
library(reshape2)
melted <- cbind(groups, melt(pca.res$x[,1:3]))
#look at your 'melted' data
ggplot(melted) +
  geom_bar(aes(x=Var1, y=value, fill=groups), stat="identity") +
  facet_wrap(~Var2)


###############################################################################################
#Cleaning up and graphing your data
###############################################################################################
#for creating graphics, we'll use the ggplot2 and ggvis packages which employ a 'grammar of graphics' approach
#load the packages
library(ggplot2)
library(dplyr)
library(ggvis)

#clean-up data table
exprs.filtered.dataframe <- as.data.frame(exprs.filtered)
colnames(exprs.filtered.dataframe)
head(exprs.filtered.dataframe)
colnames(exprs.filtered.dataframe) <- sampleLabels
geneID <- rownames(exprs.filtered.dataframe)
#use the dplyr 'mutate' command to get averages and fold changes for all your replicates
myData <- mutate(exprs.filtered.dataframe,
                 #insert columns that average your replicates
                 RH.tachy.AVG = (RH.tachy.rep1 + RH.tachy.rep2 + RH.tachy.rep3)/2,
                 PLK.tachy.AVG = (PLK.tachy.rep1 + PLK.tachy.rep2 + PLK.tachy.rep3)/2,
                 CTG.tachy.AVG = (CTG.tachy.rep1 + CTG.tachy.rep2 + CTG.tachy.rep3)/2,
                 RH.brady.AVG = (RH.brady.rep1 + RH.brady.rep2 + RH.brady.rep3)/2,
                 PLK.brady.AVG = (PLK.brady.rep1 + PLK.brady.rep2 + PLK.brady.rep3)/2,
                 CTG.brady.AVG = (CTG.brady.rep1 + CTG.brady.rep2 + CTG.brady.rep3)/2,
                 #now add fold-change columns based on the averages calculated above
                 PLK.vs.RH.tachy = (PLK.tachy.AVG - RH.tachy.AVG),
                 CTG.vs.RH.tachy = (CTG.tachy.AVG - RH.tachy.AVG),
                 PLK.vs.CTG.tachy = (PLK.tachy.AVG - CTG.tachy.AVG),
                 brady.vs.tachy.RH = (RH.brady.AVG - RH.tachy.AVG),
                 brady.vs.tachy.PLK = (PLK.brady.AVG - PLK.tachy.AVG),
                 brady.vs.tachy.CTG = (CTG.brady.AVG - CTG.tachy.AVG),
                 geneID)

#take a look at your new spreadsheet
head(myData)
#use dplyr "arrange" and "select" functions to sort by LogFC column of interest (arrange) 
#and then display only the columns of interest (select) to see the most differentially expressed genes
myData.sort <- myData %>%
  arrange(desc(brady.vs.tachy.PLK)) %>%
  select(geneID, PLK.tachy.AVG, PLK.brady.AVG)
head(myData.sort)


#use dplyr "filter" and "select" functions to pick out genes of interest (filter)
#and again display only columns of interest (select)
#filter based on specific Toxo gene IDs
myData.filter <- myData %>%
  filter(geneID=="TGME49_207130" | geneID=="TGME49_208130") %>%
  select(geneID, PLK.tachy.AVG, PLK.brady.AVG)
head(myData.filter)


#filtering based on expression level or fold change
myData.filter <- myData %>%
  filter((abs(Ecdysone.vs.PBS_18hr_gut) >= 1) | (abs(Ecdysone.vs.PBS_5hr_gut) >= 1))%>%
  select(geneID, PLK.tachy.AVG, PLK.brady.AVG)
head(myData.filter)

#create a basic scatterplot using ggplot
ggplot(myData, aes(x=B1a.Ph0.AVG, y=Mac.Ph0.AVG)) +
  geom_point(shape=1) +
  geom_point(size=4)

#define a tooltip that shows gene symbol and Log2 expression data when you mouse over each data point in the plot
tooltip <- function(data, ...) {
  paste0("<b>","Symbol: ", data$geneSymbols, "</b><br>",
         "B1a.Ph0.AVG: ", data$B1a.Ph0.AVG, "<br>",
         "B1a.nonPh0.AVG: ", data$B1a.nonPh0.AVG)
}

#plot the interactive graphic
myData %>% 
  ggvis(x= ~B1a.Ph0.AVG, y= ~B1a.nonPh0.AVG, key := ~geneSymbols) %>% 
  layer_points(fill = ~LogFC.B1a.Ph0.vs.B1a.nonPh0) %>%
  add_tooltip(tooltip)


###############################################################################################
# use Limma to find differentially expressed genes between two or more conditions
###############################################################################################
# fit the linear model to your filtered expression data
library(limma)
fit <- lmFit(exprs.filtered, design)
#fit2 <- lmFit(exprs.filtered, design2)

# set up a contrast matrix based on the pairwise comparisons of interest
contrast.matrix.strains <- makeContrasts(RHvsCTG = RH.tachy - CTG.tachy, RHvsPLK = RH.tachy - PLK.brady, CTGvsPLK = CTG.tachy - PLK.brady, levels=design)

# check each contrast matrix
contrast.matrix.strains

# extract the linear model fit for the contrast matrix that you just defined above
fits <- contrasts.fit(fit, contrast.matrix.strains)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)


###############################################################################################
# use the topTable and decideTests functions to see the differentially expressed genes
###############################################################################################

# use topTable function to take a look at the hits
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")
myTopHits

# use the 'decideTests' function to show Venn diagram for all diffexp genes for up to three comparisons
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=1)
#stats <- write.fit(ebFit)
vennDiagram(results, include="both") #all pairwise comparisons on a B6 background


# take a look at what the results of decideTests looks like
results

# now pull out probeIDs from selected regions of the Venn diagram.  In this case, I want all genes in the venn.
diffGenes <- which(results[,1] !=0 | results[,2] !=0 | results[,3] !=0)

# retrieve expression data for the probes from above
diffData <- exprs.filtered[results[,1] !=0 | results[,2] !=0 | results[,3] !=0]

#combine probeIDs, gene symbols and expression data for differentially expressed genes into one file
write.table(cbind(diffGenes, diffData),"DiffGenes.xls", sep="\t", quote=FALSE)

# take a look at each expression matrix
dim(diffData)


# ##############################################################################################################################
# #annotate your normalized data using the organism-specific database package
# ##############################################################################################################################
# library(org.Tgondii.eg.db)
# library(AnnotationDbi)
# #If we want to know what kinds of data are retriveable via the 'select' command, look at the columns of the annotation database
# columns(org.Tgondii.eg.db)
# #If we want to know what kinds of fields we could potentially use as keys to query the database, use the 'keytypes' command
# keytypes(org.Tgondii.eg.db)
# #transform you identifiers to entrezIDs
# myAnnot.unfiltered <- select(org.Tgondii.eg.db, keys=rownames(exprs.matrix.unfiltered), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL", "GENENAME"))
# myAnnot.filtered <- select(org.Tgondii.eg.db, keys=rownames(exprs.matrix.filtered), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL", "GENENAME"))
# resultTable.unfiltered <- merge(myAnnot.unfiltered, exprs.matrix.unfiltered, by.x="ENSEMBL", by.y=0)
# resultTable.filtered <- merge(myAnnot.filtered, exprs.matrix.filtered, by.x="ENSEMBL", by.y=0)
# head(resultTable.unfiltered)
# #add more appropriate sample names as column headers
# colnames(resultTable.unfiltered) <- sampleLabels
# colnames(resultTable.filtered) <- sampleLabels
# #now write these annotated datasets out
# write.table(resultTable.unfiltered, "normalizedUnfiltered.txt", sep="\t", quote=FALSE)
# write.table(resultTable.filtered, "normalizedFiltered.txt", sep="\t", quote=FALSE)
# #pull out your rownames of your dataset to use as one of the filters
# myEnsemblIDs.filtered <- rownames(rpkm.filtered)
# myEnsemblIDs.unfiltered <- rownames(rpkm.unfiltered)
# 
# #transform your identifiers to entrezIDs
# resultTable.filtered <- merge(myAnnot.filtered, rpkm.filtered, by.x="ensembl_gene_id", by.y=0)
# resultTable.unfiltered <- merge(myAnnot.unfiltered, rpkm.unfiltered, by.x="ensembl_gene_id", by.y=0)

