##############################################################################################################################
#This script carries out the steps involved in analysis of RNAseq data.  
#depending on your data and interests, different parts of this script may not apply to you
##############################################################################################################################
#begin by loading the packages required for RNAseq data
library(Rsubread)
library(limma)
library(edgeR)
library(ShortRead)
options(digits=2)
#library(Biostrings)

#read in your study design file
targets <- read.delim("studyDesign.alt.txt", row.names=NULL)
targets
groups <- factor(paste(targets$strain, targets$stage, sep="."))
pair <- factor(targets$batch)
#create some more human-readable labels for your samples using the info in this file
sampleLabels <- paste(targets$strain, targets$stage, targets$replicate, sep=".")

#set-up your experimental design
design <- model.matrix(~0+groups)
design2 <- model.matrix(~0+groups+pair)

colnames(design) <- levels(groups)
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
load("DGEList")
#retrieve all your gene/transcript identifiers from this DGEList object
myGeneIDs <- DGEList$genes$GeneID

##############################################################################################################################
#Normalize unfiltered data using 'voom' function in Limma package
#This will normalize based on the mean-variance relationship
#will also generate the log2 of counts per million based on the size of each library (also a form of normalization)
##############################################################################################################################
normData.unfiltered <- voom(DGEList, design, plot=TRUE)
exprs.unfiltered <- normData.unfiltered$E
exprs.matrix.unfiltered <- as.matrix(exprs.unfiltered)
#note that because you're now working with Log2 CPM, much of your data will be negative number (log2 of number smaller than 1 is negative) 
head(exprs.matrix.unfiltered)

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

normData.filtered <- voom(DGEList.filtered, design, plot=TRUE)
exprs.filtered <- normData.filtered$E
exprs.matrix.filtered <- as.matrix(exprs.filtered)
#note that because you're now working with Log2 CPM, much of your data will be negative number (log2 of number smaller than 1 is negative) 
head(exprs.matrix.filtered)

rpkm.filtered <- rpkm(DGEList.filtered, DGEList.filtered$genes$Length) #if you prefer, can use 'cpm' instead of 'rpkm' here
rpkm.filtered <- log2(rpkm.filtered + 1)

###############################################################################################
#explore your data using some standard approaches
###############################################################################################
#choose color scheme for graphs
cols.ALL <- topo.colors (n=18, alpha=1)
hist(rpkm.filtered, xlab = "log2 expression", main = "normalized data - histograms", col=cols.ALL)
boxplot(rpkm.filtered, ylab = "normalized log2 expression", main = "non-normalized data - boxplots", col=cols.ALL)

distance <- dist(t(rpkm.filtered),method="euclidean") # options for computing distance matrix are: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". 
clusters <- hclust(distance, method = "average") #options for clustering method: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
plot(clusters, label=sampleLabels)

library(limma)
library(statmod)
batch <- factor(targets$batch)
corfit <- duplicateCorrelation(rpkm.filtered, design=design, block = batch)
fit <- lmFit(rpkm.filtered, block = batch, cor = corfit$consensus)
fit <- eBayes(fit)
topTable(fit, adjust = "BH")
##############################################################################################################################
#annotate your normalized data using the organism-specific database package
##############################################################################################################################
library(org.Tgondii.eg.db)
library(AnnotationDbi)
#If we want to know what kinds of data are retriveable via the 'select' command, look at the columns of the annotation database
columns(org.Tgondii.eg.db)
#If we want to know what kinds of fields we could potentially use as keys to query the database, use the 'keytypes' command
keytypes(org.Tgondii.eg.db)
#transform you identifiers to entrezIDs
myAnnot.unfiltered <- select(org.Tgondii.eg.db, keys=rownames(exprs.matrix.unfiltered), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL", "GENENAME"))
myAnnot.filtered <- select(org.Tgondii.eg.db, keys=rownames(exprs.matrix.filtered), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL", "GENENAME"))
resultTable.unfiltered <- merge(myAnnot.unfiltered, exprs.matrix.unfiltered, by.x="ENSEMBL", by.y=0)
resultTable.filtered <- merge(myAnnot.filtered, exprs.matrix.filtered, by.x="ENSEMBL", by.y=0)
head(resultTable.unfiltered)
#add more appropriate sample names as column headers
colnames(resultTable.unfiltered) <- sampleLabels
colnames(resultTable.filtered) <- sampleLabels
#now write these annotated datasets out
write.table(resultTable.unfiltered, "normalizedUnfiltered.txt", sep="\t", quote=FALSE)
write.table(resultTable.filtered, "normalizedFiltered.txt", sep="\t", quote=FALSE)
#pull out your rownames of your dataset to use as one of the filters
myEnsemblIDs.filtered <- rownames(rpkm.filtered)
myEnsemblIDs.unfiltered <- rownames(rpkm.unfiltered)

#transform your identifiers to entrezIDs
resultTable.filtered <- merge(myAnnot.filtered, rpkm.filtered, by.x="ensembl_gene_id", by.y=0)
resultTable.unfiltered <- merge(myAnnot.unfiltered, rpkm.unfiltered, by.x="ensembl_gene_id", by.y=0)


#this script walks thorough some basic data wrangling for organizing expression data spreadsheets and ends with
#how to create publication quality graphics from transcriptomic data generated (regardless of platform used)
#to start this script you need a file with all your expression data and some non-redundant identifiers as row names (usually gene symbols)
#you also need a study design file

#for creating graphics, we'll use the ggplot2 and ggvis packages which employ a 'grammar of graphics' approach
#load the packages
library(ggplot2)
library(dplyr)
library(ggvis)

#clean-up data table
colnames(resultTable.filtered) <- c("geneID", "symbol", "description", sampleLabels)
colnames(resultTable.filtered)
#use the dplyr 'mutate' command to get averages and fold changes for all your replicates
myData <- mutate(resultTable.filtered,
                 Ecdysone.vs.PBS_5hr_carcass = (ecdysone.5hr.carcass - PBS.5hr.carcass),
                 Ecdysone.vs.PBS_18hr_carcass = (ecdysone.18hr.carcass - PBS.18hr.carcass),
                 Ecdysone.vs.PBS_5hr_gut = (ecdysone.5hr.midgut - PBS.5hr.midgut),
                 Ecdysone.vs.PBS_18hr_gut = (ecdysone.18hr.midgut - PBS.18hr.midgut))


#use dplyr "arrange" and "select" functions to sort by LogFC column of interest (arrange) 
#and then display only the columns of interest (select) to see the most differentially expressed genes
myData.sort <- myData %>%
  arrange(desc(Ecdysone.vs.PBS_18hr_gut)) %>%
  select(geneID, description, Ecdysone.vs.PBS_18hr_gut, Ecdysone.vs.PBS_18hr_carcass)
head(myData.sort)

#use dplyr "filter" and "select" functions to pick out genes of interest (filter)
#and again display only columns of interest (select)
myData.filter <- myData %>%
  filter(geneID=="AAEL010434" | geneID=="AAEL009600" | geneID == "AAEL001414") %>%
  select(geneID, Ecdysone.vs.PBS_5hr_carcass, Ecdysone.vs.PBS_18hr_carcass, Ecdysone.vs.PBS_5hr_gut, Ecdysone.vs.PBS_18hr_gut)
head(myData.filter)

myData.filter <- myData %>%
  filter((abs(Ecdysone.vs.PBS_18hr_gut) >= 1) | (abs(Ecdysone.vs.PBS_5hr_gut) >= 1))%>%
  select(geneID, Ecdysone.vs.PBS_5hr_carcass, Ecdysone.vs.PBS_18hr_carcass, Ecdysone.vs.PBS_5hr_gut, Ecdysone.vs.PBS_18hr_gut)
head(myData.filter)

#make simple line or bar graph
#first reorder and clean up the data you want to graph
library(dplyr)
row.names(myData.filter) <- myData.filter[,1]
myData.filter <- select(myData.filter, -geneSymbols)
myData.filter.transpose <- as.data.frame(t(myData.filter))
treatments <- row.names(myData.filter.transpose)
myData.filter.transpose <- mutate(myData.filter.transpose, treatments, 
                                  genotype=c("WT", "het","KO"))
myData.filter.transpose



