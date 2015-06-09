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
targets <- read.delim("Knoll_studyDesign.txt")
targets
groups <- factor(paste(targets$timepoint, sep="."))
#create some more human-readable labels for your samples using the info in this file
sampleLabels <- paste(targets$timepoint, targets$replicate, sep=".")

#set-up your experimental design
design <- model.matrix(~0+groups)
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
buildindex(basename="mouse",reference="Mus_musculus.GRCm38.dna.primary_assembly.fa")

##############################################################################################################################
#align your reads (in the fastq files) to your indexed reference genome that you created above
#expect this to take about 45min for a single fastq file containing 25 million reads
#the output from this is a .bam file for each of your original fastq files
##############################################################################################################################
reads_F <- targets$fastq[,1:9]
reads_R <- targets$fastq[,10:18]
align(index="TgME49", readfile1=reads1, readfile2=reads2, input_format="gzFASTQ",output_format="BAM",
      output_file="alignmentResultsPE_sample9.BAM", tieBreakHamming=TRUE, unique=TRUE, indels=5, nthreads=8)

##############################################################################################################################
#use the 'featureCounts' function to summarize read counts to genomic features (exons, genes, etc)
#will take about 1-2min per .bam file.
#for total transcriptome data summarized to mouse/human .gtf, expect about 50-60% of reads summarize to genes (rest is non-coding)
##############################################################################################################################
#read in text file with bam file names
MyBAM <- targets$Output
#summarize aligned reads to genomic features (i.e. exons)
fc <- featureCounts(files=MyBAM, annot.ext="Mus_musculus.GRCm38.79.gtf", isGTFAnnotationFile=TRUE, GTF.featureType = "exon",
                    GTF.attrType="gene_id", useMetaFeatures=TRUE, isPairedEnd=TRUE, requireBothEndsMapped=TRUE, strandSpecific=2, nthreads=8)
#use the 'DGEList' function from EdgeR to make a 'digital gene expression list' object
DGEList <- DGEList(counts=fc$counts, genes=fc$annotation)
save(DGEList, file="DGEList")
#retrieve all your gene/transcript identifiers from this DGEList object
myEnsemblIDs <- DGEList$genes$GeneID
dim(fc$counts)
tail(fc$counts)

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
tail(rpkm)

##############################################################################################################################
#Filtering your dataset
#Only keep in the analysis those genes which had >10 reads per million mapped reads in at least two libraries.
##############################################################################################################################
cpm.matrix.filtered <- rowSums(cpm(DGEList) > 10) >= 2
DGEList.filtered <- DGEList[cpm.matrix.filtered,]
dim(DGEList.filtered)
#Use Voom again to normalize this filtered dataset
normData <- voom(DGEList.filtered, design, plot=TRUE)
exprs <- normData$E
exprs.matrix.filtered <- as.matrix(exprs)
head(exprs.matrix.filtered)

##############################################################################################################################
#annotate your normalized data using the organism-specific database package
##############################################################################################################################
library(org.Mm.eg.db)
library(AnnotationDbi)
ls("package:org.Mm.eg.db")
#If we want to know what kinds of data are retriveable via the 'select' command, look at the columns of the annotation database
columns(org.Mm.eg.db)
#If we want to know what kinds of fields we could potentially use as keys to query the database, use the 'keytypes' command
keytypes(org.Mm.eg.db)
#transform you identifiers to entrezIDs
myAnnot.unfiltered <- select(org.Mm.eg.db, keys=rownames(exprs.matrix.unfiltered), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL", "GENENAME"))
myAnnot.filtered <- select(org.Mm.eg.db, keys=rownames(exprs.matrix.filtered), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL", "GENENAME"))
resultTable.unfiltered <- merge(myAnnot.unfiltered, exprs.matrix.unfiltered, by.x="ENSEMBL", by.y=0)
resultTable.filtered <- merge(myAnnot.filtered, exprs.matrix.filtered, by.x="ENSEMBL", by.y=0)
head(resultTable.unfiltered)
#add more appropriate sample names as column headers
colnames(resultTable.unfiltered) <- sampleLabels
colnames(resultTable.filtered) <- sampleLabels
#now write these annotated datasets out
write.table(resultTable.unfiltered, "normalizedUnfiltered.txt", sep="\t", quote=FALSE)
write.table(resultTable.filtered, "normalizedFiltered.txt", sep="\t", quote=FALSE)

###############################################################################################
#explore your data using some standard approaches
###############################################################################################
#choose color scheme for graphs
cols.ALL <- topo.colors (n=12, alpha=1)
hist(exprs, xlab = "log2 expression", main = "normalized data - histograms", col=cols.ALL)
boxplot(exprs, ylab = "normalized log2 expression", main = "non-normalized data - boxplots", col=cols.ALL)

distance <- dist(t(exprs.matrix),method="maximum") # options for computing distance matrix are: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". 
clusters <- hclust(distance, method = "complete") #options for clustering method: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
plot(clusters, label=sampleLabels)

###############################################################################################
#Principal component analysis of the filtered data matrix
###############################################################################################
pca.res <- prcomp(t(exprs.matrix), scale.=F, retx=T)
head(pca.res$rotation) #$rotation gives you the eigenvectors or 'loadings'
head(pca.res$x) #$x gives you the 'scores'
plot(pca.res, las=1)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per
plot.pch <- (1:length(levels(groups)))[as.numeric(groups)]
plot(pca.res$x, col=1, pch=plot.pch, las=1, cex=3, xlab=paste("PC1 (", pc.per [1], "%)", sep=""),ylab=paste("PC2 (", pc.per [2], "%)", sep=""), ylim=c(-30,30), xlim=c(-150,150))
text(pca.res$x, sampleLabels, pos= 1 )
grid()
legend(-40, 0, levels(groups), pch=1:length(levels(groups)), pt.cex=1.5)
summary(pca.res) # Prints variance summary for all principal components.
#3d PCA plot of PC1, PC2 and PC3
library(rgl)
rgl.open(); offset <- 50; par3d(windowRect=c(offset, offset, 640+offset, 640+offset)); rm(offset); rgl.clear(); rgl.viewpoint(theta=45, phi=30, fov=60, zoom=1); spheres3d(pca.res$x[,1], pca.res$x[,2], pca.res$x[,3], color=cols.ALL, radius=5, alpha=1, shininess=20); aspect3d(1, 1, 1); axes3d(col='black'); title3d("", "", "PC1", "PC2", "PC3", col='black'); bg3d("white") 

###############################################################################################
# use Limma to find differentially expressed genes between two or more conditions
###############################################################################################
# fit the linear model to the filtered expression set
library(limma)
fit <- lmFit(exprs.matrix, design)

# set up a contrast matrix based on the pairwise comparisons of interest
contrast.matrix.untxt <- makeContrasts(Ripk3_vs_B6 = Ripk3.untreated - B6.untreated, doubleKO_vs_B6 = Ripk3_Casp8.untreated - B6.untreated, doubleKO_vs_Ripk3 = Ripk3_Casp8.untreated - Ripk3.untreated, levels=design)
contrast.matrix.LPS6hr <- makeContrasts(Ripk3_vs_B6 = Ripk3.LPS_6hr - B6.LPS_6hr, doubleKO_vs_B6 = Ripk3_Casp8.LPS_6hr - B6.LPS_6hr, doubleKO_vs_Ripk3 = Ripk3_Casp8.LPS_6hr - Ripk3.LPS_6hr, levels=design)
contrast.matrix.LPSeffect <- makeContrasts(LPS_6hr = B6.LPS_6hr - B6.untreated, levels=design)

# check each contrast matrix
contrast.matrix.untxt
contrast.matrix.LPS6hr
contrast.matrix.LPSeffect

# extract the linear model fit for the contrast matrix that you just defined above
fits.untxt <- contrasts.fit(fit, contrast.matrix.untxt)
fits.LPS6hr <- contrasts.fit(fit, contrast.matrix.LPS6hr)
fits.LPSeffect <- contrasts.fit(fit, contrast.matrix.LPSeffect)

#ebFit.ALL <- eBayes(fits.ALL)
ebFit.untxt <- eBayes(fits.untxt)
ebFit.LPS6hr <- eBayes(fits.LPS6hr)
ebFit.LPSeffect <- eBayes(fits.LPSeffect)


###############################################################################################
# use the topTable and decideTests functions to see the differentially expressed genes
###############################################################################################

# use topTable function to take a look at the hits
probeList.untxt <- topTable(ebFit.untxt, adjust ="BH", coef=3, number=50, sort.by="logFC")
probeList.LPS6hr <- topTable(ebFit.LPS6hr, adjust ="BH", coef=3, number=200, sort.by="logFC")
probeList.LPSeffect <- topTable(ebFit.LPSeffect, adjust ="BH", coef=1, number=200, sort.by="logFC")

# use the 'decideTests' function to show Venn diagram for all diffexp genes for up to three comparisons
#results.ALL <- decideTests(ebFit.ALL, method="global", adjust.method="BH", p.value=0.01, lfc=1)
results.untxt <- decideTests(ebFit.untxt, method="global", adjust.method="BH", p.value=0.01, lfc=0.58)
results.LPS6hr<- decideTests(ebFit.LPS6hr, method="global", adjust.method="BH", p.value=0.01, lfc=0.58)
results.LPSeffect<- decideTests(ebFit.LPSeffect, method="global", adjust.method="BH", p.value=0.01, lfc=0.58)

vennDiagram(results.untxt, include="both") 
vennDiagram(results.LPS6hr, include="both") 
vennDiagram(results.LPSeffect, include="both") 

diffData.LPS6hr <- normData[results.LPS6hr[,1] !=0 | results.LPS6hr[,2] !=0 | results.LPS6hr[,3] !=0]
diffData.LPS6hr <- diffData.LPS6hr$E
head(diffData.LPS6hr)
dim(diffData.LPS6hr)
# print out results to a table
write.table(results.ALL, "diffGenes_6hr_2fold.txt", sep="\t", quote=FALSE)


# take a look at each expression matrix
dim(diffData.untxt)
dim(diffData.LPS2hr)
dim(diffData.LPS6hr)
dim(diffData.LPSeffect)

###################################################################################################
# generate a heatmap of differentially expressed transcripts using the entire dataset
###################################################################################################
#averaging the replicate arrays for all the data, then heatmap script below
head(diffData.untxt)
library(limma)
colnames(diffData.untxt) <- groups.untxt
colnames(diffData.LPS2hr) <- groups.LPS2hr
colnames(diffData.LPS6hr) <- groups.LPS6hr
colnames(diffData.LPSeffect) <- groups
diffData.untxt.AVG <- avearrays(diffData.untxt)
diffData.LPS2hr.AVG <- avearrays(diffData.LPS2hr)
diffData.LPS6hr.AVG <- avearrays(diffData.LPS6hr)
diffData.LPSeffect.AVG <- avearrays(diffData.LPSeffect)
head(diffData.untxt.AVG)
dim(diffData.untxt.AVG)

#cluster rows using pearson correlations
hr.untxt <- hclust(as.dist(1-cor(t(diffData.untxt.AVG), method="pearson")), method="complete") 
hr.LPS2hr <- hclust(as.dist(1-cor(t(diffData.LPS2hr.AVG), method="pearson")), method="complete") 
hr.LPS6hr <- hclust(as.dist(1-cor(t(diffData.LPS6hr.AVG), method="pearson")), method="complete") 
hr.LPSeffect <- hclust(as.dist(1-cor(t(diffData.LPSeffect.AVG), method="pearson")), method="complete") 

#cluster columns using Spearman correlation
hc.untxt <- hclust(as.dist(1-cor(diffData.untxt.AVG, method="spearman")), method="complete") 
hc.LPS2hr <- hclust(as.dist(1-cor(diffData.LPS2hr.AVG, method="spearman")), method="complete") 
hc.LPS6hr <- hclust(as.dist(1-cor(diffData.LPS6hr.AVG, method="spearman")), method="complete") 
hc.LPSeffect <- hclust(as.dist(1-cor(diffData.LPSeffect.AVG, method="spearman")), method="complete") 

# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
mycl.untxt <- cutree(hr.untxt, k=4)
mycl.LPS2hr <- cutree(hr.LPS2hr, k=4)
mycl.LPS6hr  <- cutree(hr.LPS6hr, k=4)
mycl.LPSeffect  <- cutree(hr.LPSeffect, k=4)
mycolhc.untxt <- rainbow(length(unique(mycl.untxt)), start=0.1, end=0.9) 
mycolhc.LPS2hr <- rainbow(length(unique(mycl.LPS2hr)), start=0.1, end=0.9) 
mycolhc.LPS6hr <- rainbow(length(unique(mycl.LPS6hr)), start=0.1, end=0.9) 
mycolhc.LPSeffect <- rainbow(length(unique(mycl.LPSeffect)), start=0.1, end=0.9) 
mycolhc.untxt <- mycolhc.untxt[as.vector(mycl.untxt)] 
mycolhc.LPS2hr <- mycolhc.LPS2hr[as.vector(mycl.LPS2hr)] 
mycolhc.LPS6hr <- mycolhc.LPS6hr[as.vector(mycl.LPS6hr)] 
mycolhc.LPSeffect <- mycolhc.LPSeffect[as.vector(mycl.LPSeffect)] 
#load the gplots package for plotting the heatmap
library(gplots) 
#assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcol <- greenred(75)
library(RColorBrewer)
#plot the hclust results as a heatmap
heatmap.2(diffData.untxt.AVG, Rowv=as.dendrogram(hr.untxt), Colv=NA, col=myheatcol, scale="row", labRow=diffSymbols.untxt, density.info="none", trace="none", RowSideColors=mycolhc.untxt, cexRow=0.75, cexCol=1, margins=c(8,35)) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
heatmap.2(diffData.LPS2hr.AVG, Rowv=as.dendrogram(hr.LPS2hr), Colv=NA, col=myheatcol, scale="row", labRow=NA, density.info="none", trace="none", RowSideColors=mycolhc.LPS2hr, cexRow=0.75, cexCol=1, margins=c(8,35)) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
heatmap.2(diffData.LPS6hr.AVG, Rowv=as.dendrogram(hr.LPS6hr), Colv=NA, col=myheatcol, scale="row", labRow=NA, density.info="none", trace="none", RowSideColors=mycolhc.LPS6hr, cexRow=0.75, cexCol=1, margins=c(8,35)) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
heatmap.2(diffData.LPSeffect.AVG, Rowv=as.dendrogram(hr.LPSeffect), Colv=NA, col=myheatcol, scale="row", labRow=diffSymbols.LPSeffect, density.info="none", trace="none", RowSideColors=mycolhc.LPSeffect, cexRow=0.75, cexCol=1, margins=c(4,30)) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
x11(height=6, width=2); names(mycolhc.untxt) <- names(mycl.untxt); barplot(rep(10, max(mycl.untxt)), col=unique(mycolhc.untxt[hr.untxt$labels[hr.untxt$order]]), horiz=T, names=unique(mycl.untxt[hr.untxt$order])) # Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.
x11(height=6, width=2); names(mycolhc.LPS2hr) <- names(mycl.LPS2hr); barplot(rep(10, max(mycl.LPS2hr)), col=unique(mycolhc.LPS2hr[hr.LPS2hr$labels[hr.LPS2hr$order]]), horiz=T, names=unique(mycl.LPS2hr[hr.LPS2hr$order])) # Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.
x11(height=6, width=2); names(mycolhc.LPS6hr) <- names(mycl.LPS6hr); barplot(rep(10, max(mycl.LPS6hr)), col=unique(mycolhc.LPS6hr[hr.LPS6hr$labels[hr.LPS6hr$order]]), horiz=T, names=unique(mycl.LPS6hr[hr.LPS6hr$order])) # Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.
x11(height=6, width=2); names(mycolhc.LPSeffect) <- names(mycl.LPS6hr); barplot(rep(10, max(mycl.LPS6hr)), col=unique(mycolhc.LPS6hr[hr.LPS6hr$labels[hr.LPS6hr$order]]), horiz=T, names=unique(mycl.LPS6hr[hr.LPS6hr$order])) # Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.

###############################################################################################
#select sub-clusters from 2HR TIMEPOINT
###############################################################################################
clid <- c(2); ysub <- diffData.LPS2hr.AVG[names(mycl.LPS2hr[mycl.LPS2hr%in%clid]),]; hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="complete") 
#NAMES FOR THE SUBCLUSTER ARE NOT PRINTING OUT CORRECTLY....NO STAT5a
x11(); heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=NA, col=myheatcol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc.LPS2hr[mycl.LPS2hr%in%clid]) # Create heatmap for chosen sub-cluster.
#print out row labels in same order as shown in the heatmap
clusterIDs <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterIDs <- as.vector(t(clusterIDs))
#retrieve gene symbols and entrezIDs for selected cluster and print out to an excel spreadsheet for downstream applications (i.e. GO enrichment in DAVID)
myCluster <- cbind(getSYMBOL(clusterIDs, "lumiMouseAll.db"), getEG(clusterIDs, "lumiMouseAll.db"))
write.table(myCluster, "Cluster2.green_2hr_casp8dep.xls", sep="\t", quote=FALSE)

###############################################################################################
#select sub-clusters from 6HR TIMEPOINT
###############################################################################################
clid <- c(3); ysub <- diffData.LPS6hr.AVG[names(mycl.LPS6hr[mycl.LPS6hr%in%clid]),]; hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="complete") 
x11(); heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=NA, col=myheatcol, scale="row", labRow=diffSymbols.LPS6hr, density.info="none", trace="none", RowSideColors=mycolhc.LPS6hr[mycl.LPS6hr%in%clid]) # Create heatmap for chosen sub-cluster.
#print out row labels in same order as shown in the heatmap
clusterIDs <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterIDs <- as.vector(t(clusterIDs))
#retrieve gene symbols and entrezIDs for selected cluster and print out to an excel spreadsheet for downstream applications (i.e. GO enrichment in DAVID)
myCluster <- cbind(getSYMBOL(clusterIDs, "lumiMouseAll.db"), getEG(clusterIDs, "lumiMouseAll.db"))
write.table(myCluster, "Cluster3.blue_LPS6hr_casp8dep.xls", sep="\t", quote=FALSE)
##########################################################################
#Prepare data for running GSEA analysis
#need total data set (not differentially expressed genes) with only entrez IDs as the row names
##########################################################################
#first, take filtered data matrix and make entrez IDs as the rownames
rownames(filtered.matrix) <- myEntrezAll
head(filtered.matrix)
#now go to the GSEA script and use the filtered.matrix object you just created