---
layout: page
---

### Introduction
This is the GitHub repository for a graduate level workshop in transcriptomics. Goals of this workshop include:

* learn to use the R programming language to manage and analyze data
* understand basic statistical approaches used in the analysis of gene expression data
* learn best practices for reproducible data analysis

***Course description.***  The goal of this workshop is to provide a hands-on training environment in which PenVet students and postdocs develop the skills necessary to analyze, summarize and present the results of transcriptomic experiments.  Trainees will work through each module using their own microarray or RNAseq data.  They will learn to import their data into the R/bioconductor environment, process and apply quality control metrics, explore their data using multivariate statistical approaches, and identify differentially expressed genes.  In the second part of the workshop, lessons will focus on interpreting, distilling and presenting results.  Various clustering and visualization representation methods will be used to identify co-regulated modules of genes.  Students will learn to test association of these modules with metadata, and carry out functional enrichment analyses using Gene Ontology (GO) and Gene Set Enrichment Analysis (GSEA).  Students will learn how to summarize the results of GO and GSEA, to compare these results across treatments and species, and to design graphics that present this information in a clear and concise manner.  If time permits, lessons on how to mine and query large public databases of genomic data could also be included. The workshop will conclude with each student producing a comprehensive analysis report in R markdown, and presenting the results and interpretation of their results to the rest of the class.  

***Why is this workshop needed?***  As access to high throughput technology increases, the bottle-neck in biomedical research has shifted from generating data, to analyzing data and integrating diverse data types.  Addressing these needs requires that students and postdocs equip themselves with powerful tools for data mining, manipulation and interrogation.  This workshop focuses specifically on studying global gene expression (transcriptomics) through the use of the R programming environment and the bioconductor suite of bioinformatics packages -- a versatile and robust collection tools for statistics and graphics.  While there are a number of excellent courses in computational biology and statistics offered at UPenn, the material covered is often too general, focuses more on theory or principles, and/or requires a foundation in statistics or computer programming that many biologists lack.  In addition, short workshops have been held in the past, primarily in association with the School of Medicine, but these tend to be one to two days in length and are heavily subscribed. This workshop is designed specifically for PennVet labs that have transcriptomics data in-hand (B.Y.O.D.).  

***Course structure.***  This workshop will consist of one class/week for two hours.  Each class will be a mix of lecture and hands-on exercises, so please bring your laptop and a power supply to each class.  __All classes will start promptly at 3pm and will be held in the Electronic Classroom (room H234) in the PennVet Library.__

***Office hours.***  I will have open office hours every Monday from 10am - 3pm.  __Sign-up for 45min slots using.__
----


### Preparing for the workshop

Everyone is expected to bring their own internet-enabled laptop equiped with either the recent Mac or Windows OS.

Before the first class, each of you will need to attempt to complete the following steps.<br/>
Please email me if you run into problems with any of these tasks.<br/>
All software listed below is free and open-source

* Download and install the appropriate version of the [R Programming Language](http://lib.stat.cmu.edu/R/CRAN/) for your operating system
* Download and install the graphical user interface for R, called [RStudio](http://www.rstudio.com/products/rstudio/download/)
* Download a powerfull text editor. I use [TextWrangler](http://www.barebones.com/products/textwrangler/) (Mac only) and [Sublime](http://www.sublimetext.com/) (Mac and PC)
* Download and install the version control system, called [Git](http://git-scm.com/downloads)
* Sign-up for a free account on [GitHub](https://github.com/), and email me your username.
* If you have a Mac, download the latest version of [XCode](https://developer.apple.com/xcode/) developer tools
* Download the Java-based program for running [GSEA](http://www.broadinstitute.org/gsea/index.jsp). This requires that you sign-up for a free account. You may also need to update to the latest version of [Java](https://www.java.com/en/) for this application to run properly. 
* Download the network analysis platform, [Cytoscape](http://www.cytoscape.org/)


----


### Data Assignments

Each student will be analyzing their own datasets related to their Ph.D. or Postdoc project.  See below for assignment:

Student	|	Lab | Data	|	Topic	|
:------:|:---------:|:-----------:|-----------:
Arindam Basu	|	Atchison	| Array	|	Role of YY1 in B cell fate
Corbett Berry	|	Freedman 	| Array |	Role of CRAC channels in T cell biology
Gautami Das	|	Aguire	|	RNA-seq	|	Genetic basis for retinal cone dysfunction in the dog
Gretchen Harms	|	Hunter	| RNA-seq	|	Tbet-dependent effects of TCR signaling in T cells
Yulyia Katlinskaya	|	Fuchs	| Array |	Role of Type I interferon signaling in the colon
Kelly McCorkell |	May | Array	|	non-canonical NFkB signaling
Katie Morrison  |	Bale  | Array	|	Role of pre-pubertal stress on the aging brain
Fernanda Novais	|	Scott | Array	|	Healing vs. non-healing L. braziliensis skin lesions
Naomi Philip	|	Brodsky | RNA-seq	|	Regulation of TLR signaling by Casp8
Sarah Sneed	|	Povelones | RNA-seq	|	Mosquito transcriptional response to hormone
Raghavi Sudharsan	|	Beltran | RNA-seq	|	Biomarkers of canine retinal degeneration
Camille Syrett  |	Anguera |	Array	| 	LincRNA regulation of HSC gene expression
Brain Gaudette	|	Allman	|	Array	|	central vs peripheral plasma cell programs	
Dan Beiting	|	Beiting	|	Array	|	Innate immunity and phagocytosis in B1 cells


----


### Final Project

Your final project will consist of generating a comprehensive analysis report in R Markdown, and submitting this document as a pull request to the following GitHub repository:

* [Final Project](https://github.com/transcriptomicsworkshop/finalProject): Due 7/1/15


----


### Lectures

As lectures are given, I will include links to both Apple Keynote and PDF files


----


Class	|	Date	|	Topic	|	Reading	|
:------:|:---------:|-----------|:---------:|
1	|	4/1	|	__Introduction to transcriptomics data and platforms__ ([keynote](lectures/Lecture1_TechIntro.key),[pdf](lectures/Lecture1_TechIntro.pdf))
2	|	4/8	|	__Introduction to [R](http://www.r-project.org/), [RStudio](http://www.rstudio.com/), and [Git](http://git-scm.com/)__ ([keynote](lectures/Lecture2_ToolsIntro.key),[pdf](lectures/Lecture2_ToolsIntro.pdf)) | [1](https://training.github.com/kit/downloads/github-git-cheat-sheet.pdf)
3	|	4/15	|	__Version control with Git, [GitHub](https://github.com/), and RStudio__ ([keynote](lectures/Lecture3_Git.key),[pdf](lectures/Lecture3_Git.pdf))	| [1](http://TranscriptomicsWorkshop.github.io/lectures/chocolateMuffin.txt), [2](http://TranscriptomicsWorkshop.github.io/lectures/Git_instructions.pdf)
4	|	4/22	|	__Preprocessing and normalizing microarray and RNAseq data__ ([keynote](lectures/Lecture4_preprocessing.key),[pdf](lectures/Lecture4_preprocessing.pdf))
5	|	4/29	|	__Exploratory analysis of expression data (gene-agnostic)__	| [1](http://shinyapps.stat.ubc.ca/r-graph-catalog/)
	|			|	* _defining variables in a study design file_
	|			|	* _viewing sample relationships with dist and hclust functions_ 
	|			|	* _using the prcomp function to carry out a principal component analysis_
	|			|	* _using [ggplot2](http://ggplot2.org/) to visualize your PCA analysis_
	|			|	* _using PCA 'loadings' to examine the relationship between samples and principal components_
6	|	5/6	|	__Visual analysis of gene expression data__| [1](http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)
	|		|	* _using [dplyr's](http://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html) mutate, sort and filter commands to edit data tables_
	|		|	* _static and interactive scatter plots using ggplot2 and [ggvis](http://ggvis.rstudio.com/)_
7 | 5/20 | __Identifying differentially expressed genes__
  |		| * _using the [Limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html) package_ | [1](http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf), [2](http://www.statsci.org/smyth/pubs/VoomPreprint.pdf)
  |		| * _P values, false discovery rates, and volcano plots_
8 | 5/27  | __Visualizing and dissecting your differentially expressed gene list__
9 | 6/10  | __Understanding and leveraging [Gene Ontology](http://geneontology.org/) for function enrichment analyses__ | [1](http://david.abcc.ncifcrf.gov/), [2](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2615629/)
10 | 6/16?  | __Gene Set Enrichment Analysis ([GSEA](http://www.broadinstitute.org/gsea/index.jsp)) and gene [signatures](http://www.broadinstitute.org/gsea/msigdb/index.jsp)__ | [1](http://mootha.med.harvard.edu/PubPDFs/Subramanian2005.pdf), [2](http://bioinf.wehi.edu.au/software/MSigDB/)
11  | TBD | __Beyond gene lists, to network topology using [Cytoscape](http://www.cytoscape.org/)__ | [1](http://healthsciences.ucsd.edu/som/medicine/research/labs/ideker/publications/Documents/saito_NatMethod_2012.pdf)
12  | TBD  | __Deploying results as interactive graphics using D3.j and [Shiny](http://shiny.rstudio.com/)__ | [1](http://shiny.rstudio.com/articles/cheatsheet.html)
13  | TBD  | __Making your anaysis pipeline transparent and reproducible using [R Markdown](http://rmarkdown.rstudio.com/) and [Knitr](http://yihui.name/knitr/)__ | [1](http://rmarkdown.rstudio.com/RMarkdownCheatSheet.pdf)
14  | TBD | __Final analysis reports due (in Markdown)__


----

### R Scripts

The lectures will be paired with R scripts which will 'step' through the process of analyzing transcriptional profiling data from microarray and RNAseq experiments

----

Script	|	Name	|	Relevant R packages (and functions)	|
:------:|---------|:-----------:|
1	|	[Step1_preprocessAffy.R](R_scripts/Step1_preprocessAffy.R) 	|	Oligo, annotate, platform-specific annotation database
2	|	[Step1_preprocessLumi.R](R_scripts/Step1_preprocessLumi.R)	|	Lumi, annotate, platform-specific annotation database
3	|	[Step1_preprocessRNAseq.R](R_scripts/Step1_preprocessRNAseq.R) 	|	ShortRead (buildIndex, align, featureCounts), edgeR, Limma, biomaRt
4	|	[Step2_dataExploration_part1.R](R_scripts/Step2_dataExploration_part1.R) 	|	base R (dist, hclust, prcomp), ggplot2, reshape2
5	|	[Step3_dataExploration_part3.R](Step3_dataExploration_part3.R)  	|	dplyr (filter, arrange, select), ggplot2, ggviz
6	|	[Step4_diffGenes.R](Step4_diffGenes.R)  	|	Limma (topTable and DecideTests),
7	|	[Step5_heatmap.R](Step5_heatmap.R) 	|	heatmap.2


----


### Additional Optional Reading

* **Introductory Statistics with R** by Peter Dalgaard [ [pdf](http://www.academia.dk/BiologiskAntropologi/Epidemiologi/PDF/Introductory_Statistics_with_R__2nd_ed.pdf), [Amazon](http://www.amazon.com/Introductory-Statistics-R-Computing/dp/0387954759) ]  
* **The Art of R programming** by Norman Matloff [ [pdf](http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CCAQFjAA&url=http%3A%2F%2Fsens.tistory.com%2Fattachment%2Fcfile8.uf%402375DC3D515423F9110CA1.pdf&ei=E-8FVO6dAYmnggSttoD4Bg&usg=AFQjCNE1UmWRG3i9ugNDSXN2WjRSTkkUjA&sig2=U958L8LG42vuhHdPKKBHHw&bvm=bv.74115972,d.eXY), [Amazon](http://www.amazon.com/Art-Programming-Statistical-Software-Design/dp/1593273843/ref=sr_1_1?s=books&ie=UTF8&qid=1409674972&sr=1-1&keywords=the+art+of+r+programming) ]  
* **Bioinformatics data skills** by Vincent Buffalo [ pdf, [O'Reilly](http://shop.oreilly.com/product/0636920030157.do) ]  
* **ProGit** by Scott Chacon [ [pdf](http://git-scm.com/book), [Amazon](http://www.amazon.com/Pro-Git-Scott-Chacon/dp/1430218339) ]  
