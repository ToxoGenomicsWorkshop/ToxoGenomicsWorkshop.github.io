---
layout: page
---

### Introduction
This is the GitHub repository for the Toxoplasma Conference RNAseq Workshop. Goals of this workshop include:

* gain an understanding of how RNAseq data is generated and how the raw data is processed (for practical reasons, we will dicuss but not demonstrate this)
* learn to use the R programming language to manage and analyze data (hands-on with example dataset)

***Rationale for this Workshop***  Students and postdocs training in the biomedical sciences today are faced with the challenge  analyzing and interpreting large-scale datasets, as well data being continuously made in the public domain.  Research in Toxoplasma has mine large data sets.  Once a problem specific to genomics, the issue of 'big data' is becoming ubiquitous, as everything from microscopy to flow cytometry to gene expression profiling produce  This workshop provides students and postdocs with an opportunity to begin developing the skills necessary to analyze, summarize and present the results of RNAseq experiments.  The workshop will be run as a hands-on tutorial using a real data-set  


----


### Workshop Schedule

Time	|	Description	|
:------:|---------|
8:30-9:00am	|	Meet in College lecture room for introductions and overview
9:00-10:30am	|	Analyzing RNAseq data - Part I
10:30-11:00am	|	Break
11:00-12:30am	|	Analyzing RNAseq data - Part II
12:30-2:00pm	|	Lunch with discussion of alternative career paths
2:00-5:00pm	|	Understanding RNAseq data in the context of the genome browser
5:00-7:00pm	|	Cocktails, lite bites, etc


----


### Preparing for the workshop

**There are really three main things you will need for the workshop**

***1. Your laptop computer***<br/>
Everyone is encouraged to bring their own internet-enabled laptop equiped with either the recent Mac or Windows OS.

***2. A "toolbox" of free software***<br/>
* Download and install the appropriate version of the [R Programming Language](http://lib.stat.cmu.edu/R/CRAN/) for your operating system
* Download and install the graphical user interface for R, called [RStudio](http://www.rstudio.com/products/rstudio/download/)
* Download a text editor. I use [TextWrangler](http://www.barebones.com/products/textwrangler/) (Mac only) and [Sublime](http://www.sublimetext.com/) (Mac and PC)
* If you have a Mac, download the latest version of [XCode](https://developer.apple.com/xcode/) developer tools


***3. some sample data***<br/>
* We will be analyzing RNAseq data from this [2014 BMC Genomics paper](materials/Pittman_BMCgenomics_TgBrady.pdf) published by [Laura Knoll's lab](http://www.medmicro.wisc.edu/people_faculty_profile.php?id=ljknoll&view=intro).  
* The original raw data is available on the Short Read Archive [here](http://www.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1156954), but you __DO NOT__ need to download this.  
* Please download the 'digital gene expression list' (DGEList) from [here](materials/DGEList)
* In addition to the data, please download this simple [text file]() that describes the design of the study

***4. An analysis script*** - the R script we'll use in the workshop can be downloaded [Here]().

-----

### After the workshop

***if you find this workshop interesting, you may want to consider downloading and learning about the following tools***

* the network analysis platform, [Cytoscape](http://www.cytoscape.org/)
* A version control system, such as SVN or [Git](http://git-scm.com/downloads)
* Sign-up for a free account on [GitHub](https://github.com/), and email me your username.
* Download the Java-based program for running [GSEA](http://www.broadinstitute.org/gsea/index.jsp). This requires that you sign-up for a free account. You may also need to update to the latest version of [Java](https://www.java.com/en/) for this application to run properly. 


----


### Lecture Slides and Code

A PDF copy of the lecture slides can be downloaded [here](materials/ToxoRNAseqWorkshop.pdf).

----


### Optional Reading

* **Introductory Statistics with R** by Peter Dalgaard [ [pdf](http://www.academia.dk/BiologiskAntropologi/Epidemiologi/PDF/Introductory_Statistics_with_R__2nd_ed.pdf), [Amazon](http://www.amazon.com/Introductory-Statistics-R-Computing/dp/0387954759) ]  
* **The Art of R programming** by Norman Matloff [ [pdf](http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CCAQFjAA&url=http%3A%2F%2Fsens.tistory.com%2Fattachment%2Fcfile8.uf%402375DC3D515423F9110CA1.pdf&ei=E-8FVO6dAYmnggSttoD4Bg&usg=AFQjCNE1UmWRG3i9ugNDSXN2WjRSTkkUjA&sig2=U958L8LG42vuhHdPKKBHHw&bvm=bv.74115972,d.eXY), [Amazon](http://www.amazon.com/Art-Programming-Statistical-Software-Design/dp/1593273843/ref=sr_1_1?s=books&ie=UTF8&qid=1409674972&sr=1-1&keywords=the+art+of+r+programming) ]  
* **Bioinformatics data skills** by Vincent Buffalo [ pdf, [O'Reilly](http://shop.oreilly.com/product/0636920030157.do) ]  
* **ProGit** by Scott Chacon [ [pdf](http://git-scm.com/book), [Amazon](http://www.amazon.com/Pro-Git-Scott-Chacon/dp/1430218339) ]  
