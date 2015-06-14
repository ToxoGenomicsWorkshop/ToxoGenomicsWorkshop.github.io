---
layout: page
---

### Introduction
This is the GitHub repository for the Toxoplasma Conference RNAseq Workshop. Goals of this workshop include:

* gain an understanding of how RNAseq data is generated and how the raw data is processed (for practical reasons, we will discuss but not demonstrate this)
* learn to use the R programming language to manage and explore the data (hands-on with example dataset)
* use ToxoDB to gain an understanding of how to 'read' RNAseq data in a genome browser (also hands-on)

***Rationale for this Workshop***  Students and postdocs training in the biomedical sciences today are faced with the challenge  analyzing and interpreting large-scale datasets, as well data being continuously made available in the public domain.  Research in Toxoplasma has mine large data sets.  Once a problem specific to genomics, the issue of 'big data' is becoming ubiquitous, as everything from microscopy to flow cytometry to gene expression profiling produce  This workshop provides students and postdocs with an opportunity to begin developing the skills necessary to analyze, summarize and present the results of RNAseq experiments.  The workshop will be run as a hands-on tutorial using a real data-set  


----


### Workshop Schedule

Time	|	Description	|
:------:|---------|
8:30-9:00am	|	Meet in College lecture room for introductions and overview
9:00-10:30am	|	Introduction to RNAseq technology, data and tools
10:30-10:45am	|	Short break
10:45-12:15am	|	Using R/bioconductor to mine RNAseq data 
12:15-2:15pm	|	Lunch with discussion of alternative career paths
2:15-3:45pm	|	Interpreting RNAseq in the genome browser - Part I
3:45-4:00pm	|	Short break
4:00-5:30pm |	Interpreting RNAseq in the genome browser - Part II
5:30-7:00pm	|	Cocktails, lite bites, etc


----


### Preparing for the workshop

The workshop will have two basic parts: 1) analysing RNAseq data in R/bioconductor; and 2) using the genome browser to view RNAseq data.  We encourage all attendees to bring their laptops to allow participation in the hands-on exercises.  **There are three main things you will need to do to particpate in the hands-on aspects of the workshop**

***1. Your laptop computer***<br/> - Everyone is encouraged to bring their own internet-enabled laptop equiped with either the recent Mac or Windows OS.

***2. A "toolbox" of free software***<br/>

* Download and install the appropriate version of the [R Programming Language](http://lib.stat.cmu.edu/R/CRAN/) for your operating system
* Download and install the graphical user interface for R, called [RStudio](http://www.rstudio.com/products/rstudio/download/)
* Download a text editor. I use [TextWrangler](http://www.barebones.com/products/textwrangler/) (Mac only) and [Sublime](http://www.sublimetext.com/) (Mac and PC)


***3. Sample dataset***<br/>

* We will be analyzing RNAseq data from this [2014 PLOS One](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0111297) published by [Kami Kim's lab](http://www.einstein.yu.edu/faculty/4972/kami-kim/).  
* The original raw data is available on the Short Read Archive [here](http://www.ebi.ac.uk/ena/data/view/SRP045423), but you __DO NOT__ need to download this.  Instead, please download the 'digital gene expression list' (DGEList) from [here](materials/DGEList)
* In addition to the data, please download this simple [text file](materials/studyDesign.txt) that describes the design of the study

***4. An analysis script*** - the R script we'll use in the workshop can be downloaded [Here](materials/Toxo_RNAseq_analysis.R).

-----

### After the workshop

***if you find this workshop interesting, you may want to consider downloading and learning about the following tools***

* the network analysis platform, [Cytoscape](http://www.cytoscape.org/)
* A version control system, such as SVN or [Git](http://git-scm.com/downloads)
* A free account on [GitHub](https://github.com/)
* The Java-based program for running [GSEA](http://www.broadinstitute.org/gsea/index.jsp). This requires that you sign-up for a free account. You may also need to update to the latest version of [Java](https://www.java.com/en/) for this application to run properly. 


----


### Lecture Slides

After the workshop is over, you'll be able to access the lecture slides [here](materials/ToxoRNAseqWorkshop.pdf).

----


### Optional Reading

* **Introductory Statistics with R** by Peter Dalgaard [ [pdf](http://www.academia.dk/BiologiskAntropologi/Epidemiologi/PDF/Introductory_Statistics_with_R__2nd_ed.pdf), [Amazon](http://www.amazon.com/Introductory-Statistics-R-Computing/dp/0387954759) ]  
* **The Art of R programming** by Norman Matloff [ [pdf](http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CCAQFjAA&url=http%3A%2F%2Fsens.tistory.com%2Fattachment%2Fcfile8.uf%402375DC3D515423F9110CA1.pdf&ei=E-8FVO6dAYmnggSttoD4Bg&usg=AFQjCNE1UmWRG3i9ugNDSXN2WjRSTkkUjA&sig2=U958L8LG42vuhHdPKKBHHw&bvm=bv.74115972,d.eXY), [Amazon](http://www.amazon.com/Art-Programming-Statistical-Software-Design/dp/1593273843/ref=sr_1_1?s=books&ie=UTF8&qid=1409674972&sr=1-1&keywords=the+art+of+r+programming) ]  
* **Bioinformatics data skills** by Vincent Buffalo [ pdf, [O'Reilly](http://shop.oreilly.com/product/0636920030157.do) ]  
* **ProGit** by Scott Chacon [ [pdf](http://git-scm.com/book), [Amazon](http://www.amazon.com/Pro-Git-Scott-Chacon/dp/1430218339) ]  
