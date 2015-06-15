---
layout: page
---

### Introduction
This is the GitHub repository for the Toxoplasma Pre-Conference Workshop on data analysis and carrer perspectives. Goals of this workshop include:

* to interact with scientists who have pursued careers outside of the traditional academic tenure-track 
* to gain an understanding of how RNAseq data is generated and how the raw data is processed (for practical reasons, we will discuss but not demonstrate this)
* to learn to use the R programming language to manage and explore the data (hands-on with example dataset)
* to learn to use ToxoDB to interpret RNAseq data in the context of a genome browser (also hands-on)

The workshop will have two main components: 1) career perspectives from scientists outside of the traditional tenure-track.  2) hands-on workshop on analysing RNAseq data.  The latter will be broken a morning and afternoon workshop. The morning workshop will cover RNAseq analysis in R/bioconductor (Dan Beiting); and the afternoon will focus on using the genome browser to view RNAseq data (Omar Harb).  We encourage all attendees to bring their laptops to allow participation in the hands-on exercises. 

***Rationale for this Workshop***  Students and postdocs training in the biomedical sciences today are faced with the challenge of analyzing their own large-scale datasets, as well data being continuously made available in the public domain.  The ability to carry out forward and reverse genetic screens in _Toxoplasma_, together with the excellent genome resources for this pathogen, makes _Toxoplasma_ research particularly adept at genome-scale experiments in both the parasite and the host.  This workshop serves an a hands-on introduction to the methods and challenges assocated analzing and visualizing RNAseq data.  The workshop will be run as a hands-on tutorial using a real data-set  


----


### Workshop Schedule

**Tuesday** - Gettysburg Hotel, Lincoln Bar & Atrium

Time	|	Description	|
:------:|:--------
7:00pm	|	Informal introductions & discussions (food and drink courtesy of the Burroughs Wellcome Fund)<br/>


**Wednesday** - Gettysburg College, Science Center 200 

Time	|	Description	|
:------:|:--------
8:30-9:00am	|	Meet in College lecture room for data analysis workshop - introductions and overview
9:00-10:30am	|	Introduction to RNAseq technology, data and tools
10:30-10:45am	|	Short break
10:45-12:15am	|	Using R/bioconductor to mine RNAseq data 
12:15-2:15pm	|	Lunch with panel discussion of alternative reasearch careers
2:15-3:45pm	|	Interpreting RNAseq in the genome browser - Part I
3:45-4:00pm	|	Short break
4:00-5:30pm |	Interpreting RNAseq in the genome browser - Part II
5:30-7:00pm	|	Cocktails, lite bites, etc


----


### Panel discussion of alternative careers in basic research (over lunch)

Name	|	affiliation	|	topic	|
:--------|:--------|:--------
Eric Villegas	|	EPA	|	Research opportunities in the Federal Agency
Robert Donald	|	Pfizer	|	Commercial R&D careers
Omar Harb	|	UPenn	|	Alternative academic careers
Dan Beiting	|	UPenn	|	Research-track academic careers


----


### Preparing for the data analysis workshop

**There are four main things you will need in order to particpate in the hands-on aspects of the workshop**

**1. Your laptop computer**<br/> 

* Everyone is encouraged to bring their own internet-enabled laptop equiped with either the recent Mac or Windows OS.

**2. A "toolbox" of free software**<br/>

* Download and install the appropriate version of the [R Programming Language](http://lib.stat.cmu.edu/R/CRAN/) for your operating system
* Download and install the graphical user interface for R, called [RStudio](http://www.rstudio.com/products/rstudio/download/)
* Download a text editor. I use [TextWrangler](http://www.barebones.com/products/textwrangler/) (Mac only) and [Sublime](http://www.sublimetext.com/) (Mac and PC)


**3. Sample dataset**<br/>

* We will be analyzing RNAseq data from this [2014 PLOS One](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0111297) published by [Kami Kim's lab](http://www.einstein.yu.edu/faculty/4972/kami-kim/).  
* The original raw data is available on the Short Read Archive [here](http://www.ebi.ac.uk/ena/data/view/SRP045423), but you __DO NOT__ need to download this.  Instead, please download the 'digital gene expression list' (DGEList) from [here](materials/DGEList)
* In addition to the data, please download this simple [text file](materials/studyDesign.txt) that describes the design of the study

**4. An analysis script**<br/>

* the R script we'll use in the workshop can be downloaded [here](materials/Toxo_RNAseq_analysis.R).

-----

### After the workshop

**Although it is beyond the scope of the current workshop, there are a number of other tools -- essential in any data scientists repertoire -- that you should consider downloading and learning more about.**

* the network analysis platform, [Cytoscape](http://www.cytoscape.org/)
* A version control system such as [Git](http://git-scm.com/downloads)
* A free account on [GitHub](https://github.com/)
* The Java-based program for running [GSEA](http://www.broadinstitute.org/gsea/index.jsp). This requires that you sign-up for a free account. You may also need to update to the latest version of [Java](https://www.java.com/en/) for this application to run properly. 


----


### Lecture Slides

After the workshop is over, you'll be able to access the lecture slides for the [panel discussion of alternative careers](), as well as for the data analysis workshop: Dan's slides [here](materials/ToxoRNAseqWorkshop.pdf) and Omar's portion [here](materials/ToxoRNAseqWorkshop.pdf).

----


### Optional Reading for data analysis 

* **Introductory Statistics with R** by Peter Dalgaard [ [pdf](http://www.academia.dk/BiologiskAntropologi/Epidemiologi/PDF/Introductory_Statistics_with_R__2nd_ed.pdf), [Amazon](http://www.amazon.com/Introductory-Statistics-R-Computing/dp/0387954759) ]  
* **The Art of R programming** by Norman Matloff [ [pdf](http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CCAQFjAA&url=http%3A%2F%2Fsens.tistory.com%2Fattachment%2Fcfile8.uf%402375DC3D515423F9110CA1.pdf&ei=E-8FVO6dAYmnggSttoD4Bg&usg=AFQjCNE1UmWRG3i9ugNDSXN2WjRSTkkUjA&sig2=U958L8LG42vuhHdPKKBHHw&bvm=bv.74115972,d.eXY), [Amazon](http://www.amazon.com/Art-Programming-Statistical-Software-Design/dp/1593273843/ref=sr_1_1?s=books&ie=UTF8&qid=1409674972&sr=1-1&keywords=the+art+of+r+programming) ]  
* **Bioinformatics data skills** by Vincent Buffalo [ pdf, [O'Reilly](http://shop.oreilly.com/product/0636920030157.do) ]  
* **ProGit** by Scott Chacon [ [pdf](http://git-scm.com/book), [Amazon](http://www.amazon.com/Pro-Git-Scott-Chacon/dp/1430218339) ]  
