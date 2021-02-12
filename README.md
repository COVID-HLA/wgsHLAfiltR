# Whole Genome Sequence HLA Filter (wgsHLAfiltR)  
v0.3.1 (February 11, 2021)

***
The *wgsHLAfiltR* R package extracts reads that map to classical HLA loci (HLA-A, -C, -B, -DRB1, -DRB3/4/5, -DQA1, -DQB1, -DPA1 and -DPB1) from paired or individual whole-genome or whole-exome sequencing (WGS/WES) fastg.gz files, and writes a new set of fastq.gz files that contain reads that map to the classical HLA loci (using hg19 coordinates). 

The package was developed for use with the [HLA|COVID-19 Database](https://database-hlacovid19.org)'s [Omixon CLI Explore Genotyping Portal](https://database-hlacovid19.org/shiny/Omixon-Genotyping-Portal/), with the aim of minimizing FASTQ upload time and read processing time.

***
### Package Requirements
*wgsHLAfiltR* (version 0.3.1) currently functions only on Unix, Linux or MacOS systems. 

*wgsHLAfiltR* requires [R v4.0.0 or higher](https://cran.r-project.org) to run. 

[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (versions 2.3 - 2.4) must be installed on the system running R in order for *wgsHLAfiltR* to function.

***
### Package Installation
To install the *wgsHLAfiltR* package in the R environment, first install the [*devtools*](https://CRAN.R-project.org/package=devtools) package using the command *install.packages("devtools")* in the R console.

    > install.packages("devtools")

With *devtools* installed, install *wgsHLAfiltR* using the command *devtools::install_github(repo="COVID-HLA/wgsHLAfiltR/wgsHLAfiltRpackage",ref="main")* in the R console.

    > devtools::install_github(repo="COVID-HLA/wgsHLAfiltR/wgsHLAfiltRpackage",ref="main")

Please note that the package includes > 130MB of reference alignment data, which may result in longer than expected installation times.

***
### Using the Package
While the package includes several functions, *filterHLA()* is the main function for filtering HLA reads. This function takes two arguments, *inputDirectory*, which identifies the path to the directory containing the WGS/WES FASTQ files, and *outputDirectory*, which specifies the directory into which the HLA-only FASTQ files should be written. 

    > readFilteringData <- filterHLA(inputDirectory=inputDir,outputDirectory=outputDir)

A value for *inputDirectory* is required, but a value for *outputDirectory* is optional. If *outputDirectory* is not specified, HLA-only FASTQ files will be written into a directory named "Results" in the R working directory. If the "Results" directory is not present in the R working directory, one will be created.

Note that Bowtie 2 requires file path and file names that do not include whitespaces. See the [Bowtie 2 Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output) for additional details.

<br>
 
#### Package Details
The *filterHLA()* function extracts a "name" from each FASTQ file based on the position of the first underscore in the FASTQ filename. 

For example, if a pair of fastq.gz files are named "ABC123-45DE_67FG_R1.fastq.gz" and "ABC123-45DE_67FG_R2.fastq.gz", *filterHLA()* will identify the name for that pair of files as "ABC123-45DE", and the resulting HLA-only fastq.gz files generated will be named "ABC123-45DE_HLA_1.fastq.gz" and "ABC123-45DE_HLA_2.fastq.gz".

The *filterHLA()* function returns a list object of named list elements for each FASTQ "name" that was processed. Each named list identifies parameters and file paths used in the read extraction process.

***
### Questions
For additional information about the *wgsHLAfiltR* package, the [HLA|COVID-19 Database](https://database-hlacovid19.org) or the [COVID-19|HLA & Immunogenetics Consortium](http://www.hlacovid19.org/about/), email <covid.hla@gmail.com>.