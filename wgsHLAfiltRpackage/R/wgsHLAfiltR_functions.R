## Main wrapper for extracting reads -- v 1.3.0 03/22/2021 SJM

#' Main Function For Extracting Classical HLA Reads From fastq.gz Files
#' 
#' This function calls all of the wgsHLAfiltR package's secondary functions, applying them to extract HLA reads from fastq.gz files.
#' 
#' This function identifies whole-genome and whole-exome reads that map to the classical HLA loci (HLA-A, -C, -B, -DRB1, -DRB3/4/5, -DQA1, -DQB1, -DPA1, and -DPB1) in fastq.gz files in a user-specified directory, and writes fastq.gz files containing only classical HLA-specific reads in a second output directory. Reads for the HLA genes are mapped using IPD-IMGT/HLA reference sequences.
#' @param inputDirectory A path to the directory that contains the fastq.gz files to be filtered. This parameter is required. 
#' @param outputDirectory A path to an existing directory where the HLA-specific fastq.gz files will be written. If no directory is specified, a 'Results' directory will be created in the working directory. If a 'Results' directory is already present in the working directory, output files will be written in that directory.
#' @keywords filter HLA fastq reads
#' @return A list object describing all of the pre-and post-filtering reads for each subject. 
#' @note This function requires that bowtie2 be installed in the local environment. This version only functions on unix/linux/macOS installations of R.
#' @examples 
#' # Filter non-HLA reads in fastq.gz files in a single directory
#' # inputDir <- paste(getwd(),"unfilteredReads",sep="/")
#' # outputDir <- paste(getwd(),"HLAfilteredReads",sep="/")
#' # readFilteringData <- filterHLA(inputDirectory=inputDir,outputDirectory=outputDir) 
#' @export

filterHLA <- function(inputDirectory,outputDirectory="Results"){
  if(.Platform$OS.type == "windows") {return(cat("Currently, wgsHLAFiltR only functions on Unix, Linux, and macOS systems.\n\nSupport for Windows systems is pending.","\n"))}
  
  if(missing(inputDirectory)) {return(cat("Please use the inputDirectory parameter to identify a directory that contains the fastq.gz files to be filtered.","\n"))}
  
  if(!dir.exists(outputDirectory)) {dir.create(outputDirectory,recursive = TRUE)}
  
  rawFastqDirectory <- inputDirectory # can be set to raw sequence or extractedFastq directory
  fastqPattern <- 'fastq' # use '_KIR_' to find already extracted files, otherwise use 'fastq' or whatever fits your data
  threads <- 50
  resultsDirectory <- outputDirectory # Set the master results directory (all pipeline output will be recorded here)
  shortNameDelim <- '_' # can set a delimiter to shorten sample ID's (ID will be characters before delim)
  minDP <- 10
  
  ## Check to make sure samtools is accessible
  bowtie2 <- system2('which', c('bowtie2'), stdout=T, stderr=T)
  check.system2_output(bowtie2, 'bowtie2 not found. Please install bowtie2.')
  
  # Build up a list of sample objects
  sampleList <- general.paired_sample_objects(rawFastqDirectory, fastqPattern, resultsDirectory, shortNameDelim) # no need to change
  
  # PING2 extractor ---------------------------------------------------------
  cat('\n\n----- Performing HLA extraction -----')
  # Define the extracted fastq directory
  extractedFastqDirectory <- file.path(resultsDirectory)
  # Run PING2 extractor
  sampleList <- extractor.run(sampleList,threads,extractedFastqDirectory,forceRun=T) # set forceRun=T if you want to force alignments
  
  return(sampleList)
}

## This function checks to make sure the output of system2 is valid

#' System Check
#' 
#' Checks the exteral system for the avilability of external resources. If the resource is not available, execution of the package is halted.
#' 
#' For internal wgsHLAfiltR use.
#' @param system2_output The output of a command passed to the system2() function.
#' @param system2_error The message returned if the external resource is unavailable.
#' @note This function is for interal use. 
#' @export

check.system2_output <- function(system2_output, system2_error){
  
  ## Checking if the attributes of system2_command are NULL, if not the command was not found
  if(!is.null(attributes(system2_output))){
    cat('\n',system2_output)
    stop(system2_error, '. Stopping program.')
  }
}


## This function identifies pairs of fastq files.

#' General Paired Sample Objects
#' 
#' Finds pairs of fastq files in inputDirectory for each subject, and returns a list of sample objects that contain the paired fastq file names.
#' 
#' For internal wgsHLAfiltR use.
#' @param rawFastqDirectory The directory containing unfiltered fastq files.
#' @param fastqPattern A sub-string of the fastq file names that identifies the fastq files to be filtered. A value of 'fastq' specifies all fastq files.
#' @param resultsDirectory The directory into which the filtered fastq files should be written. 
#' @param shortNameDelim A delimeter used to shorten fastq filenames.
#' @return A list of sample objects that contain the paired fastq file names for each subject, identifying each fastq pair, the path to each paired file, and a logical indicator of gzip status.
#' @importFrom stringr str_split fixed
#' @importFrom methods setRefClass
#' @note This function is for internal use. 
#' @export

general.paired_sample_objects <- function(rawFastqDirectory,fastqPattern,resultsDirectory,shortNameDelim){
  #####
  ## This function takes in a directory and a file name pattern and attempts to pair fastq files
  ## Returns a list of sample objects that contain the paired fastq file names
  #####
  
  cat("\nAttempting automatic fastq pairing in", rawFastqDirectory, "using", fastqPattern)
  
  ## Find all the files in sampleDirectory that match fastqPattern
  unpairedFastqList <- list.files(path=rawFastqDirectory, pattern=fastqPattern)
  
  ## To pair reads, we will split the file names by fastqPattern, then continuously chop a
  ## character off the end of each name until the number of unique names is exactly half of the total names
  
  ## Setting up an initial fastq list that splits the files names by fastqPattern
  strList <- sapply(unpairedFastqList, function(x) str_split(x, fastqPattern)[[1]][1])
  
  ## Shorten names based on delim
  if(shortNameDelim != ''){
    strList <- sapply(strList, function(x) str_split(x, fixed(shortNameDelim))[[1]][1])
  }
  
  ## Setting the maximum number of times to chop off the last character to the length of the shortest name
  maxChop <- min(sapply(strList, nchar))
  
  ## Iterate from 0 to maxChop, cutting off i characters from the file names each time
  for(i in 0:maxChop){
    
    ## In each iteration, the file names are reset. There is no particular reason I implemented it this way
    subStrList <- strList
    
    ## Cut off i characters from the end of each fastq file name
    subStrList <- sapply(subStrList, function(x) substr(x, 1, nchar(x)-i))
    
    ## After cutting, determine the unique names
    uniqueFastqList <- unique(subStrList)
    
    ## If the number of unique names is exactly half of the total names, then cutting should be finished
    ## and it is time to move on to matching!
    if(length(uniqueFastqList) == (length(subStrList)/2)){
      break
    }
    
    ## Pairing failed if i reaches maxChop. This will raise a fatal error.
    if(i == maxChop){
      stop("Was not able to pair fastq file names, please check that fastqPattern and sampleDirectory are set correctly.")
    }
  }
  
  ## Initialize a list for storing the paired fastq file names
  pairedFastqList <- list()
  
  ## Iterate through the unique fastq names to make pairings
  for(fastqName in uniqueFastqList){
    
    ## Pull out the matches for fastqName in the subStrList
    fastqMatches <- subStrList[fastqName == subStrList]
    
    ## Determine how many matches there are for fastqName in the subStrList
    matchCount <- length(fastqMatches)
    
    ## Stop the program if more or less than 2 matches are found for fastqName
    if(matchCount != 2){
      cat('\n',names(fastqMatches))
      stop('Auto fastq matching failed due to an improper number of matches for ',fastqName)
    }
    
    ## Save the file names in the pairedFastqList under the unique name
    pairedFastqList[[fastqName]] <- names(fastqMatches)
  }
  
  cat("\nFound", length(uniqueFastqList), "samples in", rawFastqDirectory)
  
  ## Initializing a sample object list. This will be returned
  output.sampleList <- list()
  
  for(i in 1:length(pairedFastqList)){
    
    ## Pulling the current working element out of the list
    pairedFastq <- pairedFastqList[i]
    
    ## Creating an absolute path to the first fastq file
    rawfastq1path <- normalizePath(file.path(rawFastqDirectory, pairedFastq[[1]][1]), mustWork=T)
    
    ## Creating a absolute path to the second fastq file
    rawfastq2path <- normalizePath(file.path(rawFastqDirectory, pairedFastq[[1]][2]), mustWork=T)
    
    ## Checking if the first fastq file is gzipped
    gzip <- substr(rawfastq1path, nchar(rawfastq1path)-2, nchar(rawfastq1path)) == '.gz'
    
    ## Building a sample object and adding it to sampleList
    output.sampleList[[names(pairedFastq)]] <- seqSample(name=names(pairedFastq),
                                                      rawfastq1path=rawfastq1path,
                                                      rawfastq2path=rawfastq2path,
                                                      gzip=gzip,
                                                      failed=FALSE)
  }
  
  cat("\nAll samples were successfully paired")
  return(output.sampleList)
}

## Read Alignment

#' Bowtie2 Alignment of Reads
#'
#' Performs bowtie2 alignments of reads for a subject.
#' 
#' For internal wgsHLAfiltR use.
#' @param bowtie2_command The command being sent to bowtie2 in the external system.
#' @param threads The number of threads specified for bowtie2. 
#' @param currentSample A specific element of the sampleList object returned by the general.paired_sample_objects() function.
#' @param extractedFastqDirectory The directory where filtered fastq.gz files are written. 
#' @return The aligned reads for an individual subject.
#' @note This function is for interal use. 
#' @export

extractor.bowtie2_align <- function(bowtie2_command, threads, currentSample, extractedFastqDirectory){
  currentSample$samPath <- file.path(extractedFastqDirectory,paste0(currentSample$name,'.sam'))
  optionsCommand <- c( #'-x Resources/extractor_resources/reference/output',
                      paste0('-x ',paste(find.package("wgsHLAfiltR"),"reference/output",sep="/")),
                      '-5 3','-3 7','-L 20','-i S,1,0.5','--score-min L,0,-0.187',
                      '-I 75','-X 1000', 
                      paste0('-p ',threads),
                      paste0('-1 ',currentSample$rawfastq1path),
                      paste0('-2 ',currentSample$rawfastq2path),
                      paste0('-S ',currentSample$samPath),
                      paste0('--al-conc-gz ',
                             file.path(extractedFastqDirectory,paste0(currentSample$name,'_HLA_%.fastq.gz'))
                      ),
                      '--un delete.me'
  )
  cat('\n\n',bowtie2_command, optionsCommand)
  output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  
  if(!is.null(attributes(output.sampleAlign))){
    cat('\nBowtie2 alignment failed, retrying...')
    output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  }
  
  check.system2_output(output.sampleAlign, 'Bowtie2 HLA extraction alignment failed.')
  
  cat('\n',paste0(output.sampleAlign, collapse='\n'))
  
  currentSample[['HLAfastq1path']] <- normalizePath(file.path(extractedFastqDirectory,paste0(currentSample$name,'_HLA_1.fastq.gz')))
  currentSample[['HLAfastq2path']] <- normalizePath(file.path(extractedFastqDirectory,paste0(currentSample$name,'_HLA_2.fastq.gz')))
  
  cat('\n\nSuccessfully extracted HLA reads for',currentSample$name)
  
  cat('\nCleaning up alignment files.')
  file.remove('delete.me')
  file.remove(currentSample$samPath)
  
  return(currentSample)
}

## Run the extractor 

#' Extract HLA Reads
#'
#' Extracts HLA reads based on bowtie2 alignment for each subject.
#' 
#' For internal wgsHLAfiltR use.
#' @param sampleList The sampleList object returned by the general.paired_sample_objects() function.
#' @param threads The number of threads specified for bowtie2. 
#' @param extractedFastqDirectory The directory into which filtered fastq files will be writen. 
#' @param forceRun A logical value. If forceRun = TRUE, the alignment is forced to run, even if there are previous results.
#' @return The populated sampleList object.
#' @note This function is for interal use. 
#' @export

extractor.run <- function(sampleList, threads, extractedFastqDirectory, forceRun=F){
  
  # Create the directory if it does not yet exist
  if(!file.exists(extractedFastqDirectory)){
    dir.create(extractedFastqDirectory,recursive=T)
    cat('\n',extractedFastqDirectory,'created.\n')
  }else{
    cat('\n',extractedFastqDirectory,'found.\n')
  }
  
  ## Check to make sure samtools is accessible
  bowtie2 <- system2('which', c('bowtie2'), stdout=T, stderr=T)
  check.system2_output(bowtie2, 'bowtie2 not found. Please install bowtie2.')
  
  for(currentSample in sampleList){
    cat('\n\nProcessing sample',currentSample$name,'----------')
    
    ## forceRun=T forces the alignment to run, even if there are previous results
    if(!forceRun){
      previousResultVect <- grep(currentSample$name, list.files(extractedFastqDirectory,full.names=T),value=T)
      
      ## If there are two files found and they both have _HLA_ in them, 
      #  add the file names to the sample object and skip alignment
      if(length(previousResultVect) == 2 & length(grep('_HLA_',previousResultVect,fixed=T)) == 2){
        cat('\nExtracted files found, skipping alignment.')
        currentSample[['HLAfastq1path']] <- grep('HLA_1.fastq.gz',previousResultVect,value=T)
        currentSample[['HLAfastq2path']] <- grep('HLA_2.fastq.gz',previousResultVect,value=T)
        
        next
      }
    }
    
    currentSample <- extractor.bowtie2_align(bowtie2, threads, currentSample, extractedFastqDirectory)
  }
  
  cat("\n\n----- HLA read extraction is complete. Extracted reads are in",extractedFastqDirectory,'-----\n')
  
  return(sampleList)
}

## SeqSample Class

#' Sequenced Samples Object Class
#'
#' An object class function required for the package.
#' @field name A charactrer string describing the name of the fastq file
#' @field rawfastq1path A charactrer string describing the path to the first unfiltered fastq
#' @field rawfastq2path A charactrer string describing the path to the paired unfiltered fastq
#' @field kirfastq1path A charactrer string; not used in this package
#' @field kirfastq2path A charactrer string; not used in this package
#' @field geneContent A list; not used in this package
#' @field kffHits A list; not used in this package
#' @field copyNumber A list; not used in this package
#' @field failed A logical
#' @field haploType A list; not used in this package
#' @field filterType A list; not used in this package
#' @field gzip A logical identifying wether the unfiltered fastq files are gzipped or not
#' @field samPath A character string describing the path to the sam file for the extracted fastq
#' @field bamPath A character string; not used in this package
#' @import methods
#' @export seqSample
#' @exportClass seqSample
seqSample <- setRefClass("seqSample",
                         fields=list(name='character',
                                     rawfastq1path='character',
                                     rawfastq2path='character',
                                     kirfastq1path='character',
                                     kirfastq2path='character',
                                     geneContent='list',
                                     kffHits='list',
                                     copyNumber='list',
                                     failed='logical',
                                     haploType='list',
                                     filterType='list',
                                     gzip='logical',
                                     samPath='character',
                                     bamPath='character'))
