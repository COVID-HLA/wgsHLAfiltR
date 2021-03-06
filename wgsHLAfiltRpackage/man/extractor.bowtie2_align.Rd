% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wgsHLAfiltR_functions.R
\name{extractor.bowtie2_align}
\alias{extractor.bowtie2_align}
\title{Bowtie2 Alignment of Reads}
\usage{
extractor.bowtie2_align(
  bowtie2_command,
  threads,
  currentSample,
  extractedFastqDirectory
)
}
\arguments{
\item{bowtie2_command}{The command being sent to bowtie2 in the external system.}

\item{threads}{The number of threads specified for bowtie2.}

\item{currentSample}{A specific element of the sampleList object returned by the general.paired_sample_objects() function.}

\item{extractedFastqDirectory}{The directory where filtered fastq.gz files are written.}
}
\value{
The aligned reads for an individual subject.
}
\description{
Performs bowtie2 alignments of reads for a subject.
}
\details{
For internal wgsHLAfiltR use.
}
\note{
This function is for interal use.
}
