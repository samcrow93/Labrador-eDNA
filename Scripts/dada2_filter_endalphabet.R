# Load dada2 and check version
library(dada2)
packageVersion("dada2")
library(here)

# set path to tutorial files
pathF <- here("Small dataset","Forward")
pathR <- here("Small dataset", "Reverse")
filtpathF <- file.path(pathF, "filtered2") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered2") # ...
fastqFs <- sort(list.files(pathF, pattern=c("(^12S-ALX21-L|^12S-11).*\\_R1.fastq.gz")))
fastqRs <- sort(list.files(pathR, pattern=c("(^12S-ALX21-L|^12S-11).*\\_R2.fastq.gz")))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(100,100), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)