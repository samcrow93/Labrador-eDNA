# Load dada2 and check version
library(dada2)

# set path to tutorial files for 12S teleo marker
pathF <- "/home/samcrow/scratch/eDNA/LAB_EDNA/"
pathR <- "/home/samcrow/scratch/eDNA/LAB_EDNA/"
filtpathF <- file.path(pathF, "filtered_teleo12S_test") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered_teleo12S_test") # ...
fastqFs <- sort(list.files(pathF, pattern="^12S.*\\_R1.fastq.gz")) #select just 12Steleo samples
fastqRs <- sort(list.files(pathR, pattern="^12S.*\\_R2.fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: used parameters specified in CEGA publication
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              maxEE=2, minQ=2, truncQ=2, maxN=0, compress=TRUE, multithread=TRUE)