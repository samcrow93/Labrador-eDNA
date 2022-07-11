# Load dada2 and check version
library(dada2)
packageVersion("dada2")
library(here)

# set path to tutorial files
pathF <- here("Small dataset","Forward")
pathR <- here("Small dataset", "Reverse")
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern=".fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern=".fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(100,100), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

# Infer Sequence Variants
filtpathF <- here("Small dataset","Forward","filtered") # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- here("Small dataset","Reverse","filtered") # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern=".fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern=".fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, here("Small dataset", "seqtab.Rds")) # CHANGE ME to where you want sequence table saved

# Merge multiple runs (if necessary)
st1 <- readRDS("path/to/run1/output/seqtab.rds")
st2 <- readRDS("path/to/run2/output/seqtab.rds")
st3 <- readRDS("path/to/run3/output/seqtab.rds")
st.all <- mergeSequenceTables(st1, st2, st3)
# Remove chimeras
seqtab2 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
