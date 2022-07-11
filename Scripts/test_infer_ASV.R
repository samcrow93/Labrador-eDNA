# Infer Sequence Variants
filtpathF <- "/home/samcrow/scratch/eDNA/LAB_EDNA/filtered_teleo12S_test/" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/home/samcrow/scratch/eDNA/LAB_EDNA/filtered_teleo12S_test/" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="12S-W.*\\_R1.fastq.gz", full.names = TRUE) #for test run, choose only files/sites that start with "W"
filtRs <- list.files(filtpathR, pattern="12S-W.*\\_R2.fastq.gz", full.names = TRUE)
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
# Construct sequence table
seqtab <- makeSequenceTable(mergers)
# Remove chimeras
seqtab2 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
# save sequence table
saveRDS(seqtab2, "/home/samcrow/scratch/eDNA/LAB_EDNA/filtered_teleo12S_test/seqtabW_12S.rds")