## 16S Metagenomics Analysis DADA 2 #######
###########################################

library(dada2); packageVersion("dada2")

path<-("~/Documents/Metagenomics/metagenome seq data-MMME/RAW-FASTQ/")

# file Import 

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles
# visualizing the quality p#rofiles of the forward reads:

plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnRs[1:10])

#Filter and trim

#Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,130),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=T)


#Learn the Error Rates
#The DADA2 algorithm makes use of a parametric error model (err) and 
#every amplicon dataset has a different set of error rates. 
#The learnErrors method learns this error model from the data, by alternating estimation of the error rates 
#and inference of sample composition until they converge on a jointly consistent solution. As in many 
#machine-learning problems, the algorithm must begin with an initial guess, for which the maximum possible 
#error rates in this data are used (the error rates if only the most abundant sequence is correct and all 
#the rest are errors).

errF <- learnErrors(filtFs, multithread=TRUE)

#Dereplication

#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding 
#“abundance” equal to the number of reads with that unique sequence. Dereplication substantially reduces 
#computation time by eliminating redundant comparisons.

#derepFs <- derepFastq(filtFs, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
#names(derepFs) <- sample.names
#names(derepRs) <- sample.names
