##  ###################################################  ##
##  Processing raw ITS seqs into a phyloseq object       ##
##  This script processes HN00158973 w/o pitcher plant   ##
##                                                       ##
##  Nudibranch fungi microbiome                          ##
##                                                       ##
##  Author: Ming Sheng - June 15, 2022                   ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.2.0                                            ##
##  dada2 v 1.24.0                                       ##
##  purrr v 0.3.4                                        ##
##  tidyverse v 1.3.1                                    ##
##  readxl v 1.4                                         ##
##  decontam v 1.16                                      ##
##  phyloseq v 1.40.0                                    ##
##                                                       ##
##  ###################################################  ##


# Load packages ####
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(tidyverse); packageVersion("tidyverse")
library(readxl); packageVersion("readxl")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")

# Load metadata ####
fungi_meta <- readxl::read_xlsx("MetaData_fung.xlsx")

# remove extra cols
full_meta <- fungi_meta %>% select(c("Sample ID","Species", "Location"))
names(full_meta)[1:1] <- 'SampleID'

# Add column identifying PCR negatives
full_meta$PCR_Negative <- FALSE
full_meta$PCR_Negative[grep("Blank",full_meta$SampleID)] <- TRUE

# Find raw fastq files and prepare workspace ####
path <- "./HN00158973" 
fqs <- list.files(path, full.names = TRUE, recursive = FALSE, pattern = "_1.fastq.gz$|_2.fastq.gz$")

# Parse fwd and rev reads
fnFs <- fqs[grep("_1.fastq.gz",fqs)]
fnRs <- fqs[grep("_2.fastq.gz",fqs)]

# Get Sample Names
sample.names <- str_remove(basename(fnFs),"_1.fastq.gz")

# subset metadata to this run's samples
meta <- full_meta %>% filter(SampleID %in% sample.names)
rm(fungi_meta); rm(full_meta) # quick cleanup of environment

# Peek at quality profiles
plotQualityProfile(fnFs[c(1,30)]) # fwd reads # drops below Q30 at around 270
plotQualityProfile(fnRs[c(1,30)]) # rev reads # drops below Q30 at around 160

# Make filtered outfile names
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# make new directory for filtered files
if(!dir.exists(file.path(path,"filtered"))){
  dir.create(file.path(path,"filtered"))
}

# check for duplicated sample names
sum(duplicated(sample.names))

# CHECK FOR AND REMOVE PRIMER SITES WITH CUTADAPT ####

############# You will need to change these two values to match your data #############

# Here, you should supply the primer sequences used during PCR
FWD <- "TTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC" # Sequence of FWD primer
REV <- "CGTTCTTCATCGATGCVAGARCCAAGAGATC"  # Sequence of REV primer

######################################################################################################

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients; REV.orients

# Prefilter to remove reads with ambiguous (N) bases ####
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE) # on Windows, set multithread = FALSE

# Discover primer matches, regardless of orientation ####
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Run cutadapt ####
# If the following command returns an error, you do not have cutadapt installed correctly
system2("py -m cutadapt", args = "--version") # running in windows

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2("py -m cutadapt", args = c(R1.flags, R2.flags, "-n", 2, "--minimum-length 100", # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# sanity check
# This should show no occurences in any of the orientations now
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Filter and trim ####

# cut fwd reads at 300 and rev reads at 200
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,200),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,verbose = TRUE)

#Check trimmed sequence quality
plotQualityProfile(filtFs[c(20,30)])
plotQualityProfile(filtRs[c(20,30)])

saveRDS(out,"./Output/out.RDS") # save object
out <- readRDS("./Output/out.RDS") # reload point
##HG06 has 9 reads.in but 0 reads.out

filtpath <- file.path(path,"filtered")

# reassign filts for any potentially lost samples
filtFs <- list.files("./HN00158973/filtered", pattern = "_F_filt", full.names = TRUE)
filtRs <- list.files("./HN00158973/filtered", pattern = "_R_filt", full.names = TRUE)

# Get Sample Names (again, just in case)
sample.names <- str_remove(basename(fnFs),"_1.fastq.gz")

# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF,"./Output/errF.RDS") # save object
saveRDS(errR,"./Output/errR.RDS") # save object
errF <- readRDS("./Output/errF.RDS") # reload point
errR <- readRDS("./Output/errR.RDS") # reload point

# plot error rates for sanity
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Infer sequence variants ####

# add names to filts
sample.names <- sample.names[-13] #HG06 has no reads through
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Dereplication, sample inferrence, and merging ####

# loop through each pair of fwd and rev reads, one file at a time
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

saveRDS(mergers,"./Output/mergers.RDS")

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table ####
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab,"./Output/seqtab.RDS") # save object
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim,"./Output/seqtab.nochim.RDS")
seqtab.nochim = readRDS("./Output/seqtab.nochim.RDS")

# Assign taxonomy -  / 80% bootstrap min ####
taxa <- assignTaxonomy(seqtab.nochim,"./fungal_database/sh_general_release_dynamic_s_10.05.2021.fasta", minBoot = 80,multithread = TRUE)
saveRDS(taxa,"./Output/taxa.RDS")

taxa <- addSpecies(taxa, "./fungal_database/sh_general_release_dynamic_s_10.05.2021.fasta")
saveRDS(taxa,"./Output/taxa_with_spp.RDS")

# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.nochim)
row.names(seqtab.df)

# Create phyloseq object ####

# subset to remove missing samples
in.meta <- which(names(seqtab.nochim[,1]) %in% meta$SampleID == TRUE)
seqtab.nochim <- seqtab.nochim[in.meta,]
dim(seqtab.nochim)
in.seqtab <- which(meta$SampleID %in% names(seqtab.nochim[,1]))
meta <- meta[in.seqtab,]

# re-order
meta <- meta[order(meta$SampleID),]
row.names(meta) <- meta$SampleID
seqtab.nochim <- (seqtab.nochim[row.names(meta),])
identical(row.names(seqtab.nochim), as.character(meta$SampleID))

# make phyloseq object
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
met <- sample_data(meta)
tax <- tax_table(taxa)

sample_names(met) <- met$SampleID
ps <- phyloseq(otu,met,tax)

# save it
saveRDS(ps, "./Output/raw_ps_object.RDS")

# Identify and remove contaminants
blanks = which(ps@sam_data$PCR_Negative == TRUE)
contamdf.prev <- isContaminant(ps, neg=blanks, threshold = 0.01)
table(contamdf.prev$contaminant)

ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)

# save contam-free phyloseq object
saveRDS(ps.noncontam, "./Output/noncontam_ps_object.RDS")
