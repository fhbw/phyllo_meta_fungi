#!/usr/bin/Rscript

rm(list=ls())

start_time <- Sys.time()

getwd()

# Import libraries
library(tidyverse);packageVersion("tidyverse")
library(dada2);packageVersion("dada2")
library(ShortRead);packageVersion("ShortRead")

#######################
## Input information ##
#######################

path <- "./" # Directory of fastq files
for_pattern = "_R1" # Forward reads file pattern
FWD <- "ACCWGCGGARGGATCATTA" # Forward primer sequence
REV <- "TTYRCKRCGTTCTTCATCG" # Reverse primer sequence
unite.ref <- "../sh_general_release_dynamic_all_04.02.2020.fasta" # Path to taxonomy database
maxEE_for = 1 # Maximum expected errors for FORWARD reads
cpu = 40 # No. of threads
bootstrap = 80 # Minimum bootstrapping support required to return a taxonomic classification

##################
## Run analysis ##
##################
print("")
print("######################################################################################")
print("## Script partially adapted from https://benjjneb.github.io/dada2/ITS_workflow.html ##")
print("######################################################################################")  
print("")

# Select and check fastq files
list.files(path)
fnFs <- sort(list.files(path, pattern = for_pattern,
                        full.names = TRUE))
head(fnFs)

# Make directory for figures
dir.create("Figures")

# Plot quality of raw reads
pq1 <- plotQualityProfile(fnFs)
ggsave("./Figures/plot_quality_raw.pdf", pq1, width = 20, height = 15)

## Trimming and quality filtering
# Filter N containing reads
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = cpu)

# Define all read directions of library primers
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
FWD.orients
REV.orients

# Check for primer sequences in raw reads N-filtered files
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
} 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[3]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[3]]))

cutadapt <- "/home/falk.behrens/.local/bin/cutadapt"
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt") # Put trimmed files in cutadapt/ subdirectory
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))

# Settings for cutadapt
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
threads <- paste("-j", cpu)
untrimmed <- paste("--discard-untrimmed")

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(threads, R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], # output files
                             fnFs.filtN[i], untrimmed)) # input files
}

# Check for priemr sequences in raw reads N-filterd files
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[3]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[3]]))

# filenames of trimmed files
cutFs <- sort(list.files(path.cut, pattern = for_pattern, full.names = TRUE))


# Quality filtering of trimmed reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))# Put filterd files in filtered/ subdirectory

# Perform quality filtering 
out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = maxEE_for, truncQ = 2, rm.phix=TRUE,
                     minLen = 100, compress = TRUE, multithread = cpu, verbose = TRUE)
out

# filtering zero count samples
df <- data.frame(out)
df$path <- filtFs
df <- df[df$reads.out != 0,] 
filtFs <- df$path
out <- out[out[,2] != 0,] 

# Plot quality of raw reads  
pq2 <- plotQualityProfile(filtFs)
ggsave("./Figures/plot_quality_filtered.pdf", pq2,  width = 20, height = 15)


## Denoise reads with DADA2
# dereplication
derepFs <- derepFastq(filtFs, verbose = TRUE)

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(filtFs, get.sample.name))
head(sample.names)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

# Learn errors
errF <- learnErrors(filtFs, nbases = 1e+09,
                    multithread = cpu, randomize = TRUE,
                    verbose = TRUE)

# Plot errors
perr <- plotErrors(errF)
ggsave("./Figures/plot_errors.pdf", perr,  width = 20, height = 15)

# Denoise reads
dadaFs <- dada(derepFs, err = errF, multithread = cpu, pool = FALSE)

# Make sequence table
seqtabFs <- makeSequenceTable(dadaFs)
dim(seqtabFs)

# Remove chimeric sequences
seqtab.nochimFs <- removeBimeraDenovo(seqtabFs, method="consensus", multithread = cpu, verbose = TRUE)
dim(seqtab.nochimFs)

# Combine duplicated sequences
seqtab.nochimFs <- collapseNoMismatch(seqtab.nochimFs)
dim(seqtab.nochimFs)

# Summarize sequence processing
hist_tab <- data.frame(table(nchar(getSequences(seqtab.nochimFs))))
length_dist <- ggplot(hist_tab, aes(x=Var1, y=Freq))+
  geom_histogram(stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
ggsave("./Figures/hist_length_dist.pdf", length_dist,  width = 20, height = 5)

getN <- function(x) sum(getUniques(x))
trackFs <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochimFs))
colnames(trackFs) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(trackFs) <- sample.names
trackFs

## Assign taxonomy and make final ASV table
# Assign taxonomy
set.seed(100)
taxa <- assignTaxonomy(seqtab.nochimFs, unite.ref, minBoot = bootstrap, multithread = cpu, tryRC = TRUE, verbose = TRUE)

# Prepare otutab
ASV_otumat <- t(seqtab.nochimFs) 

# Adjust format of taxa information
ASV_taxmat <- data.frame(taxa)
ASV_taxmat$Kingdom <- gsub("k__","", ASV_taxmat$Kingdom)
ASV_taxmat$Phylum <- gsub("p__","", ASV_taxmat$Phylum)
ASV_taxmat$Class <- gsub("c__","", ASV_taxmat$Class)
ASV_taxmat$Order <- gsub("o__","", ASV_taxmat$Order)
ASV_taxmat$Family <- gsub("f__","", ASV_taxmat$Family)
ASV_taxmat$Genus <- gsub("g__","", ASV_taxmat$Genus)
ASV_taxmat$Species <- gsub("s__","", ASV_taxmat$Species)
ASV_taxmat[is.na(ASV_taxmat)] <- "Unknown"
ASV_taxmat$Kingdom = gsub("\\unidentified", "Unknown", ASV_taxmat$Kingdom)
ASV_taxmat$Species = gsub("\\k__unidentified", "Unknown", ASV_taxmat$Species)
ASV_taxmat$Species <- paste(ASV_taxmat$Genus, ASV_taxmat$Species)
ASV_taxmat$Species = gsub("\\ ", "_", ASV_taxmat$Species)
ASV_taxmat$Species = gsub("\\Unknown_Unknown", "Unknown", ASV_taxmat$Species)
ASV_taxmat$Species = gsub("\\_Unknown", "_sp.", ASV_taxmat$Species)

# Make final merged.table
ASV.merged.table <- merge(ASV_otumat,ASV_taxmat, by ="row.names")
ASV.merged.table <- ASV.merged.table[c(2:(length(colnames(ASV.merged.table))),1)]
colnames(ASV.merged.table)[length(colnames(ASV.merged.table))] <- "Sequence"
ASV.merged.table$abundance<-rowSums(ASV.merged.table[,c(1:(as.numeric(length(colnames(ASV.merged.table)))-9))])
ASV.merged.table <- ASV.merged.table[order(ASV.merged.table$abundance,decreasing = TRUE),]
ASV.merged.table$abundance <- NULL
rownames(ASV.merged.table) <- paste0("ASV_", 1:nrow(ASV.merged.table))
ASV.merged.table$ID <- rownames(ASV.merged.table)
ASV.merged.table <- ASV.merged.table[c(1:(length(colnames(ASV.merged.table))-2),length(colnames(ASV.merged.table)),(length(colnames(ASV.merged.table))-1))]
write.table(ASV.merged.table,"ASV.merged.table.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE,fileEncoding="UTF-8")

finish_time <- Sys.time()

print(paste("Start ", start_time))
print(paste("Finish ", finish_time))
