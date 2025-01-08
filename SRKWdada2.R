### QAQC pipeline for Wild Orca SRKW sequences
### Sofia Kaiaua and Amy Van Cise
### Fall 2024

#### Set up environment --------------------------------------------------------
library(dada2)
library(tidyverse)
library(kableExtra)
library(phyloseq)
library(Biostrings)
library(ggplot2)

theme_set(theme_bw())

#### Get raw sequence data -----------------------------------------------------
path <- "M:\\SRKW Diet Metabarcoding\\Sequence data\\16SP1"
list.files(path)

#Separate Forward and Reverse Fastq Files
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect Read Quality Profiles
#Plot Forward Reads
#plotQualityProfile(fnFs[25])
#Plot Reverse Reads
#plotQualityProfile(fnRs[1:2])

### Filter and trim raw reads --------------------------------------------------
#Create subdirectory for filtered reads
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 40, truncLen=c(280,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)

head(out)

### Generate sequences ---------------------------------------------------------
#Learn Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#Merge Paired Reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])
tail(mergers[[1]])

#Construct Sequence Table
seqtab <- makeSequenceTable(mergers)

#sanity check
dim(seqtab)
#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
#sanity check
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

### Track Reads Through Pipeline -----------------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

#set row and column names
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

### Assign Taxonomy ------------------------------------------------------------
taxa <- assignTaxonomy(seqtab.nochim, 
                       "G:/My Drive/00 UW/00.5 W.A.D.E. lab resources/Intern Projects/SRKW diet metabarcoding/04 Data analysis/16S_salmon_groundfish_reference_database_2022.fasta")


### Combine into phyloseq object -----------------------------------------------
# Get sample metadata
samdf <- read.csv("G:/My Drive/00 UW/00.5 W.A.D.E. lab resources/Intern Projects/SRKW diet metabarcoding/04 Data analysis/SRKW_WO_SDZWA_EditedforR.csv") %>%
  column_to_rownames(var = "LabID")

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

save(ps,taxa,samdf,seqtab.nochim,track,file="G:/My Drive/00 UW/00.5 W.A.D.E. lab resources/Intern Projects/SRKW diet metabarcoding/04 Data analysis/SRKW_WO_16Sdiet_dada2out.Rdata")

