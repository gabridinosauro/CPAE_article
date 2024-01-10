# First, create barcode mapping file in text editor with 2 columns that are the barcodes and the sample ids
# Read file into R with the "barcode" variable
###### Make sure to change folder path on line 9
###### Make sure to check and change column names throughout
###### Make sure to change what file name you want in the end on line 19

library(seqinr) # needed for function below

setwd('/Users/gibugs/Box/Peds-CRC/Illumina/CPAE_Lundborg3151472/CPAE_all_samples_R/')
setwd('/Users/daniellaubitz/Box/Peds-CRC/Illumina/CPAE_Lundborg3151472/CPAE_all_samples_R/')

barcode<-read.table('CPAE_lundborg3151472_barcodes.txt', sep='\t', header=T)
barcode[,3]<-as.character(barcode[,3])

barcode_revcomp_vec<-sapply(barcode[,3], function(x){toupper(c2s(rev(comp(s2c(x)))))})

barcode_revcomp<-data.frame(barcode_revcomp_vec, barcode[,1])
colnames(barcode_revcomp)<-c("Barcode", "SampleID")

write.table(barcode_revcomp, "barcodes_revcomp.txt", sep='\t', quote=F, row.names=F)

# Packages ----------------------------------------------------------------



library(dada2)
library(tidyverse)
library(vegan)
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(plyr)
library(RColorBrewer)

theme_set(theme_bw())


# DADA2 -------------------------------------------------------------------


# Define working folder path 
path <- "/Users/gibugs/Box/Peds-CRC/MiSeqAnalysis/200213_M03190_0099_000000000-CW88V/Alignment_1/20200214_133917/Fastq/demux/"
list.files(path, pattern = "*")

# Original tutorial comment: Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
# Claire: do list.files(pattern="*") to see files. Our file formats are different than above.
fnFs <- sort(list.files(path, pattern="_R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_", full.names = TRUE))

# Claire: Remove unassigned sequences (from samples were not interested in) and those belongs to another project
fnFs = grep("CPAE[0-9]", fnFs, value = T, perl = T)
fnRs = grep("CPAE[0-9]", fnRs, value = T, perl = T)

# Check for matching lengths
length(fnFs) == length(fnRs) # TRUE

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# Claire: Changed 1 to 2 to get part of file paths after "_"
# Daniel: changed to 6, it's 6th element of the file name and the parts are separated with "_"
sample.names0 <- sapply(strsplit(basename(fnFs), "_"), `[`, 6)

# Claire: Once again, need to deviate from tutorial to get names we want
sample.names = gsub(".fastq.gz", "", sample.names0)

# Inspect quality of sequences in files
plotQualityProfile(fnFs[11:30]) # forward
ggsave("qualityF.pdf")
plotQualityProfile(fnRs[11:30]) # reverse
ggsave("qualityR.pdf")

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Trim forwards at 140bp and reverses at 145bp based on quality plots; also trim first 10 bp as recommended
# After this completes, check "out" to see if filtering was too stringent
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(145,145),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft = 10) # On Windows set multithread=FALSE


tail(out)

#read mapping file
md = read.delim("CPAE_md.txt", header = T, sep = "\t")
row.names(md) <- md$SampleID


# Claire: How much of the data was retained in each sample? The % column is added to out
out = as.data.frame(out)
out$Retained = out$reads.out/out$reads.in*100
out$samplename = sample.names

# Model and learn error rates in filtered sequences
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


# Plot error rates that are based on nucleotide transition probabilities and quality scores; observed scores (black) should match red lines (expected)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Reduce computation time by dereplicating sequences
# /Users/gibugs/Box/Peds-CRC/MiSeqAnalysis/200213_M03190_0099_000000000-CW88V/
# Alignment_1/20200214_133917/Fastq/demux//filtered/
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


# Join forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Claire: Inspect the merger data.frame from the first sample
head(mergers[[2]])
mergers[[2]]$nmatch

# Claire: Make amplicon sequence variants table (like OTU table)
seqtab <- makeSequenceTable(mergers)
dim(seqtab) #[1]   94 5414

# Inspect sequence lengths: Most amplicon sequence variants are 232 to 234bp
table(nchar(getSequences(seqtab)))

# remove non-target-length sequences with base R manipulations of the sequence table:
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(232,234)]
dim(seqtab2) #[1]   94 5316. This table contains 5316. ASVs
# Inspect sequence lengths: Most amplicon sequence variants are 232 to 234bp
table(nchar(getSequences(seqtab2)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) # 94 2921

# How much of the dataset is left after chimera identification and removal?
sum(seqtab.nochim)/sum(seqtab)*100 # ~97.68% left

# Dada2 analysis pipeline stats
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "Retained %", "SampleID","denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names


# Claire: Total percent of dataset remaining
track$TtRetained = track$nonchim/track$input*100
head(track)

row.names(md) == row.names(track)

track$status <- md$status
track$sample_type <- md$sample_type
track1 <- track[order(track$sample_type,track$status),]
track1$type_status <- paste0(track1$sample_type,"_" ,track1$status)

# does numebr of reads relates to sample type and/or patient status?
ggplot(track1,aes(x=type_status,y=nonchim))+
  geom_point(aes(colour = status), cex=2)+
  xlab("sample type / status") +
  ylab("number of reads (after QC)")+
  theme(text = element_text(size=8), axis.text.x=element_text(angle=45, hjust=1))
ggsave("number of reads per sample type.pdf", width = 6, height = 4)

# Went to website (https://benjjneb.github.io/dada2/training.html) and downloaded silva taxonomy data with following commands in terminal: 
# wget -N "https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1"  
# wget -N "https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1"
# both file have to be in working directory

# Assign taxonomy using the RDP classifier to ASVs.
#Files are in main 'Illumina" folder
taxa <- assignTaxonomy(seqtab.nochim, "/Users/gibugs/Box/Peds-CRC/Illumina/silva_nr_v132_train_set.fa.gz?download=1", multithread=TRUE)
taxa.no.species <- assignTaxonomy(seqtab.nochim, "/Users/gibugs/Box/Peds-CRC/Illumina/silva_nr_v132_train_set.fa.gz?download=1", multithread=TRUE)


# Assign species level taxonomy to ASVs with exact 100% matching
taxa <- addSpecies(taxa, "/Users/gibugs/Box/Peds-CRC/Illumina/silva_species_assignment_v132.fa.gz?download=1")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Save taxonomy table and asv table
write.table(taxa, "taxa_table.txt", sep = "\t", quote = F)
write.table(seqtab.nochim, "asv_table.txt", sep = "\t", quote = F)
write.table(track, "sequence_pipeline_stats.txt", sep = "\t", quote = F)


# POST DADA2 Gabri --------------------------
setwd("/Users/gabri/Library/CloudStorage/Box-Box/Illumina/CPAE_Lundborg3151472/CPAE_all_samples_R")
seqtab.nochim = read.table("asv_table.txt")
taxa = read.table("taxa_table.txt", row.names = 1)
md = read.delim("CPAE_md.txt", header = T, row.names = 1)
unique(md$status)
md$status = ifelse(md$status == "CNTL", "CTRL", md$status)
md$status = ifelse(md$status == "PANDAS", "CPAE", md$status)
unique(md$status)

# Same number of samples in metadata and asv table?
#because I do not want to change everything I just saved seqtab.nochim as seqtab.nochim2
seqtab.nochim2 = seqtab.nochim
nrow(md) == nrow(seqtab.nochim2) #TRUE

md[,1] # to check, it hsould be 94 levels

# Do rownames match?
rownames(md) == rownames(seqtab.nochim2) # all TRUE
hist(rowSums(seqtab.nochim2))
plot(rowSums(seqtab.nochim2))
rownames(seqtab.nochim2)
sort(rowSums(seqtab.nochim2))
stnc = as.data.frame(seqtab.nochim2)

# Visualize and save sequence counts that have made it through quality control
sort <- sort(rowSums(stnc))
sort.md <- merge(sort, md, by = "row.names", all=T)
sort.md <- sort.md[order(sort.md$x),]
sort.md
#write.table(sort, "seq count per sample sorted.txt", sep = "\t", quote = F)
ggplot(stnc, aes(rowSums(stnc))) + geom_histogram() 



## separate samples based on the sample_type
nasal.md <- droplevels(md[md$sample_type %in% "NASAL",])
nasal.stnc <- droplevels(stnc[row.names(stnc) %in% row.names(nasal.md),])

throat.md <- droplevels(md[md$sample_type %in% "THROAT",])
throat.stnc <- droplevels(stnc[row.names(stnc) %in% row.names(throat.md),])

stool.md <- droplevels(md[md$sample_type %in% "STOOL",])
stool.stnc <- droplevels(stnc[row.names(stnc) %in% row.names(stool.md),])


# Rarefy for nasal samples
sort(rowSums(nasal.stnc))
nasal.rar.otu = rrarefy(nasal.stnc, 1254)
nasal.rar.otu.df = as.data.frame(nasal.rar.otu)
rowSums(nasal.rar.otu.df)
# Remove unrarefied samples (should lose 0)
nasal.asv = nasal.rar.otu.df[rowSums(nasal.rar.otu.df) == 1254,]
dim(nasal.rar.otu.df)[1] - dim(nasal.asv)[1] # 4 samples were removed, both are blanks 
nasal.md = nasal.md[rownames(nasal.md) %in% rownames(nasal.asv),] # match the numbver of samples in md with the number of samples in ASVs table

# Rarefy for throat samples
sort(rowSums(throat.stnc))
throat.rar.otu = rrarefy(throat.stnc, 11707)
throat.rar.otu.df = as.data.frame(throat.rar.otu)
rowSums(throat.rar.otu.df)
# Remove unrarefied samples (should lose 0)
throat.asv = throat.rar.otu.df[rowSums(throat.rar.otu.df) == 11707,]
dim(throat.rar.otu.df)[1] - dim(throat.asv)[1] # 0 samples were removed, both are blanks 
throat.md = throat.md[rownames(throat.md) %in% rownames(throat.asv),] # match the numbver of samples in md with the number of samples in ASVs table

# Rarefy for stool samples
sort(rowSums(stool.stnc))
stool.rar.otu = rrarefy(stool.stnc, 23769)
stool.rar.otu.df = as.data.frame(stool.rar.otu)
rowSums(stool.rar.otu.df)
# Remove unrarefied samples (should lose 0)
stool.asv = stool.rar.otu.df[rowSums(stool.rar.otu.df) == 23769,]
dim(stool.rar.otu.df)[1] - dim(stool.asv)[1] # 0 samples were removed, both are blanks 
stool.md = stool.md[rownames(stool.md) %in% rownames(stool.asv),] # match the numbver of samples in md with the number of samples in ASVs table





############### ALPHA DIVERISTY on rarefied otu (data frame) #############


nasal.md$Richness <- specnumber(nasal.asv)
nasal.md$Shannon <- diversity(nasal.asv, index="shannon")
nasal.md$Simpson <- diversity(nasal.asv, index="simpson")
colnames(nasal.md)

stool.md$Richness <- specnumber(stool.asv)
stool.md$Shannon <- diversity(stool.asv, index="shannon")
stool.md$Simpson <- diversity(stool.asv, index="simpson")
colnames(stool.md)

throat.md$Richness <- specnumber(throat.asv)
throat.md$Shannon <- diversity(throat.asv, index="shannon")
throat.md$Simpson <- diversity(throat.asv, index="simpson")
colnames(throat.md)



#### Plot Alpha Diversity 
library(ggprism)
library(patchwork)
library(magrittr)
library(rstatix)
library(ggpubr)
library(rstatix)

### Nasal MD -----------
result <- t.test(Richness ~ status, data = nasal.md)$p.value
result <- signif(result, digits = 3)

df_p_val <- data.frame(
  group1 = "CPAE",
  group2 = "CTRL",
  label = result,
  y.position = max(nasal.md$Richness) + 5,
  status = "NA"
)


a = ggplot(nasal.md, aes(x=status, y=Richness, fill=status)) +
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = status))+ 
  xlab(NULL) + 
  ylab("Number of ASVs")+
  scale_fill_manual(values=c('darkred',"darkblue" , "white"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Nasal")+
  add_pvalue(df_p_val,
             xmin = "group1",
             xmax = "group2",
            label = "label",
            y.position = "y.position") 
a
## throat
result <- t.test(Richness ~ status, data = throat.md)$p.value
result <- signif(result, digits = 3)

df_p_val <- data.frame(
  group1 = "CPAE",
  group2 = "CTRL",
  label = result,
  y.position = max(throat.md$Richness) + 5,
  status = "NA"
)


b = ggplot(throat.md, aes(x=status, y=Richness, fill=status)) +
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = status))+ 
  xlab(NULL) + 
  ylab(NULL)+
  scale_fill_manual(values=c('darkred',"darkblue" , "white"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Throat")+
  add_pvalue(df_p_val,
             xmin = "group1",
             xmax = "group2",
             label = "label",
             y.position = "y.position") 

b

## stool
result <- t.test(Richness ~ status, data = stool.md)$p.value
result <- signif(result, digits = 3)

df_p_val <- data.frame(
  group1 = "CPAE",
  group2 = "CTRL",
  label = result,
  y.position = max(stool.md$Richness) + 5,
  status = "NA"
)


c = ggplot(stool.md, aes(x=status, y=Richness, fill=status)) +
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = status))+ 
  xlab(NULL) + 
  ylab(NULL)+
  scale_fill_manual(values=c('darkred',"darkblue" , "white"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Stool")+
  add_pvalue(df_p_val,
             xmin = "group1",
             xmax = "group2",
             label = "label",
             y.position = "y.position") 


### Shannon --------
## Nasal MD -----------
result <- t.test(Shannon ~ status, data = nasal.md)$p.value
result <- signif(result, digits = 3)

df_p_val <- data.frame(
  group1 = "CPAE",
  group2 = "CTRL",
  label = result,
  y.position = max(nasal.md$Shannon) + 0.2,
  status = "NA"
)


d = ggplot(nasal.md, aes(x=status, y=Shannon, fill=status)) +
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = status))+ 
  xlab(NULL) + 
  ylab("Shannon H'")+
  scale_fill_manual(values=c('darkred',"darkblue" , "white"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  #ggtitle("Nasal diversity")+
  add_pvalue(df_p_val,
             xmin = "group1",
             xmax = "group2",
             label = "label",
             y.position = "y.position") 

## throat
result <- t.test(Shannon ~ status, data = throat.md)$p.value
result <- signif(result, digits = 0.2)

df_p_val <- data.frame(
  group1 = "CPAE",
  group2 = "CTRL",
  label = result,
  y.position = max(throat.md$Shannon) + 0.2,
  status = "NA"
)


e = ggplot(throat.md, aes(x=status, y=Shannon, fill=status)) +
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = status))+ 
  xlab(NULL) + 
  ylab(NULL)+
  scale_fill_manual(values=c('darkred',"darkblue" , "white"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  #ggtitle("throat diversity")+
  add_pvalue(df_p_val,
             xmin = "group1",
             xmax = "group2",
             label = "label",
             y.position = "y.position") 



## stool
result <- t.test(Shannon ~ status, data = stool.md)$p.value
result <- signif(result, digits = 3)

df_p_val <- data.frame(
  group1 = "CPAE",
  group2 = "CTRL",
  label = result,
  y.position = max(stool.md$Shannon) + 0.2,
  status = "NA"
)


f = ggplot(stool.md, aes(x=status, y=Shannon, fill=status)) +
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = status))+ 
  xlab(NULL) + 
  ylab(NULL)+
  scale_fill_manual(values=c('darkred',"darkblue" , "white"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
#ggtitle("stool diversity")+
  add_pvalue(df_p_val,
             xmin = "group1",
             xmax = "group2",
             label = "label",
             y.position = "y.position") 

## 
ggarrange(a,b,c,d,e,f, heights = c(1.05, 0.95), labels = "AUTO")








#################### BETA DIVERSITY ######################
#NMDS ordination for nasal
nasal.micro.dist<-vegdist(nasal.asv, method="bray")
nasal.micro.dist.nmds<-metaMDS(nasal.micro.dist, k=2)
nasal.micro.dist.nmds$stress #0.161
stressplot(nasal.micro.dist.nmds)
nasal.md$NMDS001<-nasal.micro.dist.nmds$points[,1]
nasal.md$NMDS002<-nasal.micro.dist.nmds$points[,2]
stress_nasal = paste("stress =", round(nasal.micro.dist.nmds$stress, digits = 2))



#NMDS ordination for throat
throat.micro.dist<-vegdist(throat.asv, method="bray")
throat.micro.dist.nmds<-metaMDS(throat.micro.dist, k=2)
throat.micro.dist.nmds$stress #0.1196257
stressplot(throat.micro.dist.nmds)
throat.md$NMDS001<-throat.micro.dist.nmds$points[,1]
throat.md$NMDS002<-throat.micro.dist.nmds$points[,2]
stress_throat = paste("stress =", round(throat.micro.dist.nmds$stress, digits = 2))
#NMDS ordination for stool
stool.micro.dist<-vegdist(stool.asv, method="bray")
stool.micro.dist.nmds<-metaMDS(stool.micro.dist, k=2)
stool.micro.dist.nmds$stress #0.1585297
stressplot(stool.micro.dist.nmds)
stool.md$NMDS001<-stool.micro.dist.nmds$points[,1]
stool.md$NMDS002<-stool.micro.dist.nmds$points[,2]
stress= paste("stress =", round(stool.micro.dist.nmds$stress, digits = 2))


library(ggforce)
#ordination plots
ord1 = ggplot(nasal.md, aes(NMDS001, NMDS002))+
  geom_point(aes(color=status), size=4)+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=status, color = status, label = NA), con.type  = "none")+
  scale_color_manual(values=c('darkred', "darkblue"))+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  xlab("NMDS1") + 
  ylab("NMDS2")+
  theme_bw() +ggtitle("Nasal")+
  annotate("text",x=max(nasal.md$NMDS001),y=max(nasal.md$NMDS002),hjust=1,label= stress_nasal,  size = 3.5)
ord1

ord2 = ggplot(throat.md, aes(NMDS001, NMDS002))+
  geom_point(aes(color=status), size=4)+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=status, color = status, label = NA), con.type  = "none")+
  scale_color_manual(values=c('darkred', "darkblue"))+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  xlab("NMDS1") + 
  ylab("NMDS2")+
  theme_bw() +ggtitle("throat")+
  annotate("text",x=max(throat.md$NMDS001),y=max(throat.md$NMDS002),hjust=1,label= stress_throat,  size = 3.5)
ord2

ord3 = ggplot(stool.md, aes(NMDS001, NMDS002))+
  geom_point(aes(color=status), size=4)+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=status, color = status, label = NA), con.type  = "none")+
  scale_color_manual(values=c('darkred', "darkblue"))+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  xlab("NMDS1") + 
  ylab("NMDS2")+
  theme_bw() +ggtitle("stool")+
  annotate("text",x=max(stool.md$NMDS001),y=max(stool.md$NMDS002),hjust=1,label= stress,  size = 3.5)
ord3

ggarrange(ord1,ord2,ord3, common.legend = TRUE, nrow = 1)

adonis2(nasal.asv ~ status, nasal.md, permutations = 999,distance = "bray", strata = NULL)
adonis2(throat.asv ~ status, throat.md, permutations = 999,distance = "bray", strata = NULL)
adonis2(stool.asv ~ status, stool.md, permutations = 999,distance = "bray", strata = NULL)





adonis2(throat.asv ~ status, throat.md, permutations = 999,distance = "bray", strata = NULL)

ggplot(stool.md, aes(Axis01, Axis02))+
  geom_point(aes(fill=status), size=6, alpha=1, pch=21, col="black")+
  #geom_polygon(data = nasal.hull, aes(colour=status, fill=status), alpha = 0.5)+
  stat_ellipse(aes(color=status), type = "norm", linetype = 2) +
  stat_ellipse(aes(color=status),type = "t")+
  theme(legend.position="bottom", text = element_text(size = 18))+
  ggtitle("Stool")
ggsave("stool BC NMDS.pdf", height = 5, width = 7)

adonis2(stool.asv ~ status, stool.md, permutations = 999,distance = "bray")


patient.hull <- ddply(stool.md, "patientID", find_hull)
ggplot(stool.md, aes(Axis01, Axis02, color=patientID))+
  geom_point()+
  geom_text(aes(label=patientID),hjust=0, vjust=0)+
  geom_polygon(data = patient.hull, aes(colour=patientID), alpha = 0)#+
  #stat_ellipse(aes(color=status), type = "norm", linetype = 2) +
  #stat_ellipse(aes(color=patientID),type = "t")#+
  theme(legend.position="bottom", text = element_text(size = 18))+
  ggtitle("Stool")
ggsave("stool BC NMDS.pdf", height = 5, width = 7)

stool.md$unique <- paste0(stool.md$status,stool.md$patientID)



#################### Taxonomic analysis ###################


generate.tax.summary.modified<-function(otu, tax.table) {
  tax3.summary<-as.data.frame(apply((apply(otu, 1, function(x) by(x, tax.table$Phylum, sum))), 2, function(x){x/sum(x)}))
  tax4.summary<-as.data.frame(apply((apply(otu, 1, function(x) by(x, tax.table$Class, sum))), 2, function(x){x/sum(x)}))
  tax5.summary<-as.data.frame(apply((apply(otu, 1, function(x) by(x, tax.table$Order, sum))), 2, function(x){x/sum(x)}))
  tax6.summary<-as.data.frame(apply((apply(otu, 1, function(x) by(x, tax.table$Family, sum))), 2, function(x){x/sum(x)}))
  tax7.summary<-as.data.frame(apply((apply(otu, 1, function(x) by(x, tax.table$Genus, sum))), 2, function(x){x/sum(x)}))
  tax8.summary<-as.data.frame(apply((apply(otu, 1, function(x) by(x, tax.table$Species, sum))), 2, function(x){x/sum(x)}))
  return(list(tax3=tax3.summary, 
              tax4=tax4.summary, 
              tax5=tax5.summary, 
              tax6=tax6.summary, 
              tax7=tax7.summary, 
              tax8=tax8.summary))
}


nasal.taxa.sum= generate.tax.summary.modified(nasal.asv, as.data.frame(taxa)) # generates list
stool.taxa.sum= generate.tax.summary.modified(stool.asv, as.data.frame(taxa))
throat.taxa.sum= generate.tax.summary.modified(throat.asv, as.data.frame(taxa))
taxa.sum = generate.tax.summary.modified(stnc, as.data.frame(taxa))








####################################### GENUS LEVEL ####################################### 


# Selecting family-level table
taxa.sum.genus = taxa.sum$tax7

# Check that all columns are numeric before transformation (otherwise all values will be annoyingly converted to non-numeric)
# Are any columns not numeric?
sapply(taxa.sum.genus, class)[sapply(taxa.sum.genus, class) != "numeric"] # Yes, column Taxa is not

# Removing Taxa
fil.taxa.sum.genus = taxa.sum.genus[,!names(taxa.sum.genus) %in% "Taxa"]


# Transpose and convert to dataframe
t.fil.taxa.sum.genus = as.data.frame(t(fil.taxa.sum.genus))

# Add metadata
identical(rownames(t.fil.taxa.sum.genus), rownames(md)) # TRUE
md$status_type <- paste0(md$status,"_",md$sample_type)
genus.md = as.data.frame(cbind("status_type" = md$status_type, t.fil.taxa.sum.genus))
names(genus.md)

# Group and plot

genus.md.avg = genus.md %>%
  group_by(status_type) %>%
  summarise_all(list(mean)) %>%
  gather(key = "Taxa", value = "Avg_Relabund","1174-901-12":"Xylophilus") 

#remove taxa with Rel Abu less than 0.1%
genus.md.avg = genus.md.avg[genus.md.avg$Avg_Relabund > 0.005,]

mypalette2  <- c("#40004b","#ffffbf","#762a83","#fdbf6f","#9970ab","#ccebc5","#c2a5cf","#e7d4e8","#b8e186","#bebada",
                 "#de77ae","#f46d43","#f7f7f7","#d9f0d3","#a6dba0","#4d9221","#2166ac","#5aae61","#1b7837","#d73027",
                 "#00441b","#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#b2182b","#de77ae",
                 "#bc80bd","#f5f5f5","#c7eae5","#80cdc1","#35978f","#003c30",
                 "#8e0152","#fb8072","#c51b7d","#f1b6da","#fde0ef","#fdb462","#b3de69","#7fbc41","#276419",
                 "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#ff7f00","#cab2d6","#6a3d9a","#ffff99",
                 "#8dd3c7","#ffffb3","#80b1d3","#fdb462","#01665e","#fccde5","#d9d9d9",
                 "#40004b","#762a83","#fdbf6f","#9970ab","#ccebc5","#c2a5cf","#e7d4e8","#b8e186","#bebada",
                 "#de77ae","#f7f7f7","#d9f0d3","#a6dba0","#4d9221","#5aae61","#1b7837"
                 ,"#00441b","#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#bc80bd","#f5f5f5","#c7eae5","#80cdc1","#35978f","#003c30",
                 "#8e0152","#fb8072","#c51b7d","#f1b6da","#fde0ef","#fdb462","#b3de69","#7fbc41","#276419",
                 "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#ff7f00","#cab2d6","#6a3d9a","#ffff99",
                 "#8dd3c7","#ffffb3","#80b1d3","#fdb462","#01665e","#fccde5","#d9d9d9")

pdf("Taxa genus Level All samples 0.001.pdf", width = 10, height = 6)
ggplot(genus.md.avg, aes(x = status_type, y = Avg_Relabund, fill = Taxa)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mypalette2, guide = guide_legend(nrow=20)) +
  theme(legend.position="right", text = element_text(size = 12),  axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=c("CNTL_NASAL","PANDAS_NASAL",
                            "","CNTL_STOOL","PANDAS_STOOL",
                            "","CNTL_THROAT","PANDAS_THROAT"))
dev.off()




####################################### GENUS LEVEL FOR STOOL ####################################### 


# Selecting family-level table
stool.taxa.sum.genus = stool.taxa.sum$tax7

# Check that all columns are numeric before transformation (otherwise all values will be annoyingly converted to non-numeric)
# Are any columns not numeric?
sapply(stool.taxa.sum.genus, class)[sapply(stool.taxa.sum.genus, class) != "numeric"] # Yes, column Taxa is not

# Removing Taxa
fil.stool.taxa.sum.genus = stool.taxa.sum.genus[,!names(stool.taxa.sum.genus) %in% "Taxa"]


# Transpose and convert to dataframe
t.fil.stool.taxa.sum.genus = as.data.frame(t(stool.taxa.sum.genus))

# Add metadata
identical(rownames(t.fil.stool.taxa.sum.genus), rownames(stool.md)) # TRUE
stool.md$status_type <- paste0(stool.md$status,"_",stool.md$sample_type)
stool.genus.md = as.data.frame(cbind("status" = stool.md$status, t.fil.stool.taxa.sum.genus))
names(stool.genus.md)

# Group and plot

stool.genus.md.avg = stool.genus.md %>%
  group_by(status) %>%
  summarise_all(list(mean)) %>%
  gather(key = "Taxa", value = "Avg_Relabund","1174-901-12":"Xylophilus") 

#remove taxa with Rel Abu less than 0.1%
stool.genus.md.avg = stool.genus.md.avg[stool.genus.md.avg$Avg_Relabund > 0.001,]

mypalette2  <- c("#40004b","#ffffbf","#762a83","#fdbf6f","#9970ab","#ccebc5","#c2a5cf","#e7d4e8","#b8e186","#bebada",
                 "#de77ae","#f46d43","#f7f7f7","#d9f0d3","#a6dba0","#4d9221","#2166ac","#5aae61","#1b7837","#d73027",
                 "#00441b","#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#b2182b","#de77ae",
                 "#bc80bd","#f5f5f5","#c7eae5","#80cdc1","#35978f","#003c30",
                 "#8e0152","#fb8072","#c51b7d","#f1b6da","#fde0ef","#fdb462","#b3de69","#7fbc41","#276419",
                 "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#ff7f00","#cab2d6","#6a3d9a","#ffff99",
                 "#8dd3c7","#ffffb3","#80b1d3","#fdb462","#01665e","#fccde5","#d9d9d9",
                 "#40004b","#762a83","#fdbf6f","#9970ab","#ccebc5","#c2a5cf","#e7d4e8","#b8e186","#bebada",
                 "#de77ae","#f7f7f7","#d9f0d3","#a6dba0","#4d9221","#5aae61","#1b7837"
                 ,"#00441b","#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#bc80bd","#f5f5f5","#c7eae5","#80cdc1","#35978f","#003c30",
                 "#8e0152","#fb8072","#c51b7d","#f1b6da","#fde0ef","#fdb462","#b3de69","#7fbc41","#276419",
                 "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#ff7f00","#cab2d6","#6a3d9a","#ffff99",
                 "#8dd3c7","#ffffb3","#80b1d3","#fdb462","#01665e","#fccde5","#d9d9d9")


ggplot(stool.genus.md.avg, aes(x = status, y = Avg_Relabund, fill = Taxa)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mypalette2, guide = guide_legend(nrow=30)) +
  theme(legend.position="right", text = element_text(size = 12),  axis.text.x=element_text(angle=90))
ggsave("stool genus .pdf",height = 9, width = 10)


######## find genera statisticly significan for Hypertension Chow vs S-Chow
# Transpose and add metadata
# First check for equality
sapply(stool.taxa.sum, function(o){
  identical(rownames(stool.md), names(o))
}) # TRUE
# Transpose
t.stool.taxa.sum = lapply(stool.taxa.sum, function(v){
  as.data.frame(t(v))
})
# Add metadata to split samples by groups
t.stool.taxa.sum.md = lapply(t.stool.taxa.sum, function(h){
  as.data.frame(cbind(stool.md$status, h))
})
# Split into different lists
stool.split = lapply(t.stool.taxa.sum.md, function(p){
  split(p,stool.md$status)
})
# Remove empty dataframes
sapply(stool.split, function(v){
  print(names(v))
})
fil.stool.split = lapply(stool.split, function(w){
  w[names(w) %in% c("PANDAS", "CNTL")]
})
# Remove metadata column
fil.stool.split.cleaned = lapply(fil.stool.split, function(w){
  lapply(w, function(l){
    l[,-1]
  })
})
# Mann-Whitney U tests for all taxa at each level
# Old after ABX vs young after ABX
stool.mw = lapply(fil.stool.split.cleaned, function(d){
  mapply(function(a, z){ # Run wilcoxon by column
    as.data.frame(wilcox.test(a, z, exact = F)$p.value)
  }, d$"PANDAS", d$"CNTL")
  #print(d$POST_CONTROL[1])
})


sig.stool.mw = lapply(stool.mw, function(c){
  c[c<=0.05]
})

# Convert nested lists (lists in lists) to character vectors
write.sig.stool.mw = lapply(sig.stool.mw, function(r){
  w = as.character(r) # Convert lists to character vectors
  names(w) = names(r) # Make into named character vectors
  w
})

# Check  - should be character vectors
sapply(write.sig.stool.mw, class)

# Filter out NAs and convert character vector to numeric vector

w = lapply(write.sig.stool.mw, function(w){
  w[!is.na(names(w))]
})


# Clean names 
cleaning = lapply(w, function(r){
  names(r) = gsub(".wilcox.test(a, z, exact = F)$p.value", "", names(r), fixed = T)
  r
})

# Make naming vector for files
cleaning.names = c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# check
length(cleaning.names) == length(cleaning) # TRUE

# Save lists to separate text files
mapply(function(d, v){
  write.table(d, v, quote = F, sep = "\t")
}, cleaning, paste0(cleaning.names, "sig_p_vals_stool.txt")) # Ignore what prints to screen, can do list.files(pattern = "*_p_vals_Rivera_dada2.txt") to see new files 















#### Sample 6 compare A, B, C colelctions #########
X.noNULL.md$SampleRepetition <- paste(X.noNULL.md$Collection.day, X.noNULL.md$Sample, X.noNULL.md$Repetition)
X.noNULL.md.genus.sample6day12 = as.data.frame(cbind("SampleRepetition" = X.noNULL.md$SampleRepetition, 
                                                     "lactic.acid" = X.noNULL.md$lactic.acid, 
                                                     "acetic.acid" = X.noNULL.md$acetic.acid,
                                                     "propionic.acid" = X.noNULL.md$propionic.acid,
                                                     "pH" = X.noNULL.md$pH,
                                                     "butiric.acid"  = X.noNULL.md$butyric.acid,
                                                     t.fil.XnoNULL.genus))





X.noNULL.md.genus.avg.sample6day12 = X.noNULL.md.genus.sample6day12 %>%
  group_by(SampleRepetition) %>%
  summarise_all(funs(mean)) %>%
  gather(key = "Taxa", value = "Avg_Relabund","A2":"Zoogloea") 

f.X.noNULL.md.genus.avg.sample6day12 = X.noNULL.md.genus.avg.sample6day12[X.noNULL.md.genus.avg.sample6day12$Avg_Relabund > 0.001,]

pdf("Taxa genus Sample6 day12 repetition 0.001.pdf", width = 6, height = 6)
ggplot(f.X.noNULL.md.genus.avg.sample6day12, aes(x = SampleRepetition, y = Avg_Relabund)) + 
  geom_bar(aes(y= Avg_Relabund, fill= Taxa), stat = "identity") +
  scale_fill_manual(values = mypalette3, guide = guide_legend(nrow=20)) +
  theme(legend.position="right", text = element_text(size = 12),  axis.text.x=element_text(angle=90))+
  geom_point(data = f.X.noNULL.md.genus.avg.sample6day12, aes(x = SampleRepetition, y=pH/8), colour = "red", size = 4)+
  scale_y_continuous(sec.axis = sec_axis(~.*8, name = "pH", breaks = seq(0, 8, 0.5)))+
  scale_x_discrete(limits=c("12 6 A","12 6 B","12 6 C"))
dev.off()

pdf("Taxa genus Sample6 day15 repetition 0.001.pdf", width = 6, height = 6)
ggplot(f.X.noNULL.md.genus.avg.sample6day12, aes(x = SampleRepetition, y = Avg_Relabund)) + 
  geom_bar(aes(y= Avg_Relabund, fill= Taxa), stat = "identity") +
  scale_fill_manual(values = mypalette3, guide = guide_legend(nrow=20)) +
  theme(legend.position="right", text = element_text(size = 12),  axis.text.x=element_text(angle=90))+
  geom_point(data = f.X.noNULL.md.genus.avg.sample6day12, aes(x = SampleRepetition, y=pH/8), colour = "red", size = 4)+
  scale_y_continuous(sec.axis = sec_axis(~.*8, name = "pH", breaks = seq(0, 8, 0.5)))+
  scale_x_discrete(limits=c("15 6 A","15 6 B","15 6 C"))
dev.off()

pdf("Taxa genus Sample6 day18 repetition 0.001.pdf", width = 6, height = 6)
ggplot(f.X.noNULL.md.genus.avg.sample6day12, aes(x = SampleRepetition, y = Avg_Relabund)) + 
  geom_bar(aes(y= Avg_Relabund, fill= Taxa), stat = "identity") +
  scale_fill_manual(values = mypalette3, guide = guide_legend(nrow=20)) +
  theme(legend.position="right", text = element_text(size = 12),  axis.text.x=element_text(angle=90))+
  geom_point(data = f.X.noNULL.md.genus.avg.sample6day12, aes(x = SampleRepetition, y=pH/8), colour = "red", size = 4)+
  scale_y_continuous(sec.axis = sec_axis(~.*8, name = "pH", breaks = seq(0, 8, 0.5)))+
  scale_x_discrete(limits=c("18 6 A","18 6 B","18 6 C"))
dev.off()

#to remove sample 0 0
X6 <- droplevels(X.noNULL.md6[!X.noNULL.md6$Sample %in% "0",])
ggplot(X6, aes(x=Collection.day))+
  geom_point(aes(y=lactic.acid), colour = "red", size = 4)+
  geom_point(aes(y=acetic.acid), colour = "blue", size = 4)+
  geom_point(aes(y=propionic.acid), colour = "black", size = 4)+
  geom_point(aes(y=butyric.acid), colour = "darkgreen", size = 4)+
  theme(legend.position = "right")+
  facet_grid(Repetition ~ Collection.day, labeller = label_both)
ggsave("acids.pdf", width = 10, height =8,dpi=300)



############################### HEATMAP ##########################################

library(Heatplus)
library(RColorBrewer)
library(gplots)
library(vegan3d)

# on calulated previously B-C dissmiliratiy matrix on the full X dataset micro.dist.X.noNULL and XnoNULL.genus (for genus levle)
micro.dist.X.noNULL
rownames(micro.dist.X.noNULL)
names(micro.dist.X.noNULL)

# Transpose and convert to dataframe
t.XnoNULL.genus = as.data.frame(t(XnoNULL.genus))
t.XnoNULL.genus[1:3,1:4]
rownames(t.XnoNULL.genus)
rownames(X.noNULL.md)

names(X.noNULL.md)

X.noNULL.md$SampleRepetition2 <- paste(X.noNULL.md$Sample, X.noNULL.md$Collection.day, X.noNULL.md$Repetition)

rownames(t.XnoNULL.genus) == rownames(X.noNULL.md)
rownames(t.XnoNULL.genus) = X.noNULL.md$SampleRepetition2
rownames(t.XnoNULL.genus) ==X.noNULL.md$SampleRepetition2

rownames(X.noNULL.md) <- X.noNULL.md$SampleRepetition2
rownames(X.noNULL.md) == rownames(t.XnoNULL.genus)


# colorRampPalette is in the RColorBrewer package.  This creates a colour palette that shades from light yellow to red in RGB space with 100 unique colours
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)

heatmap(as.matrix(t.XnoNULL.genus), Rowv = NA, Colv = NA, col = scaleyellowred)

# determine the maximum relative abundance for each column
maxab <- apply(t.XnoNULL.genus, 2, max)
head(maxab)
ncol(t.XnoNULL.genus) #254

# remove the genera with less than 0.1% as their maximum relative abundance
n1 <- names(which(maxab < 0.001))
t.XnoNULL.genus1 <- t.XnoNULL.genus[, -which(names(t.XnoNULL.genus) %in% n1)]
head(t.XnoNULL.genus1)
ncol(t.XnoNULL.genus1) #31 , it means that most of the taxa occur at very low relative abundance (223 taxa with rel abun lower than 0.1%)

#to remove sample "0 0 0"
t.XnoNULL.noInoculum.genus1 <- droplevels(t.XnoNULL.genus1[!rownames(t.XnoNULL.genus1) %in% "0 0 0",])

heatmap(as.matrix(t.XnoNULL.genus1), Rowv = NA, Colv = NA, col = scaleyellowred)


# Do average linkage hierarchical clustering. Other options are 'complete' or 'single'. 
# You'll need to choose the one that best fits the needs of your situation and your data.
row.clus <- hclust(micro.dist.X.noNULL, "aver")

# make the heatmap with Rowv = as.dendrogram(row.clus)
heatmap(as.matrix(t.XnoNULL.genus1), Rowv = as.dendrogram(row.clus), Colv = NA, col = scaleyellowred, margins = c(10, 3))

# You can also add a column dendrogram to cluster the genera that occur more often together. 
# Note that this one must be done on the same dataset that is used in the Heatmap (i.e. reduced number of genera).

# you have to transpose the dataset to get the genera as rows
data.dist.g <- vegdist(t(t.XnoNULL.genus1), method = "bray")
data.dist.g1 <- vegdist(t(t.XnoNULL.noInoculum.genus1), method = "bray")
col.clus <- hclust(data.dist.g, "aver")
col.clus1 <- hclust(data.dist.g1, "aver")


# make the heatmap with Rowv = as.dendrogram(row.clus)
heatmap(as.matrix(t.XnoNULL.genus1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, margins = c(10, 3))
heatmap.2(as.matrix(t.XnoNULL.genus1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, margins = c(10, 3))

heatmap.2(as.matrix(t.XnoNULL.genus1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred,
          margins = c(11, 5), trace = "none", density.info = "none", xlab = "genera", ylab = "Samples", 
          main = "Heatmap example", lhei = c(2, 8)) # this makes the colour-key legend a little thinner

# the annHeatmap2 function needs to be wrapped in the plot function in order to display the results
dd<-c("lactic.acid","pH","Sample","acetic.acid","propionic.acid","butyric.acid")
ann.dat <- X.noNULL.md[dd]


# the annHeatmap2 function needs to be wrapped in the plot function in order to display the results
dd1<-c("lactic.acid","pH","Sample","acetic.acid","propionic.acid","butyric.acid")
ann.dat1 <- X.noNULL.md[dd]

X.rar.otu.df1 <- X.rar.otu.df
rownames(X.rar.otu.df1) <- rownames(X.noNULL.md)
rownames(X.rar.otu.df1) == rownames(X.noNULL.md)
X.rar.otu.df1 <- droplevels(X.rar.otu.df1[!rownames(X.rar.otu.df1) %in% "0 0 0",])
X.noNULL.md12 <- droplevels(X.noNULL.md[!rownames(X.noNULL.md) %in% "0 0 0",])

micro.dist.X.noNULL1<-vegdist(X.rar.otu.df1, method="bray")
NMDS_X <- metaMDS(micro.dist.X.noNULL1, distance = "bray", k=2, autotransform = T)
stressplot(NMDS_X)

row.clus1 <- hclust(micro.dist.X.noNULL1, "aver")

pdf("Summary Heatmap.pdf", width = 15, height = 15)
plot(annHeatmap2(as.matrix(t.XnoNULL.genus1), col = colorRampPalette(c("lightblue", "red"), space = "rgb")(31), 
                 breaks = 30, dendrogram = list(Row = list(dendro = as.dendrogram(row.clus)), Col = list(dendro = as.dendrogram(col.clus))), 
                 legend = 2, labels = list(Col = list(nrow = 12)),
                 ann = list(Row = list(data = ann.dat)), cluster=list(cuth=0.5, label=c("Sample4","sss","erf", "ddd"))))
dev.off()


pdf("Summary Heatmap no inoculum.pdf", width = 15, height = 15)
plot(annHeatmap2(as.matrix(t.XnoNULL.noInoculum.genus1), col = colorRampPalette(c("lightblue", "red"),space = "rgb")(31),  
                 breaks = 30, dendrogram = list(Row = list(dendro = as.dendrogram(row.clus1)), Col = list(status="hide")), 
                 legend = 3, labels = list(Col = list(nrow = 18, cex = 1.5), Row = list(nrow= 3, cex=1.2)),
                 ann = list(Row = list(data = ann.dat1)), cluster=list(cuth=0.5, label=c("Sample4","sss","erf", "ddd"))))
dev.off()



plot(annHeatmap2(as.matrix(t.XnoNULL.noInoculum.genus1), col = colorRampPalette(c("lightblue", "red"),space = "rgb")(31),  
                 breaks = 30, dendrogram = list(Row = list(dendro = as.dendrogram(row.clus1)),Col = list(status="hide")), 
                 legend = 3, labels = list(Col = list(nrow = 18, cex = 1.5), Row = list(nrow= 3, cex=1.2)),
                 ann = list(Row = list(data = ann.dat1)), cluster=list(cuth=0.5, label=c("Sample4","sss","erf", "ddd"))))


############## changes 


# the annHeatmap2 function needs to be wrapped in the plot function in order to display the results
dd2<-c("pH","lactic.acid","acetic.acid","propionic.acid","butyric.acid","ethanol","Sample")
ann.dat2 <- X.noNULL.md[dd2]


pdf("Summary Heatmap no inoculum new.pdf", width = 15, height = 15)
plot(annHeatmap2(as.matrix(t.XnoNULL.noInoculum.genus1), col = colorRampPalette(c("lightblue", "red"),space = "rgb")(31),  
                 breaks = 30, dendrogram = list(Row = list(dendro = as.dendrogram(row.clus1)),Col = list(status="hide")), 
                 legend = FALSE, labels = list(Col = list(nrow = 18, cex = 1.5), Row = list(nrow= 3, cex=1.2)),
                 ann = list(Row = list(data = ann.dat2)), cluster=list(cuth=0.5)))
dev.off()


ania1<-annHeatmap2(as.matrix(t.XnoNULL.noInoculum.genus1), col = colorRampPalette(c("lightblue", "red"),space = "rgb")(31),  
                   breaks = 30, dendrogram = list(Row = list(dendro = as.dendrogram(row.clus1)),Col = list(status="hide")), 
                   legend = FALSE, labels = list(Col = list(nrow = 18, cex = 1.5), Row = list(nrow= 3, cex=1.2)),
                   ann = list(Row = list(data = ann.dat2)), cluster=list(cuth=0.5))
plot(ania1)


############################## T E S T ############## https://rdrr.io/cran/vegan3d/man/ordiplot3d.html
### Default 'ordiplot3d'
# data(dune, dune.env)
ord <- cca(t.XnoNULL.noInoculum.genus1 ~ butyric.acid + lactic.acid + acetic.acid, X.noNULL.md12)
ordiplot3d(ord)
### A boxed 'pin' version
ordiplot3d(ord, type = "h")
### More user control
pl <- ordiplot3d(ord, scaling = "symmetric", angle=15, type="n")
points(pl, "points", pch=16, col="red", cex = 0.7)
### identify(pl, "arrows", col="blue") would put labels in better positions
text(pl, "arrows", col="blue", pos=3)
text(pl, "centroids", col="blue", pos=1, cex = 1)
### Add species using xyz.convert function returned by ordiplot3d
sp <- scores(ord, choices=1:3, display="species", scaling="symmetric")
text(pl$xyz.convert(sp), rownames(sp), cex=0.9, xpd=TRUE)
### Two ways of adding fitted variables to ordination plots
ord <- cca(t.XnoNULL.noInoculum.genus1)
ef <- envfit(ord ~ butyric.acid + lactic.acid + acetic.acid, X.noNULL.md12, choices = 1:5)
### 1. use argument 'envfit'
ordiplot3d(ord, envfit = ef)
### 2. use returned envfit.convert function for better user control
pl3 <- ordiplot3d(ord)
plot(pl3$envfit.convert(ef), at = pl3$origin)
### envfit.convert() also handles different 'choices' of axes
pl3 <- with(X.noNULL.md12, ordiplot3d(ord, col = Sample, pch=16))
plot(pl3$envfit.convert(ef), at = pl3$origin)
### vegan::ordiXXXX functions can add items to the plot
ord <- cca(t.XnoNULL.noInoculum.genus1)
pl4 <- with(X.noNULL.md12, ordiplot3d(ord, col = Sample, pch=16))
with(X.noNULL.md12, ordiellipse(pl4, lactic.acid, draw = "poly", col = 1:4,
                                alpha = 60))
with(X.noNULL.md12, ordispider(pl4, lactic.acid, col = 1:4, label = TRUE))


######    TEST #####   https://rpubs.com/collnell/manova

NMDS1 <- NMDS_X$points[,1]
NMDS2 <- NMDS_X$points[,2]
NMDS.plot <- cbind(X.noNULL.md12, NMDS1, NMDS2)
p12<-ggplot(NMDS.plot, aes(NMDS1, NMDS2, color=Sample))+
  geom_point(position=position_jitter(.1), shape=3)+##separates overlapping points
  stat_ellipse(type='t',size =1)+ ##draws 95% confidence interval ellipses
  theme_minimal()
p12


#plot ordination withlactic acid 

p12<-ggplot(NMDS.plot, aes(NMDS1, NMDS2, color=Sample))+
  #geom_point(position=position_jitter(.1), shape=3)+##separates overlapping points
  stat_ellipse(type='t',size =1)+ ##draws 95% confidence interval ellipses
  theme_minimal()+
  geom_text(data=NMDS.plot,aes(NMDS1, NMDS2, label=lactic.acid), position=position_jitter(.35))+
  annotate("text", x=min(NMDS1), y=min(NMDS2), label=paste('Stress =',round(NMDS_X$stress,3)))

p12


# Fit vectors to ordination
dat12 <- X.noNULL.md12[dd2]
fit<-envfit(NMDS_X, dat12)
arrow<-data.frame(fit$vectors$arrows,R = fit$vectors$r, P = fit$vectors$pvals)
arrow$FG <- rownames(arrow)
arrow.p<-filter(arrow, P <= 0.05)

p<-ggplot(data=NMDS.plot, aes(NMDS1, NMDS2))+
  geom_point(data=NMDS.plot, aes(NMDS1, NMDS2, color=Sample),position=position_jitter(.1))+##separates overlapping points
  stat_ellipse(aes(fill=Sample), alpha=.2,type='t',size =1, geom="polygon")+ ##changes shading on ellipses
  theme_minimal()+
  geom_segment(data=arrow.p, aes(x=0, y=0, xend=NMDS1, yend=NMDS2, label=FG, lty=FG), arrow=arrow(length=unit(.2, "cm")*arrow.p$R)) ##add arrows (scaled by R-squared value)

p


p.grad <- ggplot(ordisurf(NMDS_X, dat12[,'lactic.acid'], bubble=TRUE)) + geom_point(data=NMDS.plot, aes(NMDS1, NMDS2, color=Sample))



# 2. merging samples based on the patient number  ----------------------------

seqtab.nochim.2 <- seqtab.nochim
md.2 <- md
md.2$pt_type <- paste0(md.2$patientID,".",md.2$sample_type)

rownames(seqtab.nochim.2) == rownames(md.2)
rownames(seqtab.nochim.2) <- md.2$pt_type

seqtab.nochim.2.sum <- rowsum(seqtab.nochim.2, rownames(seqtab.nochim.2))
rownames(seqtab.nochim.2.sum)
md.2.no_duplicates <- md.2 %>% distinct(pt_type, .keep_all = T) # to remove duplicated samples
nrow(md.2.no_duplicates)
nrow(md.2)
rownames(md.2.no_duplicates) <- md.2.no_duplicates$pt_type

rownames(md.2.no_duplicates) == rownames(seqtab.nochim.2.sum)
md.2.no_duplicates <- md.2.no_duplicates[order(rownames(md.2.no_duplicates)),]
seqtab.nochim.2.sum <- seqtab.nochim.2.sum[order(rownames(seqtab.nochim.2.sum)),]
rownames(md.2.no_duplicates) == rownames(seqtab.nochim.2.sum)

md2 <- droplevels(md.2.no_duplicates[md.2.no_duplicates$sample_type == "STOOL",])
table(md2$status)
seqtab.nochim.2.sum.stool <- seqtab.nochim.2.sum[row.names(seqtab.nochim.2.sum) %in% row.names(md2),]
rownames(seqtab.nochim.2.sum.stool) == row.names(md2)


stnc2 = as.data.frame(seqtab.nochim.2.sum.stool)

# Visualize and save sequence counts that have made it through quality control
sort(rowSums(stnc2))

ggplot(stnc2, aes(rowSums(stnc2))) + geom_histogram() 


# Rarefy for nasal samples
sort(rowSums(stnc2))
stool.sum.rar = rrarefy(stnc2, 38385)
stool.sum.rar.df = as.data.frame(stool.sum.rar)
rowSums(stool.sum.rar.df)
# Remove unrarefied samples (should lose 0)
stool.sum.asv = stool.sum.rar.df[rowSums(stool.sum.rar.df) == 38385,]
dim(stool.sum.rar.df)[1] - dim(stool.sum.asv)[1] # 4 samples were removed, both are blanks 


############### ALPHA DIVERISTY on rarefied otu (data frame) #############


md2$Richness <- specnumber(stool.sum.asv)
md2$Shannon <- diversity(stool.sum.asv, index="shannon")
md2$Simpson <- diversity(stool.sum.asv, index="simpson")
colnames(md2)

#### Plot Alpha Diversity 

ggplot(md2, aes(x=status, y=Richness, fill=status)) +
  geom_violin()+
  geom_jitter(position=position_jitter(0.2), size=4) +
  geom_boxplot(alpha=0.4) +
  xlab("Patient status") +
  ylab("Richness") +
  theme_bw()+
  theme(legend.position="none", text = element_text(size=16), axis.text.x=element_text(angle=0, hjust=0.5))+
  ggtitle("Stool combined")
ggsave("stool combined Richness.pdf",width = 6, height = 4)

wilcox.test(md2$Richness ~ md2$status, md2) #W = 78, p-value = 0.007936

ggplot(md2, aes(x=status, y=Shannon, fill=status)) +
  geom_violin()+
  geom_jitter(position=position_jitter(0.2), size=4) +
  geom_boxplot(alpha=0.4) +
  xlab("Patient status") +
  ylab("Shannon") +
  theme_bw()+
  theme(legend.position="none", text = element_text(size=16), axis.text.x=element_text(angle=0, hjust=0.5))+
  ggtitle("Stool combined")
ggsave("stool combined Shannon",width = 6, height = 4)

wilcox.test(md2$Shannon ~ md2$status, md2) #W = 78, p-value = 0.02793

ggplot(md2, aes(x=status, y=Simpson, fill=status)) +
  geom_violin()+
  geom_jitter(position=position_jitter(0.2), size=4) +
  geom_boxplot(alpha=0.4) +
  xlab("Patient status") +
  ylab("Simpson") +
  theme_bw()+
  theme(legend.position="none", text = element_text(size=16), axis.text.x=element_text(angle=0, hjust=0.5))+
  ggtitle("Stool combined")
ggsave("stool combined Simpson",width = 6, height = 4)

wilcox.test(Simpson ~ status, md2) #W = 78, p-value = 0.05348


#################### BETA DIVERSITY ######################


#NMDS ordination for nasal
stool.sum.micro.dist<-vegdist(stool.sum.asv, method="bray")
stool.sum.micro.dist.nmds<-metaMDS(stool.sum.micro.dist, k=2)
stool.sum.micro.dist.nmds$stress #0.161
stressplot(stool.sum.micro.dist.nmds)
md2$NMDS01<-stool.sum.micro.dist.nmds$points[,1]
md2$NMDS02<-stool.sum.micro.dist.nmds$points[,2]


#ordination plots
find_hull <- function(df) df[chull(df$NMDS01, df$NMDS02),]

stool.sum.hull <- ddply(md2, "status", find_hull)


ggplot(md2, aes(NMDS01, NMDS02))+
  geom_point(aes(fill=status), size=6, alpha=1, pch=21, col="black")+
  #geom_polygon(data = stool.sum.hull, aes(colour=status, fill=status), alpha = 0.5)+
  stat_ellipse(aes(color=status), type = "norm", linetype = 2) +
  stat_ellipse(aes(color=status),type = "t")+
  theme(legend.position="bottom", text = element_text(size = 18))+
  ggtitle("Stool combined")
ggsave("nasal BC NMDS.pdf", height = 5, width = 7)

adonis2(stool.sum.asv ~ status, md2, permutations = 999,distance = "bray")
