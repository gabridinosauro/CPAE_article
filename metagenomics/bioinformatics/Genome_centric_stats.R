### MAGs ####
# MAGs analysis 
# Header --------

#MAGs retrival from the CPAE project
#Gabriele Schiro

##Libraries ----
library("ggplot2")
library(ggtree)
library(stringr)
library(ggtreeExtra)
library("vegan")

`%nin%` = Negate(`%in%`)
##Files -------
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Meta_gen_Daniel/MAGs")
summary = read.delim("reassembled_bins.stats")
tree <- read.tree("RAxML_bestTree.reassembled_bins_refined.tre")
taxonomy <- read.delim("gtdbtk.bac120.summary.tsv")
map = read.csv("map.csv")
bin_abundance = read.delim("bin_abundance_table.tab")
Ko_mapping = read.csv("KO_table.csv", row.names = 1)
Ko_mapping = Ko_mapping[,order(ncol(Ko_mapping):1)]
count_table_genes = read.csv("gene_table.csv", row.names = 1)
count_table_genes = count_table_genes[,order(ncol(count_table_genes):1)]
taxotab= read.csv("metaphlan_results.csv")
summary$bin = gsub("^([^.]*.[^.]*)..*$", "\\1", summary$bin)
tree$tip.label =  gsub("^([^.]*.[^.]*)..*$", "\\1", tree$tip.label)
taxonomy$user_genome = gsub("^([^.]*.[^.]*)..*$", "\\1", taxonomy$user_genome)
summary$High_quality = ifelse(summary$completeness >95 & summary$contamination < 5, "High quality", "Low quality")
ggplot(summary,aes(completeness, contamination, color = High_quality)) + geom_point() + theme_bw()







###Fix taxonomy
classification = data.frame(str_split_fixed(taxonomy$classification, ";", 7))
rownames(classification) =taxonomy$user_genome
colnames(classification) =c("Kingdom", "Phylum","Class" ,"Order", "Family", "Genus", "Species")
classification = classification[-which(rownames(classification) %nin% tree$tip.label),]
classification = classification[match(tree$tip.label,rownames(classification)),]
rownames(classification) == tree$tip.label #great
classification = data.frame(apply(classification, 2, function(x) substring(x,4)))
classification$Species=ifelse(classification$Genus == "Collinsella", "Collinsella",classification$Species)
filtered_text <- classification$Species[!grepl("_[ABC]$", classification$Species)]
length(table(filtered_text))
filtered_text <-  sub("_[ABCQ]$", "", classification$Genus)
length(table(filtered_text))
table(classification$Phylum)


### Check if names are the same
unique(summary$bin == tree$tip.label) #FALSE, oh no (shorter length)
unique(summary$bin %in% tree$tip.label) #TRUE false
unique(tree$tip.label %in% summary$bin)#TRUE oh yes
summary$bin[which(summary$bin  %nin% tree$tip.label)] #somehow one bin was not included in the tree...remove it
summary = summary[-which(summary$bin  %nin% tree$tip.label),]
summary = summary[match(tree$tip.label,summary$bin),]
summary$bin == tree$tip.label #great
summary$`size MB` = summary$size/1000

# Plot tree ----------------------------
tree$tip.label  = classification$Species
tree$tip.label = ifelse(classification$Species == "", classification$Genus, classification$Species)
tree$tip.color = c(classification$Phylum)
p <- ggtree(tree, branch.length = "none",layout = "circular") + 
    theme(legend.position='none')
p$data$phylum = NA
p$data$phylum[1:157] = classification$Phylum
p$data
summary$ID=p$data$label[1:157]
p + geom_tiplab(align= T, linetype=NA, size=2.5,aes(color = phylum)) + scale_colour_brewer(palette = "Paired") + geom_treescale()
#+geom_tippoint(aes(color = phylum))

#ggsave("tree_taxo.pdf")
#+ geom_fruit(data=summary,
#             geom=geom_col,                                                                                           
#             mapping=aes(y=ID, x=`size MB`), 
#             pwidth=0.5,
#             offset = 1,
#             axis.params=list(
#             axis="x", # add axis text of the layer.
#             text.angle=-45, # the text size of axis.
#             hjust=0,text.size = 2 ,
#             #title = "Size (MB)"# adjust the horizontal position of text of axis.
#            grid.params=list() # add the grid line of the external bar plot.
#) 

## Correlations with diversity
unique(bin_abundance$Genomic.bins == rownames(classification)) #FALSE, TRUE oh no (shorter length)
unique(bin_abundance$Genomic.bins %in% rownames(classification)) # TRUE FALSE
unique(rownames(classification) %in% bin_abundance$Genomic.bins) #TRUE 
bin_abundance$Genomic.bins[which(bin_abundance$Genomic.bins %nin% rownames(classification))] #somehow one bin was not included in the tree...remove it
bin_abundance = bin_abundance[-which(bin_abundance$Genomic.bins %nin% rownames(classification)),]
bin_abundance = bin_abundance[match(rownames(classification),bin_abundance$Genomic.bins),]
unique(bin_abundance$Genomic.bins == rownames(classification)) # YEY!
rownames(bin_abundance) = bin_abundance$Genomic.bins 
bin_abundance = bin_abundance[,-1]
colnames(bin_abundance) = substring(colnames(bin_abundance), 1, 4) #OK
bin_abundance = bin_abundance[,match(map$Sample_ID,substring(colnames(bin_abundance), 2, 4))] #OK

# Select those bins that are present in every sample
bin_abundance_core = bin_abundance[apply(bin_abundance!=0, 1, all),]
bin_abundance_core$species = classification[match(rownames(bin_abundance_core),rownames(classification)),]




### ASVs to MAGs 
## First I will blast the ASVs into the MAgs
### First step is to modify all headers of the files by  going to the folder 
# /Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/MAGs_genomes_renamed

#and running the follwoing loop

#for input_file in *.fa; do
#if [ -f "$input_file" ]; then
#python modify_headers.py "$input_file"
#fi
#done


## I have a conda environment with blast installed.
# I first cat all fasta files into one

# cat *head.fasta > all_genomes.fasta

# create a database:
# makeblastdb -in all_genomes.fasta -dbtype nucl -out database


# create.a fasta file with my ASVs
#sequences = data.frame(names =taxa_stool$names, ASV = rownames(taxa))
#library(seqinr)
#write.fasta(as.list(sequences$ASV), names = sequences$names, file.out = "/Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/MAGs_genomes_renamed/ASV.fasta")
# Blast it this way blastn -query ASV.fasta -db database -out results_blast.txt -perc_identity 100 -outfmt "6 qseqid sseqid length"
### 
#Blast_results = read.delim("/Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/MAGs_genomes_renamed/results_blast.txt", header = FALSE)
#Blast_results$bin = result <- sub(".*?(bin.*)", "\\1", Blast_results$V2)
#Blast_results = Blast_results[Blast_results$V3 >150,]
#Blast_results$bin_taxonomy = taxonomy$classification[match(Blast_results$bin,taxonomy$user_genome)]
#Blast_results$ASV_tax = taxa_stool$names[match(Blast_results$V1,rownames(taxa_stool))]
#stool.md$new_names = ifelse(stool.md$status == "CPAE", "P", "C")
#stool.md$new_names = paste(stool.md$new_names, substr(stool.md$patientID, nchar(stool.md$patientID) - 3, nchar(stool.md$patientID)), sep = "-")
#stool.md_16s = stool.md[stool.md$new_names %in% map$Patient,]












### Run Omixer on isolates
### Gut metabolites 
library(omixerRpm)
setwd("/Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/new_analyses")
Ko_mapping_genomes = readRDS("KO_counts.RDS")
Ko_mapping_genomes = cbind(rownames(Ko_mapping_genomes), Ko_mapping_genomes)
colnames(Ko_mapping_genomes)[1] = "entry"
db <- loadDefaultDB()
mods = rpm(Ko_mapping_genomes, normalize.by.length = TRUE)
coverage_gmm = asDataFrame(mods, "coverage")
db <- loadDB("GBMs.v1.0")
mods = rpm(Ko_mapping_genomes,  module.db =  loadDB("GBMs.v1.0"))
mods@coverage
coverage_gbm <- asDataFrame(mods, "coverage")
t_coverage_gbm = data.frame(t(coverage_gbm))





### Find Isolates that go up or down.
### Use the results from coverM instead of metaWRAP.
cover_M = readRDS("coverm_rel_ab.RDS")
colnames(cover_M) = gsub("_S.*", "", colnames(cover_M))
cover_M =  cover_M[,match(map$Sample_ID,substring(colnames(cover_M), 2, 4))] #OK
source("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")
WC = wilcox_calc(t(cover_M), map$Diagnosis)

map$Diagnosis = factor(map$Diagnosis, levels = c("CPAE", "CTRL"))
rownames(map) = paste("X",map$Sample_ID, sep = "")
library("Maaslin2")
fit_data2 = Maaslin2(input_data     = t(cover_M), 
                     input_metadata = map, 
                     min_prevalence = 0.4, 
                     #normalization  = "CSS",
                     max_significance = 0.5,
                     transform = "LOG",
                     output         = "demo_output2", 
                     fixed_effects  = "Diagnosis",
                     reference      = c("Diagnosis,CTRL"),
                     #random_effects = "Patient",
                     plot_heatmap = FALSE,
                     plot_scatter = FALSE
                     
)
res = fit_data2$results
res$feature = sub("\\.[^.]*$", "", res$feature)
res$coef = -1 *  res$coef
res$taxonomy = classification$Species[match(res$feature,rownames(classification))]
res$final = paste(res$feature, res$taxonomy, sep = "-")
#xlsx::write.xlsx(res, "/Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/supplementary_information/Supplementary_table3.xlsx")
res_sig = res[res$pval<0.05,]
relative_abundance = sort(rowSums(cover_M), decreasing = TRUE)
names(relative_abundance) =  sub("\\.[^.]*$", "", names(relative_abundance))
res_sig_relative_abundance = relative_abundance[names(relative_abundance) %in% res_sig$feature]




#geom_vline(xintercept = 0)
coverM_t = data.frame(t(cover_M))
colnames(coverM_t) = sub("\\.[^.]*$", "", colnames(coverM_t) )
coverM_t = coverM_t[,colnames(coverM_t) %in%res_sig$feature]
relative_abundance = sort(colSums(coverM_t), decreasing = TRUE)
coverM_t$status = map$Diagnosis
matched = res_sig$final[match(names(relative_abundance),res_sig$feature)]
res_sig$final = factor(res_sig$final, matched)

## Cumulative relative abundance
a = ggplot(res_sig) +
  geom_point(aes(x = coef, y = final)) +
  theme_bw() + ylab(NULL) +ggtitle("CPAE vs CONTROL") +
  theme(plot.title = element_text(size=9)) + geom_vline(xintercept = 0) +
  scale_y_discrete(limits=rev)
#+
#labs(col = "Phylum")
a
names(relative_abundance)






melted_ASV = reshape2::melt(coverM_t)
#melted_ASV$variable = sub("\\.[^.]*$", "", melted_ASV$variable)
melted_ASV$variable
matched =sub("^(.*?)-.*$", "\\1", matched)
matched
melted_ASV$variable = factor(melted_ASV$variable,levels = matched )
melted_ASV$value = melted_ASV$value / 100

b = ggplot(melted_ASV, aes(value,variable , color = status))+
  geom_boxplot(outlier.shape = NA) +theme_bw() +
#  geom_point(position = position_jitterdodge(), alpha=0.3) +
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_y_discrete(limits=rev) + xlab("Relative abundance")+ theme(legend.position = "none") + xlim(0,0.05)
b

library(ggpubr)
ggarrange(a,b, align = "h",widths = c(0.5,0.4))


coverage_gbm1 =  coverage_gbm[,3:ncol(coverage_gbm)]
colnames(coverage_gbm1) = sub("\\.[^.]*$", "", colnames(coverage_gbm1) )
rownames(coverage_gbm1) = coverage_gbm$Description
coverage_gbm1 = coverage_gbm1[,colnames(coverage_gbm1) %in% res_sig$feature]
coverage_gbm1[coverage_gbm1<1] = 0
coverage_gbm1 = coverage_gbm1[rowSums(coverage_gbm1)>0,]
melted_c = data.table::melt(as.matrix(coverage_gbm1))
melted_c$Modules = "GBM"
coverage_gmm1=  coverage_gmm[,3:ncol(coverage_gmm)]
colnames(coverage_gmm1) = sub("\\.[^.]*$", "", colnames(coverage_gmm1) )
rownames(coverage_gmm1) = coverage_gmm$Description
coverage_gmm1 = coverage_gmm1[,colnames(coverage_gmm1) %in% res_sig$feature]
coverage_gmm1[coverage_gmm1<1] = 0
coverage_gmm1 = coverage_gmm1[rowSums(coverage_gmm1)>0,]
melted_d = data.table::melt(as.matrix(coverage_gmm1))
melted_d$Modules = "GMM"
melted_f = rbind(melted_c, melted_d)

xlsx::write.xlsx(coverage_gbm, "/Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/supplementary_information/Supplementary_table4.xlsx",
                 sheetName = "GBM")
xlsx::write.xlsx(coverage_gmm, "/Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/supplementary_information/Supplementary_table4.xlsx",
                 sheetName = "GMM", append = TRUE)

other_plot_modules = read.csv("results_diffabb_modules.csv", row.names = 1)


melted_f$final = res_sig$final[match(melted_f$Var2,res_sig$feature )]
melted_f$value = ifelse(melted_f$value == 0, "Absent", "Present")
melted_f = melted_f[melted_f$Var1 %in% other_plot_modules$x,]
c = ggplot(melted_f, aes(x=final, y=Var1, shape = factor(value))) +
  geom_point(size = 3, color = "darkblue") +
  scale_shape_manual(values = c(1, 16)) +theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +xlab(NULL) +facet_grid(~Modules, scales = "free", space = "free" ) +
  labs(shape=NULL)  + theme(legend.position="top") + ylab(NULL) +coord_flip() +scale_x_discrete(limits=rev) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

c
ggarrange(a,b, c,align = "h",widths = c(0.5,0.3, 0.3), ncol = 3, nrow = 1)
library(patchwork)

# Arrange plots horizontally
a + b + c











# Calculate richness --------

#plot(specnumber(t(taxo_tab_species))~ordered_map$Number_of_reads)
map$richness = specnumber(t(cover_M))
map$Shannon_H = diversity(t(cover_M), index = "shannon")



a = ggplot(map, aes( y=richness, x=Diagnosis,fill = Diagnosis)) +  theme_bw() +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab(NULL)+ ylab("Species richness")  +
  geom_line(aes(group=Diagnosis), position = position_dodge(0.2), alpha = 0.1) +
  geom_point(aes(group=Diagnosis), position = position_dodge(0.2)) +
  scale_fill_manual(values=c('darkred', "darkblue"))+ 
  theme(legend.position = "none", axis.text.x = element_blank())
a



b = ggplot(map, aes( y=Shannon_H, x=Diagnosis,fill = Diagnosis)) +  theme_bw() +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab(NULL)+ ylab("Shannon H'") + 
  geom_line(aes(group=Diagnosis), position = position_dodge(0.2), alpha = 0.1) +
  geom_point(aes(group=Diagnosis), position = position_dodge(0.2)) +
  scale_fill_manual(values=c('darkred', "darkblue"))
b
ggarrange(a,b, ncol = 1, align = "v",heights = c(0.8,1.1), common.legend = TRUE )







library(metagenomeSeq)


#metaSeqObject = newMRexperiment(taxo_tab_species) #create a metagenomeseq experiment
#metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) ) #CSS transformation
#taxotab_CSS = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE)) 
species_bray = vegdist(t(cover_M))
species_nmds <- metaMDS(species_bray, k=2, try = 100)
stress= paste("stress = ", round(species_nmds$stress, digits = 2))
map$NMDS001 = species_nmds$points[,1]
map$NMDS002 = species_nmds$points[,2]
bac_nmds1<-ggplot(map, aes(x=NMDS001, y=NMDS002))+
  geom_point(aes(color= Diagnosis , shape = Diagnosis), size = 3)+
  scale_color_manual(values=c('darkred', "darkblue")) +
  theme_bw()+
  annotate("text",x=min(map$NMDS001)+ 0.5,y=max(map$NMDS002),hjust=1,label= stress) 
bac_nmds1
set.seed(4)
ado = adonis(species_bray~Diagnosis, map)
ado$aov.tab



### MAGS realtive abundance plot
coverM_t = data.frame(t(cover_M))
colnames(coverM_t) = sub("\\.[^.]*$", "", colnames(coverM_t))
relative_abundance = sort(colSums(coverM_t), decreasing = TRUE)
coverM_t$status = map$Diagnosis
melted_ASV = reshape2::melt(coverM_t)
order = as.data.frame(relative_abundance[-1])
#xlsx::write.xlsx(order,"relative_ab_MAGs.xlsx")
taxonomy$names = paste(taxonomy$user_genome ,sub(".*s__","", taxonomy$classification))
melted_ASV$names = taxonomy$names[match(melted_ASV$variable,taxonomy$user_genome)]
melted_ASV = melted_ASV[!melted_ASV$variable == "unmapped",]
order$names = taxonomy$names[match(rownames(order),taxonomy$user_genome)]
