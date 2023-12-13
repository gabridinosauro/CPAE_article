###Metagenomic analysis####
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(RColorBrewer)
library(vegan)
library(lme4)
library(ggforce)
library(concaveman)
library(ggtree)
library(microbiomeViz)
library(DESeq2)
library(pheatmap)
library(ggpubr)
library(Maaslin2)
library(dplyr)
library(metagenomeSeq)

mypalette2  <-c("#40004b","#2166ac","#ffffbf","#d73027","#762a83","#de77ae","#f46d43","#f7f7f7","#d9f0d3","#a6dba0",
                "#a6cee3","#1f78b4","#b2df8a","#33a02c", "#00441b","#543005","#8c510a","#bf812d","#dfc27d","#5aae61",
                "#f6e8c3","#b2182b","#de77ae","#bc80bd","#f5f5f5","#4d9221","#c7eae5","#80cdc1","#35978f","#003c30",
                "#8e0152","#fb8072","#c51b7d","#f1b6da","#fde0ef","#fdb462","#b3de69","#7fbc41","#276419","#a6cee3",
                "#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#ff7f00","#cab2d6","#6a3d9a","#ffff99")
### read mapping file 
setwd("/Users/gabri/OneDrive - University of Arizona/Meta_gen_Daniel/results")
map = read.csv("SampleSheetUsed.csv")
source("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")

###Read in the metaphlan table
taxo_tab = read.table("metaphlan/merged_abundance_table.txt", header = TRUE)
colnames(taxo_tab)[3:15] = substr(colnames(taxo_tab)[3:15] ,2, 4)
taxo_tab_phylum=taxo_tab[-grep("c__", taxo_tab$clade_name), ]


####Lets extract first viruses and make a boxplot
taxo_tab_virus = taxo_tab[taxo_tab$clade_name =="k__Viruses",]
rownames(taxo_tab_virus) = taxo_tab_virus$clade_name
taxo_tab_virus =data.frame(t(taxo_tab_virus[,-c(1:2)]))
taxo_tab_virus$treatment = map[match(rownames(taxo_tab_virus), map$Sample_ID), 8]
ggplot(taxo_tab_virus, aes(y = k__Viruses, x = treatment, fill = treatment))+ # for label --> aes(label = metadata$sample)
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = treatment))+ 
  xlab(NULL) + 
  ylab("Relative Abundance")+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Viral Abundance")
wilcox.test(k__Viruses~treatment, data = taxo_tab_virus)#W = 19, p-value = 0.8301















####Now,  let's do a Phylum stacked barplot
taxo_tab_species=taxo_tab[-grep("Viru", taxo_tab$clade_name), ]
taxo_tab_species=taxo_tab_species[grep("s__", taxo_tab_species$clade_name), ]
rownames(taxo_tab_species) = sub(".*s__", "", taxo_tab_species$clade_name)                 # Extract characters after pattern)
taxo_tab_species =data.frame(t(taxo_tab_species[,-c(1:2)]))
taxo_tab_species<-mutate_all(taxo_tab_species, function(x) as.numeric(as.character(x)))
taxo_tab_species1 = taxo_tab_species
summs = colSums(taxo_tab_species)
roww_summs = rowSums(taxo_tab_species)
taxo_tab_species$treatment = map[match(rownames(taxo_tab_species), map$Sample_ID), 8]
taxo_tab_species$samples = rownames(taxo_tab_species)






###Select the 20 most abundant species across the whole study
top_species = names(sort(summs, decreasing = TRUE))[1:25]
taxo_tab_species_melt = melt(taxo_tab_species, id.vars = c("treatment", "samples"))
taxo_tab_species_melt$variable = ifelse(taxo_tab_species_melt$variable %in% top_species, as.character(taxo_tab_species_melt$variable), "others")


# Stacked + percent
a = ggplot(taxo_tab_species_melt, aes(fill=variable, y=value, x= samples)) + 
  geom_bar(position="fill", stat="identity") +facet_wrap(~treatment, scales = "free") +
  theme_bw() + scale_fill_manual(values = mypalette2) + labs(fill='Species') + 
  ylab("Relative Abundance") +xlab("Sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
a



ordered_map = map[match(rownames(taxo_tab_species1),map$Sample_ID),]
### taxonomic ordination ###
species_bray <- vegdist(taxo_tab_species1, method="bray")
permanova = adonis(taxo_tab_species1 ~ Diagnosis, data = ordered_map)
permanova$aov.tab
species_nmds <- metaMDS(species_bray, k=2, try = 100) 
ordered_map$NMDS001 = species_nmds$points[,1]
ordered_map$NMDS002 = species_nmds$points[,2]
stress= paste("stress = ", round(species_nmds$stress, digits = 2))
ord1 = ggplot(ordered_map, aes(NMDS001, NMDS002, label = Sample_ID))+
  geom_point(aes(color=Diagnosis), size=4)+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=Diagnosis, color = Diagnosis, label = NA), con.type  = "none")+
  scale_color_manual(values=c('darkred', "darkblue"))+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  xlab("NMDS1") + 
  ylab("NMDS2")+
  theme_bw() +ggtitle("Species")+
  annotate("text",x=max(ordered_map$NMDS001),y=max(ordered_map$NMDS002),hjust=1,label= stress,  size = 3.5)+
  geom_text(hjust=1.2, vjust=0,  size = 3)
ord1

### Richness boxplot
ordered_map$taxo_richness  = specnumber(taxo_tab_species, MARGIN = 1)
ggplot(ordered_map, aes(y = taxo_richness, x = Diagnosis, fill = Diagnosis))+ # for label --> aes(label = metadata$sample)
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = Diagnosis))+ 
  xlab(NULL) + 
  ylab("Species richness")+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Richness")
summary(lm(taxo_richness ~ Diagnosis + Number_of_reads, data = ordered_map))









## Functional profiles ---------
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Meta_gen_Daniel/MAGs")
Ko_mapping = read.csv("KO_table_coverage_with_zeros.csv", row.names = 1)
Ko_mapping = Ko_mapping[,order(ncol(Ko_mapping):1)]
count_table_genes = read.csv("gene_table_coverage.csv", row.names = 1)
count_table_genes = count_table_genes[,order(ncol(count_table_genes):1)]


### Gene diversity --------
# Maybe it is better to use CSS
#count_table_genes_million_reads = sweep(count_table_genes, 2,ordered_map$Number_of_reads, FUN = "/")
#count_table_genes_million_reads = count_table_genes_million_reads*500000
count_table_genes_zeros = count_table_genes
count_table_genes_zeros[count_table_genes_zeros < 0.8] = 0
metaSeqObject = newMRexperiment(count_table_genes_zeros) #create a metagenomeseq experiment
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) ) #CSS transformation
OTU_read_count_CSS = data.frame(t(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))) # retransofmr into data frame
gene_dissimilairty = vegdist(OTU_read_count_CSS)

permanova = adonis(OTU_read_count_CSS ~ Diagnosis, data = ordered_map)
permanova$aov.tab

species_nmds<-metaMDS(gene_dissimilairty, k=2, try=100)
species_nmds$stress #0.11
ordered_map$NMDS001 = species_nmds$points[,1]
ordered_map$NMDS002 = species_nmds$points[,2]
stress= paste("stress = ", round(species_nmds$stress, digits = 2))
ord1 = ggplot(ordered_map, aes(NMDS001, NMDS002, label = Sample_ID))+
  geom_point(aes(color=Diagnosis), size=4)+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=Diagnosis, color = Diagnosis, label = NA), con.type  = "none")+
  scale_color_manual(values=c('darkred', "darkblue"))+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  xlab("NMDS1") + 
  ylab("NMDS2")+
  theme_bw() +ggtitle("Genes")+
  annotate("text",x=max(ordered_map$NMDS001),y=max(ordered_map$NMDS002),hjust=1,label= stress,  size = 3.5)+
  geom_text(hjust=1.2, vjust=0,  size = 3)
ord1
  
### Richness boxplot
ordered_map$gene_richness  = specnumber(count_table_genes_zeros, MARGIN = 2)
gene_richness_plot = ggplot(ordered_map, aes(y = gene_richness, x = Diagnosis, fill = Diagnosis))+ # for label --> aes(label = metadata$sample)
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = Diagnosis))+ 
  xlab(NULL) + 
  ylab("Number of genes")+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Gene Richness") #+ stat_compare_means(method = "t.test")
gene_richness_plot
summary(lm(gene_richness ~ Diagnosis, data = ordered_map))
t.test(gene_richness ~ Diagnosis, data = ordered_map)


### How many genes have a KO annotation
annotations = read.delim("annotations.txt",  fill =TRUE, row.names = NULL, header = FALSE )
### reads that are annotated
df = data.frame(ortholog = c("annotated", "unannotated"), count = c((nrow(annotations)-sum(annotations$V2 == "")),sum(annotations$V2 == "")))
df$fraction = c(df$count/sum(df$count)*100)
nrow(annotations)
sum(annotations$V2 == "")
# Barplot
bp<- ggplot(df, aes(x="", y=count, fill=ortholog))+
  geom_bar(width = 1, stat = "identity")
bp
library(scales)

pie <- bp + coord_polar("y", start=0) + 
  theme_minimal() +
  geom_text(aes(y = c(750000, 250000), 
                label = paste(round(df$fraction, digits = 1), "%", sep="")))
pie



# Number of genes per sample
genes_without_annotation = annotations$V1[annotations$V2 == ""]
genes_without_annotation = substr(genes_without_annotation,1,nchar(genes_without_annotation)-2)
count_table_annotations = apply(count_table_genes, 2,function(x) ifelse(x > 0, 1, 0))
total_genes = colSums(count_table_annotations)
count_table_no_annotations = count_table_annotations[rownames(count_table_annotations) %in% genes_without_annotation,]
count_table_no_annotations = colSums(count_table_no_annotations)
fraction_no_annotation = count_table_no_annotations/total_genes
ordered_map$fraction_no_annotation = fraction_no_annotation
ggplot(ordered_map, aes(y = fraction_no_annotation, x = Diagnosis, fill = Diagnosis))+ # for label --> aes(label = metadata$sample)
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = Diagnosis))+ 
  xlab(NULL) + 
  ylab("Fraction of genes")+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Fraction of unannotated genes")
summary(lm(fraction_no_annotation ~ Number_of_reads +fraction_no_annotation, data = ordered_map))







#### KO mapping (check if it is coverage or )
### Check for specific genes -----------------
#Ko_mapping1 = sweep(Ko_mapping,2, ordered_map$Number_of_reads, FUN = "/")
#Ko_mapping1 = Ko_mapping1*500000
Ko_mapping1 = Ko_mapping
metaSeqObject = newMRexperiment(Ko_mapping1) #create a metagenomeseq experiment
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) ) #CSS transformation
OTU_read_count_CSS = data.frame(t(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))) # retranso

functional_dissimilairty = vegdist(OTU_read_count_CSS)
set.seed(1)
perm = adonis2(functional_dissimilairty~Diagnosis, ordered_map)
perm
species_nmds<-metaMDS(functional_dissimilairty, k=2, try=100)
species_nmds$stress #0.11
ordered_map$NMDS001 = species_nmds$points[,1]
ordered_map$NMDS002 = species_nmds$points[,2]
stress= paste("stress = ", round(species_nmds$stress, digits = 2))
ord2 = ggplot(ordered_map, aes(NMDS001, NMDS002, label = Sample_ID))+
  geom_point(aes(color=Diagnosis), size=4)+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=Diagnosis, color = Diagnosis, label = NA), con.type  = "none")+
  scale_color_manual(values=c('darkred', "darkblue"))+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  xlab("NMDS1") + 
  ylab("NMDS2")+
  theme_bw() +ggtitle("KO orthologs")+
  annotate("text",x=max(ordered_map$NMDS001),y=max(ordered_map$NMDS002),hjust=1,label= stress,  size = 3.5)+
  geom_text(hjust=1.2, vjust=0,  size = 3)
ord2

ordered_map$functional_richness = specnumber(t(Ko_mapping))
ko_richness = ggplot(ordered_map, aes(y = functional_richness, x = Diagnosis, fill = Diagnosis))+ # for label --> aes(label = metadata$sample)
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = Diagnosis))+ 
  xlab(NULL) + 
  ylab("Number of orthologs")+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("KO Richness")#+ stat_pvalue_manual()

ko_richness
summary(lm(functional_richness ~ Number_of_reads +Diagnosis, data = ordered_map))
t.test(functional_richness ~ Diagnosis, data = ordered_map)

ggarrange(gene_richness_plot, ko_richness)
ggarrange(ord1, ord2, common.legend = TRUE)






###different KOs categories
KO_table_description = read.delim("KO_Orthology_Feb_2022.txt", header = TRUE)
KO_table_description_metabolism = KO_table_description[KO_table_description$KEGG_1 == "Metabolism" |
                                                       KO_table_description$KEGG_2 == "Protein families: metabolism" |
                                                       KO_table_description$KEGG_2 == "Unclassified: metabolism",]


Ko_mapping_metabolism = Ko_mapping1[rownames(Ko_mapping1) %in% KO_table_description_metabolism$KO,]
metaSeqObject = newMRexperiment(Ko_mapping_metabolism) #create a metagenomeseq experiment
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) ) #CSS transformation
OTU_read_count_CSS = data.frame(t(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE)))
functional_dissimilairty = vegdist(OTU_read_count_CSS)
perm = adonis2(functional_dissimilairty~Diagnosis, ordered_map)
perm
species_nmds<-metaMDS(functional_dissimilairty, k=2, try=100)
species_nmds$stress #0.11
ordered_map$NMDS001 = species_nmds$points[,1]
ordered_map$NMDS002 = species_nmds$points[,2]
stress= paste("stress = ", round(species_nmds$stress, digits = 2))
ord1 = ggplot(ordered_map, aes(NMDS001, NMDS002, label = Sample_ID))+
  geom_point(aes(color=Diagnosis), size=4)+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=Diagnosis, color = Diagnosis, label = NA), con.type  = "none")+
  scale_color_manual(values=c('darkred', "darkblue"))+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  xlab("NMDS1") + 
  ylab("NMDS2")+
  theme_bw() +ggtitle("Metabolism related KOs")+
  annotate("text",x=max(ordered_map$NMDS001),y=max(ordered_map$NMDS002),hjust=1,label= stress,  size = 3.5)+
  geom_text(hjust=1.2, vjust=0,  size = 3)
ord1
ordered_map$functional_richness = specnumber(t(Ko_mapping_metabolism))
ggplot(ordered_map, aes(y = functional_richness, x = Diagnosis, fill = Diagnosis))+ # for label --> aes(label = metadata$sample)
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = Diagnosis))+ 
  xlab(NULL) + 
  ylab("Number of orthologs")+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("KO richness - metabolism")
summary(lm(functional_richness ~ Number_of_reads +Diagnosis, data = ordered_map))




#### KO _ 
KO_table_description_genes = KO_table_description[KO_table_description$KEGG_1 == "Genetic Information Processing" |
                                                         KO_table_description$KEGG_2 == "Protein families: genetic information processing" |
                                                         KO_table_description$KEGG_2 == "Unclassified: genetic information processing",]


Ko_mapping_genes = Ko_mapping1[rownames(Ko_mapping1) %in% KO_table_description_genes$KO,]
metaSeqObject = newMRexperiment(Ko_mapping_genes) #create a metagenomeseq experiment
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) ) #CSS transformation
OTU_read_count_CSS = data.frame(t(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE)))
functional_dissimilairty = vegdist(OTU_read_count_CSS)
perm = adonis2(functional_dissimilairty~Diagnosis, ordered_map)
perm
species_nmds<-metaMDS(functional_dissimilairty, k=2, try=100)
species_nmds$stress #0.11
ordered_map$NMDS001 = species_nmds$points[,1]
ordered_map$NMDS002 = species_nmds$points[,2]
stress= paste("stress = ", round(species_nmds$stress, digits = 2))
ord1 = ggplot(ordered_map, aes(NMDS001, NMDS002, label = Sample_ID))+
  geom_point(aes(color=Diagnosis), size=4)+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=Diagnosis, color = Diagnosis, label = NA), con.type  = "none")+
  scale_color_manual(values=c('darkred', "darkblue"))+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  xlab("NMDS1") + 
  ylab("NMDS2")+
  theme_bw() +ggtitle("Genetic information processing KOs")+
  annotate("text",x=max(ordered_map$NMDS001),y=max(ordered_map$NMDS002),hjust=1,label= stress,  size = 3.5)+
  geom_text(hjust=1.2, vjust=0,  size = 3)
ord1
ordered_map$functional_richness = specnumber(t(Ko_mapping_genes))
ggplot(ordered_map, aes(y = functional_richness, x = Diagnosis, fill = Diagnosis))+ # for label --> aes(label = metadata$sample)
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = Diagnosis))+ 
  xlab(NULL) + 
  ylab("Number of orthologs")+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("KO richness - genetic information processing")
summary(lm(functional_richness ~ Number_of_reads +Diagnosis, data = ordered_map))


#### KO_cellular_processes ---------- 
KO_table_description_cell_processes = KO_table_description[KO_table_description$KEGG_1 == "Cellular Processes" |
                                                    KO_table_description$KEGG_2 == "Protein families: genetic information processing" |
                                                    KO_table_description$KEGG_2 == "Unclassified: genetic information processing",]


Ko_mapping_cell = Ko_mapping1[rownames(Ko_mapping1) %in% KO_table_description_cell_processes$KO,]
metaSeqObject = newMRexperiment(Ko_mapping_cell) #create a metagenomeseq experiment
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) ) #CSS transformation
OTU_read_count_CSS = data.frame(t(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE)))
functional_dissimilairty = vegdist(OTU_read_count_CSS)
perm = adonis(functional_dissimilairty~Diagnosis, ordered_map)
perm$aov.tab
species_nmds<-metaMDS(functional_dissimilairty, k=2, try=100)
species_nmds$stress #0.11
ordered_map$NMDS001 = species_nmds$points[,1]
ordered_map$NMDS002 = species_nmds$points[,2]
stress= paste("stress = ", round(species_nmds$stress, digits = 2))
ord1 = ggplot(ordered_map, aes(NMDS001, NMDS002, label = Sample_ID))+
  geom_point(aes(color=Diagnosis), size=4)+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=Diagnosis, color = Diagnosis, label = NA), con.type  = "none")+
  scale_color_manual(values=c('darkred', "darkblue"))+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  xlab("NMDS1") + 
  ylab("NMDS2")+
  theme_bw() +ggtitle("Cellular processes KOs")+
  annotate("text",x=max(ordered_map$NMDS001),y=max(ordered_map$NMDS002),hjust=1,label= stress,  size = 3.5)+
  geom_text(hjust=1.2, vjust=0,  size = 3)
ord1
ordered_map$functional_richness = specnumber(t(Ko_mapping_cell))
ggplot(ordered_map, aes(y = functional_richness, x = Diagnosis, fill = Diagnosis))+ # for label --> aes(label = metadata$sample)
  #geom_label_repel()+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(aes(color = Diagnosis))+ 
  xlab(NULL) + 
  ylab("Number of orthologs")+
  scale_fill_manual(values=c('darkred', "darkblue"))+
  scale_color_manual(values=c('darkred', "darkblue"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("KO richness - Cellular processes")
summary(lm(functional_richness ~ Number_of_reads +Diagnosis, data = ordered_map))




### differential abundance analysis of KOs
Ko_mapping2 = Ko_mapping1[rowSums(Ko_mapping1 != 0) > 4,]
resilts = wilcox_calc(t(Ko_mapping2), ordered_map$Diagnosis)
sig_results = resilts[resilts$p_value<0.05,]
sig_results_no_inf = sig_results[!sig_results$l2cf == Inf,]
max(sig_results_no_inf$l2cf)
sig_results_no_inf = sig_results_no_inf[!sig_results_no_inf$l2cf == -Inf,]
min(sig_results_no_inf$l2cf) # - 10
resilts$l2cf = ifelse(resilts$l2cf == Inf, 10, resilts$l2cf)
resilts$l2cf = ifelse(resilts$l2cf == -Inf, -10, resilts$l2cf)
resilts$descrition =  KO_table_description$Description[match(rownames(resilts),KO_table_description$KO )]
resilts$Kegg3  = KO_table_description$KEGG_3[match(rownames(resilts),KO_table_description$KO )]
sig_results = resilts[resilts$p_value<0.05,]



### Vulcano plot






 
### Calculate average relative abundance 
row_means = rowMeans(replace(Ko_mapping2, Ko_mapping2 == 0, NA), na.rm = TRUE)
sig_results$rowMeans = row_means[match(rownames(sig_results),names(row_means))]




rownames(ordered_map) = ordered_map$Sample_ID
#ordered_map$Diagnosis[ordered_map$Diagnosis == "CTRL"] = "a_CTRL" 
fit_data = Maaslin2(
  input_data = Ko_mapping1, 
  input_metadata = ordered_map, 
  output = "demo_output", 
  normalization = "CSS",
  min_prevalence	= 0.3,
  fixed_effects = "Diagnosis",
  reference = "Diagnosis,CTRL")


### results
results = fit_data$results

results$descrition =  KO_table_description$Description[match(results$feature,KO_table_description$KO )]
results$Kegg3  = KO_table_description$KEGG_3[match(results$feature,KO_table_description$KO )]
results$kegg2  = KO_table_description$KEGG_2[match(results$feature,KO_table_description$KO )]
results$kegg1  = KO_table_description$KEGG_1[match(results$feature,KO_table_description$KO )]

row_means = rowMeans(replace(Ko_mapping2, Ko_mapping2 == 0, NA), na.rm = TRUE)
results$rowMeans = row_means[match(results$feature,names(row_means))]

results_sig_masalin = results[results$pval< 0.05,]
intersect(results_sig_masalin$feature,rownames(sig_results))
both_sig = intersect(results_sig_masalin$feature,rownames(sig_results))
results_sig_bioth = results_sig_masalin[results_sig_masalin$feature %in%both_sig, ] # Gonna use this
results_sig_bioth_increasing = results_sig_bioth[results_sig_bioth$coef > 1,]
results_sig_bioth_dec = results_sig_bioth[results_sig_bioth$coef < -1,]
results$differential = ifelse(results$feature %in%results_sig_bioth_increasing$feature, "CTRL", "not sig")
results$differential = ifelse(results$feature %in%results_sig_bioth_dec$feature, "CPAE", results$differential)
results$label1 = ifelse(results$differential == "not sig", NA , results$feature)


#### Make a vulcano plot
ggplot(data=results, aes(x=coef, y=-log10(pval), col=differential, label=label1)) + 
  geom_point() + 
  theme_minimal() +
  geom_text() + xlab("linear coefficient")+
  scale_fill_manual(values=c('darkred', "darkblue", "grey"))+
  scale_color_manual(values=c('darkred', "darkblue","grey"))

xlsx::write.xlsx(results, "/Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/supplementary_information/Supplementary_table2.xlsx")

