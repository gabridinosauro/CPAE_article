### Gut metabolites 
library(omixerRpm)
setwd("/Users/gabri/OneDrive - University of Arizona/Meta_gen_Daniel/results")
map = read.csv("SampleSheetUsed.csv")
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Meta_gen_Daniel/MAGs")
Ko_mapping = read.csv("KO_table_coverage_with_zeros.csv", row.names = 1)
Ko_mapping = Ko_mapping[,order(ncol(Ko_mapping):1)]
#Ko_mapping = Ko_mapping[-1,]
Ko_mapping = cbind(rownames(Ko_mapping), Ko_mapping)
colnames(Ko_mapping)[1] = "entry"
mods = rpm(Ko_mapping, normalize.by.length = TRUE)



coverage <- asDataFrame(mods, "coverage")
rownames(coverage) = coverage[,2]
coverage = coverage[,-c(1,2)]

abundance = asDataFrame(mods, "abundance")
rownames(abundance) = abundance[,2]
abundance = abundance[,-c(1,2)]
rownames(map) =paste("X", map$Sample_ID, sep = "")

for(i in 1:ncol(abundance)) { abundance[,i] = ifelse(coverage[,i] < 0.60, 0, abundance[,i]) }
abundance = abundance[,match(rownames(map),colnames(abundance))]
colnames(abundance) == rownames(map)

fit_data = Maaslin2(
  input_data = abundance, 
  input_metadata = map, 
  output = "demo_output", 
  #normalization = "none",
  min_prevalence	= 0.3,
  fixed_effects = c("Diagnosis"),
  reference = "Diagnosis,CTRL",plot_heatmap = FALSE,
  plot_scatter = FALSE,
  )


### results
results = fit_data$results
df_plots = data.frame(module = as.numeric(abundance[rownames(abundance) == "Glutamate degradation I",]) , control = map$Diagnosis)
source("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")
krusk = wilcox_calc(t(abundance), map$Diagnosis)




#### Gut brain axis
db <- loadDB("GBMs.v1.0")
mods = rpm(Ko_mapping,  module.db =  loadDB("GBMs.v1.0"))
mods@coverage
coverage <- asDataFrame(mods, "coverage")


rownames(coverage) = coverage[,2]
coverage = coverage[,-c(1,2)]

abundance1 = asDataFrame(mods, "abundance")
rownames(abundance1) = abundance1[,2]
abundance1 = abundance1[,-c(1,2)]
abundance1 = abundance1[,match(rownames(map),colnames(abundance1))]
colnames(abundance1) == rownames(map)

rownames(map) =paste("X", map$Sample_ID, sep = "")
fit_data = Maaslin2(
  input_data = abundance1, 
  input_metadata = map, 
  output = "demo_output", 
  #normalization = "none",
  min_prevalence	= 0.3,
  fixed_effects = c("Diagnosis"),
  reference = "Diagnosis,CTRL",plot_heatmap = FALSE,
  plot_scatter = FALSE,
)
results1 = fit_data$results

### Plot results boxplots
### Boxplots 
results_sig = results[which(results$pval<0.1),]
results_sig$coef = results_sig$coef *-1
results_sig = results_sig[]
results_sig$feature = gsub("\\.", " ",results_sig$feature)
results_sig$significance = ifelse(results_sig$pval<0.05, "<0.05", "<0.1")
results_sig$module = "GMM"
results_gut_functions = abundance[which(rownames(abundance) %in% results_sig$feature),]

results_sig1 = results1[which(results1$pval<0.1),]
results_sig1$coef = results_sig1$coef *-1
results_sig1 = results_sig1[]
results_sig1$feature = gsub("\\.", " ",results_sig1$feature)
results_sig1$significance = ifelse(results_sig1$pval<0.05, "<0.05", "<0.1")
results_sig1$module = "GBM"

results_sig = results_sig[results_sig$feature != "glutamate degradation I",]
results_sig$module = ifelse(results_sig$feature == "glutamate degradation II", "GBM", "GMM")

results_tot = rbind(results_sig, results_sig1)

#results_gut_functions = abundance[which(rownames(abundance) %in% results_tot$feature),]




a = ggplot(results_tot) +
  geom_point(aes(x = coef, y = feature, fill = significance), pch = 21) +
  theme_bw() + ylab(NULL) +ggtitle("CPAE vs CONTROL") + facet_grid(module~., scales = "free", space = "free") +
  theme(plot.title = element_text(size=9)) + geom_vline(xintercept = 0) +
  scale_y_discrete(limits=rev) + scale_fill_manual(values = c("darkred","orange")) + theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#+
#labs(col = "Phylum")
a

#geom_vline(xintercept = 0)
rel_ab_stool$status = stool.md$status
melted_ASV = reshape2::melt(rel_ab_stool)
melted_ASV$variable = factor(melted_ASV$Var2,rev(as.character(rownames(results_patch)) ))
melted_ASV = melted_ASV[melted_ASV$variable %in% res$feature,]
melted_ASV$variable  = factor(melted_ASV$variable, levels = res$feature)

b = ggplot(melted_ASV, aes(value, variable, color = status))+
  geom_boxplot(outlier.shape = NA) +theme_bw() +
  #geom_point(position = position_jitterdodge(), alpha=0.3)+
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab(NULL) + xlab("relative abundance") +
  scale_color_manual(values=c('darkred', "darkblue"))+
  scale_y_discrete(limits=rev) +
  xlim(0, 0.1) + theme(legend.position = "none")
b
ggarrange(a,b, align = "h",widths = c(0.75,0.3))











### Get the lines for each file
file = readLines("GMMs.v1.07.txt") 
file1 = file[grep("MF",file)]
# Specify the file path
file_path <- "path/to/your/file.txt"

# Read the file
lines <- readLines("GMMs.v1.07.txt")


# Find the indices of the starting marker
start_indices <- grep(start_marker, lines)

# Find the indices of the ending marker
end_indices <- grep(end_marker, lines)

# Extract the lines between the markers
lines_between_markers <- character()

# Iterate over each start marker index
for (start_index in start_indices) {
  # Find the first end marker index that comes after the start marker
  end_index <- min(end_indices[end_indices > start_index])
  
  # Extract the lines between the markers for this section
  section_lines <- lines[(start_index + 1):(end_index - 1)]
  
  # Append the section lines to the lines_between_markers variable
  lines_between_markers <- c(lines_between_markers, section_lines)
}

# Print the lines between the markers
print(lines_between_markers)







