#### plot metabolites
library(ggplot2)

setwd("/Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/metabolites")

#load the file

plot = xlsx::read.xlsx("plot_shabanna.xlsx",1)

#### Aminoacids
plot_aminoacids = plot[plot$Super.Pathway == "Amino Acid",]

a = ggplot(plot_aminoacids) +
  geom_point(aes(x = CPAE.CONTROL, y = Biochemical.Name)) +
  theme_bw() + ylab(NULL) +ggtitle("Amino acids") + facet_grid(Description~., scales= "free")+
  theme(plot.title = element_text(size=9)) + geom_vline(xintercept = 0) + xlab(NULL)
a

### Lipids 
plot_lipids = plot[plot$Super.Pathway == "Lipid",]

b = ggplot(plot_lipids) +
  geom_point(aes(x = CPAE.CONTROL, y = Biochemical.Name)) +
  theme_bw() + ylab(NULL) +ggtitle("Lipids") + facet_grid(Description~., scales= "free")+
  theme(plot.title = element_text(size=9)) + geom_vline(xintercept = 0) + xlab("l2fc")
b





### nucleotide 
plot_Nucleotide = plot[plot$Super.Pathway%in% c("Nucleotide", "Cofactors and Vitamins"),]

c = ggplot(plot_Nucleotide) +
  geom_point(aes(x = CPAE.CONTROL, y = Biochemical.Name)) +
  theme_bw() + ylab(NULL) +ggtitle("Others") + facet_grid(Description~., scales= "free")+
  theme(plot.title = element_text(size=9)) + geom_vline(xintercept = 0) + xlab("l2fc")
c

library(ggpubr)
ggarrange("a","b","c")
library(patchwork)

# Arrange plots horizontally
one = a / c
one | b


