library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(viridis)
library(microbiome)
library(reshape2)
library(ggpubr)
library(broom)
library(ggfortify)
library(dplyr)
library(phyloseq)
library(decontam)
library(vegan)
library(tidyverse)
library(indicspecies)
library(ggalluvial)
library(pheatmap)
library(metagenomeSeq)
library(DESeq2)
theme_set(theme_light())

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("processx")
library(devtools)

devtools::install_github("vmikk/metagMisc")
library(metagMisc)

install.packages("devtools") #Installs devtools (if not already installed)
devtools::install_github("donaldtmcknight/microDecon") #Installs microDecon

#Set working directory
setwd("/Volumes/ShadeLab/WorkingSpace/MarcoMechan_WorkingSpace/sorghum_16S_its_R_analysis/")

#Import files
otu=read.table("otu_table.txt", header = T, row.names = 1)
tax=read.delim("taxonomy.tsv",row.names = 1)
map=read.table("metadata_16S_20-21.csv", header = T, row.names = 1, sep=",")
tree=read_tree("merged_16S_rooted_tree.nwk")

muc.its.otu=read.table("its_otu-table.tsv", header = T, row.names = 1)
muc.its.tax=read.delim("its_taxonomy.tsv",row.names = 1)
muc.its.map=read.table("metadata_ITS_20-21.csv", header = T, row.names = 1, sep=",")
muc.its.tree=read_tree("its_rooted_tree.nwk")

################################################################################
####################   Bacteria community analysis  ############################
################################################################################

tax=tax[-2]
tax_df <- colsplit(tax$Taxon, '; ', names =  c("Kingdom", "Phylum", "Class", 
                                               "Order", "Family", "Genus", "Species"))
tax_df[1:7] <- lapply(tax_df[1:7], function(x) gsub(".*__", "", x))
rownames(tax_df) <- rownames(tax)

OTU=otu_table(as.matrix(otu), taxa_are_rows = T, errorIfNULL = TRUE)
TAX=tax_table(as.matrix(tax_df), errorIfNULL = TRUE)
MAP=sample_data(data.frame(map))

otuPhyloseq=phyloseq(OTU,TAX,MAP)
sample_sums(otuPhyloseq)

summarize_phyloseq(otuPhyloseq)

Muc.outphyloseq <- subset_samples(otuPhyloseq, Compartment%in%c("Mucilage")&Year%in%c("2020"))
Muc.outphyloseq <- prune_taxa(taxa_sums(Muc.outphyloseq) > 0, Muc.outphyloseq)
summarize_phyloseq(Muc.outphyloseq)
sample_sums(Muc.outphyloseq)


Wax.outphyloseq <- subset_samples(otuPhyloseq, Compartment%in%c("Epicuticular.Wax")&Year%in%c("2020"))
Wax.outphyloseq <- prune_taxa(taxa_sums(Wax.outphyloseq) > 0, Wax.outphyloseq)
summarize_phyloseq(Wax.outphyloseq)
sample_sums(Wax.outphyloseq)


#Filtering mitochondria, Chloroplast and Unclassified taxa
otuPhyloseq_filt <- otuPhyloseq %>%
  subset_taxa(Kingdom != "Archaea")

otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa(Kingdom != "Eukaryota")

otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa(Kingdom != "Unassigned")

otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa((Family != "Mitochondria") | is.na(Family))

otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa((Order != "Chloroplast") | is.na(Class))

sample_sums(otuPhyloseq_filt) 
print(otuPhyloseq_filt) 

summarize_phyloseq(otuPhyloseq_filt)

filtered_otus <- phyloseq_to_df(otuPhyloseq_filt, addtax=T, addtot=T)


Muc.outphyloseq.filt <- subset_samples(otuPhyloseq_filt, Compartment%in%c("Mucilage")&Year%in%c("2021"))
Muc.outphyloseq.filt <- prune_taxa(taxa_sums(Muc.outphyloseq.filt) > 0, Muc.outphyloseq.filt)
summarize_phyloseq(Muc.outphyloseq.filt)
sample_sums(Muc.outphyloseq.filt)

Wax.outphyloseq.filt <- subset_samples(otuPhyloseq_filt, Compartment%in%c("Epicuticular.Wax"))
Wax.outphyloseq.filt <- prune_taxa(taxa_sums(Wax.outphyloseq.filt) > 0, Wax.outphyloseq.filt)
summarize_phyloseq(Wax.outphyloseq.filt)
sample_sums(Wax.outphyloseq.filt) 


############################ Subset Mucilage samples #################################

#### Sample G2R1.1.M8 was removed due inconsistency compared with all other Mucilage samples 
Muc <- subset_samples(otuPhyloseq_filt, Compartment%in%c("Mucilage"))
Muc <- subset_samples(Muc, sample_names(Muc) != 'G2R1.1.M8')
Muc <- prune_taxa(taxa_sums(Muc) > 0, Muc)

### remove singletons
Muc <- prune_taxa(taxa_sums(Muc) > 1, Muc)
print(Muc) ####12047 taxa and 179 samples

#### Samples with less than 20K reads were removed: 21 samples in total were removed
Muc.1 <- prune_samples(sample_sums(Muc) >= 20000, Muc)
Muc.1 <- prune_taxa(taxa_sums(Muc.1) > 0, Muc.1)
print(Muc.1) ######11169 taxa and 158 samples 

### Samples were rarefied to 20519 reads/sample
set.seed(13)
Muc.rare = rarefy_even_depth(Muc.1, rngseed=1, sample.size=min(sample_sums(Muc.1)), replace=F)
print (Muc.rare) ### 9492 taxa and 158 samples

### Mucilage Alpha Diversity
### Rarefaction curves
rarecurve.out = rarecurve(t(otu_table(Muc)), step=1000, label=FALSE)
names(rarecurve.out) <- sample_data(Muc)$sample.name

protox.muc <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "rare_count"
  mydf$sample.name <- y
  mydf$Subsample_size <- attr(x, "Subsample")
  mydf
}, x = rarecurve.out, y = as.list(names(rarecurve.out)), SIMPLIFY = FALSE)

rarecurve.df = do.call(rbind, protox.muc) %>%
  left_join(data.frame(sample_data(Muc)), by = "sample.name")

Muc.rare.curves = ggplot(data=rarecurve.df, aes(x=Subsample_size, y=rare_count)) +
  geom_vline(xintercept = 20519, color= "red",  linetype='dashed') +
  scale_color_manual(values=c("#56B4E9")) +
  geom_line(aes(group=sample.name, color=Compartment)) +
  labs(x="Rarefied read depth", y="Bacterial ASVs count") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8,angle=0, hjust=0.5),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10,face="bold", vjust = 10),
        legend.title = element_text(size=0),
        legend.position = c(.75, .1),
        legend.text = element_text(colour="#56B4E9", size = 12, face = "bold"))+
  ggtitle("B) Bacterial alpha diversity")+
  theme(plot.title = element_text(size = 10, face = "bold"))
Muc.rare.curves


### Diversity metrics
Muc.alpha.df = rbind(data.frame(alpha_measure = specnumber(t(otu_table(Muc.rare)))) %>%
                       mutate(Index = "Observed taxa", Compartment="Mucilage"),
                     data.frame(alpha_measure = picante::pd(t(otu_table(Muc.rare)), tree)$PD) %>%
                       mutate(Index = "Phylogenetic", Compartment="Mucilage"))

alpha_div_muc <- ggplot(Muc.alpha.df, aes(x=Compartment, y=alpha_measure, color=Compartment)) +
  geom_boxplot(alpha=0.6) +
  scale_color_manual(values=c("#56B4E9")) +
  geom_jitter(width=0.15, alpha=0.9) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme_bw() +
  theme(axis.text.x=element_text(size=0,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size=10,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=10, face = 'bold'),
        plot.title = element_text(size = rel(2)), 
        axis.title=element_text(size=12,face="bold", vjust = 10),
        legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~Index, scales = "free_y")
alpha_div_muc

Muc.rich = estimate_richness(Muc.rare)
Muc.DAE.obs = wilcox.test(Muc.rich$Observed ~ sample_data(Muc.rare)$Developmental.Stage, data = Muc.rich, exact = FALSE, conf.int = TRUE)
Muc.DAE.Chao1 = wilcox.test(Muc.rich$Chao1 ~ sample_data(Muc.rare)$Developmental.Stage, data = Muc.rich, exact = FALSE, conf.int = TRUE)
Muc.DAE.Simp = wilcox.test(Muc.rich$Simpson ~ sample_data(Muc.rare)$Developmental.Stage, data = Muc.rich, exact = FALSE, conf.int = TRUE)
Muc.DAE.Shan = wilcox.test(Muc.rich$Shannon ~ sample_data(Muc.rare)$Developmental.Stage, data = Muc.rich, exact = FALSE, conf.int = TRUE)

Muc.DAE.alpha.wicox = data.frame(index = c("Observed", "Chao1", "Shannon", "Simpson"),
                             pvalue = c(Muc.DAE.obs$p.value,Muc.DAE.Chao1$p.value,
                                        Muc.DAE.Shan$p.value, Muc.DAE.Simp$p.value))
Muc.DAE.alpha.wicox

Muc.Fert.obs = wilcox.test(Muc.rich$Observed ~ sample_data(Muc.rare)$Fertilization.Status, data = Muc.rich, exact = FALSE, conf.int = TRUE)
Muc.Fert.Chao1 = wilcox.test(Muc.rich$Chao1 ~ sample_data(Muc.rare)$Fertilization.Status, data = Muc.rich, exact = FALSE, conf.int = TRUE)
Muc.Fert.Simp = wilcox.test(Muc.rich$Simpson ~ sample_data(Muc.rare)$Fertilization.Status, data = Muc.rich, exact = FALSE, conf.int = TRUE)
Muc.Fert.Shan = wilcox.test(Muc.rich$Shannon ~ sample_data(Muc.rare)$Fertilization.Status, data = Muc.rich, exact = FALSE, conf.int = TRUE)

Muc.fert.alpha.wicox = data.frame(index = c("Observed", "Chao1", "Shannon", "Simpson"),
                             pvalue = c(Muc.Fert.obs$p.value,Muc.Fert.Chao1$p.value,
                                        Muc.Fert.Shan$p.value, Muc.Fert.Simp$p.value))
Muc.fert.alpha.wicox

Muc.Year.obs = wilcox.test(Muc.rich$Observed ~ sample_data(Muc.rare)$Year, data = Muc.rich, exact = FALSE, conf.int = TRUE)
Muc.Year.Chao1 = wilcox.test(Muc.rich$Chao1 ~ sample_data(Muc.rare)$Year, data = Muc.rich, exact = FALSE, conf.int = TRUE)
Muc.Year.Simp = wilcox.test(Muc.rich$Simpson ~ sample_data(Muc.rare)$Year, data = Muc.rich, exact = FALSE, conf.int = TRUE)
Muc.Year.Shan = wilcox.test(Muc.rich$Shannon ~ sample_data(Muc.rare)$Year, data = Muc.rich, exact = FALSE, conf.int = TRUE)

Muc.Year.alpha.wicox = data.frame(index = c("Observed", "Chao1", "Shannon", "Simpson"),
                                  pvalue = c(Muc.Year.obs$p.value,Muc.Year.Chao1$p.value,
                                             Muc.Year.Shan$p.value, Muc.Year.Simp$p.value))
Muc.Year.alpha.wicox


alpha_meas = c("Observed", "Shannon", "Simpson")
(p_muc.bact <- plot_richness(Muc.rare, x= "Year", measures=alpha_meas, color="Year"))
p_muc.bact + geom_boxplot(alpha=0.6) +
  geom_jitter(aes(color=Year), width=0.15, alpha=0.9)+
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme(axis.text.x=element_text(size=14,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10))

### Merge Rarefaction curves and Div metrics
Figure1b <- cowplot::plot_grid(Muc.rare.curves, alpha_div_muc, ncol=1, rel_heights = c(1, 0.9), align = "v",
                   rel_widths = c(1,2))
Figure1b


### Mucilage Beta Diversity based on Bray-Curtis dissimilarity

bray_diss_Muc = phyloseq::distance(Muc.rare, method="bray")
ordination_Muc = ordinate(Muc.rare, method="PCoA", distance=bray_diss_Muc)
# ADONIS test
Muc.bact.FT.adonOut<- adonis2(bray_diss_Muc ~ sample_data(Muc.rare)$Fertilization.Status)
pVal.bact.FT=Muc.bact.FT.adonOut$`Pr(>F)`[1]
r2.bact.FT = round(Muc.bact.FT.adonOut$R2[1], digits = 4)

Muc.bact.DAE.adonOut<- adonis2(bray_diss_Muc ~ sample_data(Muc.rare)$Developmental.Stage)
pVal.bact.DAE=Muc.bact.DAE.adonOut$`Pr(>F)`[1]
r2.bact.DAE = round(Muc.bact.DAE.adonOut$R2[1], digits = 4)

Muc.Year.adonOut<- adonis2(bray_diss_Muc ~ sample_data(Muc.rare)$Year)
pVal.bact.Y=Muc.Year.adonOut$`Pr(>F)`[1]
r2.bact.Y = round(Muc.Year.adonOut$R2[1], digits = 4)


Muc.DAE.FS.adonOut<- adonis2(bray_diss_Muc ~ sample_data(Muc.rare)$Fertilization.Status*sample_data(Muc.rare)$Developmental.Stage)
Muc.DAE.Year.adonOut<- adonis2(bray_diss_Muc ~ sample_data(Muc.rare)$Developmental.Stage*sample_data(Muc.rare)$Year)
Muc.Year.FS.adonOut<- adonis2(bray_diss_Muc ~ sample_data(Muc.rare)$Year*sample_data(Muc.rare)$Fertilization.Status)

Muc.PCoA = plot_ordination(Muc.rare, ordination_Muc, color = 'Developmental.Stage', shape = "Fertilization.Status") + theme(aspect.ratio=1)+
  geom_point(na.rm = FALSE, size=4) + 
  scale_shape_manual(values = c(19, 1))+
  scale_color_manual(values=c("#56B4E9", "#336c8b")) +
  labs(x="PCoA1 (22.7% var. explained)", y="PCoA1 (8.4% var. explained)") +
  theme_bw() +
  theme(axis.text.x=element_text(size=12,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size=12,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=12, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=12,face="bold", vjust = 10),
        legend.text = element_text(size=10),
        legend.title = element_text(size =12, face="bold"))+
  labs(colour="Developmental Stage", shape="Fertilization Status")+
  annotate(geom = 'text', label=paste("PERMANOVA\nR2=", r2.bact, "\nP-val=",pVal.bact, sep=''), x=.05, y=-.3)+
  ggtitle("Bacterial Root Mucilage")

Muc.PCoA$layers <- Muc.PCoA$layers[-1]
Figure3b <- Muc.PCoA
Figure3b

### Betadisper 
Muc.ds.dispr <- betadisper(bray_diss_Muc, phyloseq::sample_data(Muc.rare)$Developmental.Stage)
Muc.ds.dispr
## Permutation test for F
per.Muc.ds.dispr <-permutest(Muc.ds.dispr, pairwise = TRUE, permutations = 999)
per.Muc.ds.dispr
## Tukey's Honest Significant Differences
(Muc.ds.dispr.HSD <- TukeyHSD(Muc.ds.dispr))
plot(Muc.ds.dispr.HSD)
## Plot the groups and distances to centroids on the first two PCoA axes
plot(Muc.ds.dispr)
## Draw a boxplot of the distances to centroid for each group
boxplot(Muc.ds.dispr)

#### Betadisper 
Muc.fs.dispr <- betadisper(bray_diss_Muc, phyloseq::sample_data(Muc.rare)$Fertilization.Status)
Muc.fs.dispr
## Permutation test for F
per.Muc.fs.dispr <-permutest(Muc.fs.dispr, pairwise = TRUE, permutations = 999)
per.Muc.fs.dispr
## Tukey's Honest Significant Differences
(Muc.fs.dispr.HSD <- TukeyHSD(Muc.fs.dispr))
plot(Muc.fs.dispr.HSD)
## Plot the groups and distances to centroids on the first two PCoA axes
plot(Muc.fs.dispr)
## Draw a boxplot of the distances to centroid for each group
boxplot(Muc.fs.dispr)

### Mucilage Bacterial Relative abundance 
Muc_phylum <- Muc.rare %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  #filter(Abundance > 0.005) %>%
  arrange(Phylum)
  
#write.csv(Muc_phylum, 'Mucilage_RA_phylum.csv') 

n <- dim(Muc_phylum)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Devl.names = c(`60DAE` = "60 Days after emergence", `90DAE` = "90 Days after emergence")

Muc.RA.Phylum = ggplot(Muc_phylum,aes(x=Sample,y=Abundance,fill=Phylum, order=as.factor(Phylum))) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1, labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Bacterial Phylum")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Phylum) \n") +
  xlab("Sorghum aerial root mucilage") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=10),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Phylum Composition")
Muc.RA.Phylum

Muc_class <- Muc.rare %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Class)
  
#write.csv(Muc_class, 'Mucilage_RA_class.csv')

n <- dim(Muc_class)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Muc.RA.Class = ggplot(Muc_class,aes(x=Sample,y=Abundance,fill=Class)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1, labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Bacterial Class")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Class > 1%) \n") +
  xlab("Sorghum aerial root mucilage") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Class Composition")
Muc.RA.Class

Muc_family <- Muc.rare %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  #filter(Abundance > 0.03) %>%
  arrange(Family)
write.csv(Muc_family, 'Mucilage_RA_03_family.csv')

n <- dim(Muc_family)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Muc.RA.Family <- ggplot(Muc_family,aes(x=Sample,y=Abundance,fill=Family)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1, labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Bacterial Family")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Family > 3%) \n") +
  xlab("Sorghum aerial root mucilage") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Family composition")
Muc.RA.Family

Muc_genus <- Muc.rare %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.03) %>%
  arrange(Genus)
#write.csv(Muc_genus, 'Mucilage_RA_genus.csv')

n <- dim(Muc_genus)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Muc.RA.Genus = ggplot(Muc_genus,aes(x=Sample,y=Abundance,fill=Genus)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1, labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Bacterial Genus")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Genus > 3%) \n") +
  xlab("Sorghum aerial root mucilage") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=14,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=14, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #ggtitle("Mucilage bacterial diversity")
Muc.RA.Genus


### Differential abundance Analysis - DESeq: Plant Development as variable
deseq_Muc1 = phyloseq_to_deseq2(Muc.rare, ~ Developmental.Stage)
sample_data(Muc.rare)$Developmental.Stage
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_Muc1), 1, gm_mean)
deseq_Muc1 = estimateSizeFactors(deseq_Muc1, geoMeans = geoMeans)
deseq_Muc1 = DESeq(deseq_Muc1, fitType="local")

res = results(deseq_Muc1)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Muc.rare)[rownames(sigtab), ], "matrix"))
head(sigtab)
# To write all OTUs that were significant different: positives and negatives
sigtab = sigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus", "Species")]
#write.csv(sigtab, 'all_Differential_abundance_Muc.ds.csv')

# To subset positives values
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(posigtab, 'Differential_abundance_Muc.ds.csv')

theme_set(theme_bw())
#sigtabgen = subset(sigtab, !is.na(Genus))
sigtabgen = sigtab %>%
  mutate(Genus = ifelse(Genus == "", paste("Unclassified", Family), Genus)) %>%
  mutate(Family = ifelse(Family == "uncultured", paste("", Class), Family))%>%
  mutate(Genus = ifelse(Genus == "uncultured", paste("Unclassified" , Family), Genus))
  
# Class order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Class, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Class = factor(as.character(sigtabgen$Class), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

Muc.DESeq = ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Class)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  scale_color_manual(values = c("#A6D854", "#377EB8","#E78AC3","#F781BF","#F0027F","#FDCDAC","#FF7F00","#E41A1C","#BEBADA"))+
  theme(axis.text.x=element_text(size=12,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size=12,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=10, face = 'bold'),
        plot.title = element_text(size = rel(2), face="bold"),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=10)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Differentially abundant ASVs: 60 vs 90 DAE ")
Muc.DESeq

Figure4 <- Muc.DESeq

### Differential abundance Analysis - DESeq: Fertilization as variable

deseq_Muc.fs = phyloseq_to_deseq2(Muc.rare, ~ Fertilization.Status)
sample_data(Muc.rare)$Fertilization.Status
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_Muc.fs), 1, gm_mean)
deseq_Muc.fs = estimateSizeFactors(deseq_Muc.fs, geoMeans = geoMeans)
deseq_Muc.fs = DESeq(deseq_Muc.fs, fitType="local")

res = results(deseq_Muc.fs)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab.fs = res[(res$padj < alpha), ]
sigtab.fs = cbind(as(sigtab.fs, "data.frame"), as(tax_table(Muc.rare)[rownames(sigtab.fs), ], "matrix"))
head(sigtab.fs)
# To write all OTUs that were significant different: positives and negatives
sigtab.fs = sigtab.fs[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

# To subset positives values
posigtab.fs = sigtab.fs[sigtab.fs[, "log2FoldChange"] > 0, ]
posigtab.fs = posigtab.fs[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(posigtab.fs, 'Differential_abundance_Muc_fs.csv')

theme_set(theme_bw())
sigtabgen.fs = subset(sigtab.fs, !is.na(Genus))
# Class order
x = tapply(sigtabgen.fs$log2FoldChange, sigtabgen.fs$Class, function(x) max(x))
x = sort(x, TRUE)
sigtabgen.fs$Class = factor(as.character(sigtabgen.fs$Class), levels=names(x))
# Genus order
x = tapply(sigtabgen.fs$log2FoldChange, sigtabgen.fs$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen.fs$Genus = factor(as.character(sigtabgen.fs$Genus), levels=names(x))

Muc.DESeq.fs = ggplot(sigtabgen.fs, aes(y=Genus, x=log2FoldChange, color=Class)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  #scale_color_manual(values = col_class) + 
  #scale_color_manual(values = c("#A6D854", "#377EB8","#E78AC3","#F781BF","#F0027F","#FDCDAC","#FF7F00","#E41A1C","#BEBADA"))+
  theme(axis.text.x=element_text(size=12,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size=12,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=10, face = 'bold'),
        plot.title = element_text(size = rel(2), face="bold"),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=10)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Differently enriched taxa in the mucilage: Non Fertilizer vs Fertilizer")
Muc.DESeq.fs



######################## Subset Epicuticular Wax samples ########################

Wax <- subset_samples(otuPhyloseq_filt, Compartment%in%c("Epicuticular.Wax"))
Wax <- prune_taxa(taxa_sums(Wax) > 0, Wax)
### remove singletons
Wax <- prune_taxa(taxa_sums(Wax) > 1, Wax)
print(Wax) ####1112 taxa and 48 samples

summarize_phyloseq(Wax)
#### Samples with less than 1K reads were removed: 12 samples in total were removed
Wax.1 <- prune_samples(sample_sums(Wax) >= 1000, Wax)
Wax.1 <- prune_taxa(taxa_sums(Wax.1) > 0, Wax.1)
print(Wax.1) ###### 1094 taxa and 42 samples 

### Samples were rarefied to 1303 reads/sample
set.seed(13)
Wax.rare = rarefy_even_depth(Wax.1, rngseed=1, sample.size=min(sample_sums(Wax.1)), replace=F)
print (Wax.rare) ### 540 taxa and 42 samples

## Wax Alpha diversity
## rarefaction curves
wax.rarecurve.out = rarecurve(t(otu_table(Wax)), step=1000, label=FALSE)
names(wax.rarecurve.out) <- sample_data(Wax)$sample.name

wax.protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "rare_count"
  mydf$sample.name <- y
  mydf$Subsample_size <- attr(x, "Subsample")
  mydf
}, x = wax.rarecurve.out, y = as.list(names(wax.rarecurve.out)), SIMPLIFY = FALSE)

wax.rarecurve.df = do.call(rbind, wax.protox) %>%
  left_join(data.frame(sample_data(Wax)), by = "sample.name")

Wax.rare.curves = ggplot(data=wax.rarecurve.df, aes(x=Subsample_size, y=rare_count)) +
  geom_vline(xintercept = 1303, color= "red",  linetype='dashed') +
  scale_color_manual(values=c("#66A61E")) +
  geom_line(aes(group=sample.name, color=Compartment)) +
  labs(x="Rarefied read depth", y="Bacterial ASVs count") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8,angle=0, hjust=0.5),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10,face="bold", vjust = 10),
        legend.position = c(.75, .1),
        legend.title = element_text(size=0),
        legend.text = element_text(colour="#66A61E", size = 12, face = "bold"))+
  ggtitle("A) Bacterial alpha diversity")+
  theme(plot.title = element_text(size = 10, face = "bold"))
Wax.rare.curves

### metrics
Wax.alpha.df = rbind(data.frame(alpha_measure = specnumber(t(otu_table(Wax.rare)))) %>%
                       mutate(Index = "Observed taxa", Compartment="Epicuticular Wax"),
                     data.frame(alpha_measure = picante::pd(t(otu_table(Wax.rare)), tree)$PD) %>%
                       mutate(Index = "Phylogenetic", Compartment="Epicuticular Wax"))

alpha_div_wax <- ggplot(Wax.alpha.df, aes(x=Compartment, y=alpha_measure, color=Compartment)) +
  geom_boxplot(alpha=0.6) +
  scale_color_manual(values=c("#66A61E")) +
  geom_jitter(width=0.15, alpha=0.9) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme_bw() +
  theme(axis.text.x=element_text(size=0,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size=10,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=10, face = 'bold'),
        plot.title = element_text(size = rel(2)), 
        axis.title=element_text(size=12,face="bold", vjust = 10),
        legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~Index, scales = "free_y")
alpha_div_wax

Wax.rich = estimate_richness(Wax.rare)
Wax.obs = wilcox.test(Wax.rich$Observed ~ sample_data(Wax.rare)$Developmental.Stage, data = Wax.rich, exact = FALSE, conf.int = TRUE)
Wax.Chao1 = wilcox.test(Wax.rich$Chao1 ~ sample_data(Wax.rare)$Developmental.Stage, data = Wax.rich, exact = FALSE, conf.int = TRUE)
Wax.Simp = wilcox.test(Wax.rich$Simpson ~ sample_data(Wax.rare)$Developmental.Stage, data = Wax.rich, exact = FALSE, conf.int = TRUE)
Wax.Shan = wilcox.test(Wax.rich$Shannon ~ sample_data(Wax.rare)$Developmental.Stage, data = Wax.rich, exact = FALSE, conf.int = TRUE)

wax.alpha.wicox = data.frame(index = c("Observed", "Chao1", "Shannon", "Simpson"),
                             pvalue = c(Wax.obs$p.value,Wax.Chao1$p.value,
                                        Wax.Shan$p.value, Wax.Simp$p.value))
wax.alpha.wicox

#### to merge figure 1a
Figure1a <- cowplot::plot_grid(Wax.rare.curves, alpha_div_wax, ncol=1, rel_heights = c(1, 0.9), align = "v",
                               rel_widths = c(1,2))
Figure1a

## Wax beta diversity based on Bray-Curtis dissimilarity

bray_diss_Wax = phyloseq::distance(Wax.rare, method="bray")
ordination_Wax = ordinate(Wax.rare, method="PCoA", distance=bray_diss_Wax)

Wax.adonOut<- adonis2(bray_diss_Wax ~ sample_data(Wax.rare)$Developmental.Stage)
Wax.pVal=Wax.adonOut$`Pr(>F)`[1]
Wax.r2 = round(Wax.adonOut$R2[1], digits = 4)

Wax.PCoA = plot_ordination(Wax.rare, ordination_Wax, color = 'Developmental.Stage') + theme(aspect.ratio=1)+
  geom_point(size=3) + 
  #scale_color_manual(values=c("#66A61E", "#3d6312")) +
  scale_color_manual(labels = c("60DAE", "90DAE"), values = c("#66A61E", "#33530f")) +
  labs(x="PCoA1 (16.5% var. explained)", y="PCoA1 (15.3% var. explained)") +
  theme_bw() +
  theme(axis.text.x=element_text(size=12,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size=12,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=12, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=12,face="bold", vjust = 10),
        legend.text = element_text(size=10),
        legend.title = element_text(size =12, face="bold"))+
  labs(colour="Developmental Stage")+
  annotate(geom = 'text', label=paste("PERMANOVA\nR2=", Wax.r2, "\nP-val=",Wax.pVal, sep=''), x=-.2, y=-.4)+
  ggtitle("Bacterial Epicuticular Wax")
Wax.PCoA
Figure3a <- Wax.PCoA


#### Betadisper 
Wax.ds.dispr <- betadisper(bray_diss_Wax, phyloseq::sample_data(Wax.rare)$Developmental.Stage)
Wax.ds.dispr
## Permutation test for F
per.Wax.ds.dispr <-permutest(Wax.ds.dispr, pairwise = TRUE, permutations = 999)
per.Wax.ds.dispr
## Tukey's Honest Significant Differences
(Wax.ds.dispr.HSD <- TukeyHSD(Wax.ds.dispr))
plot(Wax.ds.dispr.HSD)
## Plot the groups and distances to centroids on the first two PCoA axes
plot(Wax.ds.dispr)
## Draw a boxplot of the distances to centroid for each group
boxplot(Wax.ds.dispr)


### Wax Relative abundance
Wax_phylum <- Wax.rare %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  #filter(Abundance > 0.01) %>%
  arrange(Phylum)
#write.csv(Wax_phylum, 'Wax_RA_phylum.csv') 

n <- dim(Wax_phylum)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Devl.names = c(`60DAE` = "60 Days after emergence", `90DAE` = "90 Days after emergence")

Wa.RA.Phylum = ggplot(Wax_phylum,aes(x=Sample,y=Abundance,fill=Phylum)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1,labeller = as_labeller(Devl.names) )+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Bacterial Phylum")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Phylum) \n") +
  xlab("Sorghum Epicuticular Wax") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Phylum Composition")
Wa.RA.Phylum

Wax_class <- Wax.rare %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Class)
#write.csv(Wax_class, 'Wax_RA_class.csv')

n <- dim(Wax_class)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Wa.RA.Class = ggplot(Wax_class,aes(x=Sample,y=Abundance,fill=Class)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1, labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Bacterial Class")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Class > 1%) \n") +
  xlab("Sorghum Epicuticular Wax") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Class Composition")
Wa.RA.Class

Wax_family <- Wax.rare %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.03) %>%
  arrange(Family)
#write.csv(Wax_family, 'Wax_RA_family.csv')

n <- dim(Wax_family)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Wax.RA.family <-ggplot(Wax_family,aes(x=Sample,y=Abundance,fill=Family)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1, labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) +
  labs(fill="Bacterial Family")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Family > 3%) \n") +
  xlab("Sorghum Epicuticular Wax") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Family Composition")
Wax.RA.family

Wax_genus <- Wax.rare %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.03) %>%
  arrange(Genus)
#write.csv(Wax_genus, 'Wax_RA_genus.csv')

n <- dim(Wax_genus)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Wax.RA.genus = ggplot(Wax_genus,aes(x=Sample,y=Abundance,fill=Genus)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1, labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Bacterial Genus")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Genus > 2%) \n") +
  xlab("Sorghum Epicuticular Wax") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=7),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  #ggtitle("Genus Composition")
Wax.RA.genus

## Differential abundance with DESeq: Collection Date as variable

deseq_wax1 = phyloseq_to_deseq2(Wax.rare, ~ Developmental.Stage)
sample_data(Wax.rare)$Developmental.Stage
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_wax1), 1, gm_mean)
deseq_wax1 = estimateSizeFactors(deseq_wax1, geoMeans = geoMeans)
deseq_wax1 = DESeq(deseq_wax1, fitType="local")

res.wax = results(deseq_wax1)
res.wax = res.wax[order(res.wax$padj, na.last=NA), ]
alpha = 0.01
sigtab.wax = res.wax[(res.wax$padj < alpha), ]
sigtab.wax = cbind(as(sigtab.wax, "data.frame"), as(tax_table(Wax.rare)[rownames(sigtab.wax), ], "matrix"))
head(sigtab.wax)
# To write all OTUs that were significant different: positives and negatives
sigtab.wax = sigtab.wax[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus", "Species")]
#write.csv(sigtab, 'all_Differential_abundance_Muc.ds.csv')

# To subset positives values
posigtab.wax = sigtab.wax[sigtab.wax[, "log2FoldChange"] > 0, ]
posigtab.wax = posigtab.wax[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(posigtab, 'Differential_abundance_Muc.ds.csv')

theme_set(theme_bw())
#sigtabgen = subset(sigtab, !is.na(Genus))
sigtabgen.wax = sigtab.wax %>%
  mutate(Genus = ifelse(Genus == "", paste("Unclassified", Family), Genus)) %>%
  mutate(Family = ifelse(Family == "uncultured", paste("", Class), Family))%>%
  mutate(Genus = ifelse(Genus == "uncultured", paste("Unclassified" , Family), Genus))

# Class order
x = tapply(sigtabgen.wax$log2FoldChange, sigtabgen.wax$Class, function(x) max(x))
x = sort(x, TRUE)
sigtabgen.wax$Class = factor(as.character(sigtabgen.wax$Class), levels=names(x))
# Genus order
x = tapply(sigtabgen.wax$log2FoldChange, sigtabgen.wax$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen.wax$Genus = factor(as.character(sigtabgen.wax$Genus), levels=names(x))

wax.DESeq = ggplot(sigtabgen.wax, aes(y=Genus, x=log2FoldChange, color=Class)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  scale_color_manual(values = c("#A6D854", "#377EB8","#E78AC3","#F781BF","#F0027F","#FDCDAC","#FF7F00","#E41A1C","#BEBADA"))+
  theme(axis.text.x=element_text(size=12,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size=12,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=10, face = 'bold'),
        plot.title = element_text(size = rel(2), face="bold"),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=10)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Differentially abundant ASVs: 60 vs 90 DAE ")
wax.DESeq



##################################################################################################
##### Assign the same colors to the taxa in different plots (e.g. Phyllum in mucilage vs wax)####
###### Unique Family colors #####
unique.Family <- union(unique(Muc_family$Family), unique(Wax_family$Family))
length(unique.Family)
n <- length(unique.Family)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_family = sample(col_vector, n)
names(col_family) <- unique.Family
col_family

Devl.names = c(`60DAE` = "60 Days after emergence", `90DAE` = "90 Days after emergence")

Wax_family_plot <- ggplot(Wax_family,aes(x=Sample,y=Abundance,fill=Family)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1,labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_family) + 
  labs(fill="Bacterial Family")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Family  > 3%) \n") +
  xlab("Epicuticular Wax - Bacterial Microbiome") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 16),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=12),
        legend.title = element_text(size =12, face="bold"),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Family Composition")
Wax_family_plot

Muc_family_plot <- ggplot(Muc_family,aes(x=Sample,y=Abundance,fill=Family)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1,labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_family) + 
  labs(fill="Bacterial Family")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Family  > 3%) \n") +
  xlab("Aerial Root Mucilage - Bacterial Microbiome") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 16),
        strip.text.x = element_text(size=0,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=12),
        legend.title = element_text(size =12, face="bold"),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Muc_family_plot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
bact_legend = g_legend(Wax_family_plot)
fung_legend = g_legend(muc.its.RA.Family)
Figure2 <- cowplot::plot_grid(Wax_family_plot + theme(legend.position = "none"),
                              Muc_family_plot+ theme(legend.position = "none"),
                              muc.its.RA.Family + theme(legend.position = "none"), 
                              ncol=1, align = "v",
                              axis = "lr") 
Fig2_leg = cowplot::plot_grid(bact_legend, fung_legend, ncol=1)
Fig2_wleg = cowplot::plot_grid(Fig2, Fig2_leg, nrow=1, label_y ="Relative abundance (Family > 3%)", rel_widths = c(1,0.4))

Figure2_wleg


######## Unique Genus colors #####################
unique.Genus <- union(unique(Muc_genus$Genus), unique(Wax_genus$Genus))
length(unique.Genus)
n <- length(unique.Genus)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_genus = sample(col_vector, n)
names(col_genus) <- unique.Genus
col_genus

Devl.names = c(`60DAE` = "60 Days after emergence", `90DAE` = "90 Days after emergence")

Wax_genus_plot <- ggplot(Wax_genus,aes(x=Sample,y=Abundance,fill=Genus)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1,labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_genus) + 
  labs(fill="Bacterial Genus")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Genus  > 3%) \n") +
  xlab("Epicuticular Wax - Bacterial Microbiome") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 16),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        legend.title = element_text(size =12, face="bold"),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  #ggtitle("Genus Composition")
Wax_genus_plot

Muc_genus_plot <- ggplot(Muc_genus,aes(x=Sample,y=Abundance,fill=Genus)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1,labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_genus) + 
  labs(fill="Bacterial Genus")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Genus  > 3%) \n") +
  xlab("Aerial Root Mucilage - Bacterial Microbiome") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 16),
        strip.text.x = element_text(size=0,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        legend.title = element_text(size =12, face="bold"),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Muc_genus_plot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
bact_genus_legend = g_legend(Wax_genus_plot)
fung_genus_legend = g_legend(muc.its.RA.Genus)
Supp_Figure1 <- cowplot::plot_grid(Wax_genus_plot + theme(legend.position = "none"),
                              Muc_genus_plot+ theme(legend.position = "none"),
                              muc.its.RA.Genus + theme(legend.position = "none"), 
                              ncol=1, align = "v",
                              axis = "lr") 
Supp_Fig1_leg = cowplot::plot_grid(bact_genus_legend, fung_genus_legend, ncol=1)
Supp_Fig1_wleg = cowplot::plot_grid(Supp_Fig1, Supp_Fig1_leg, nrow=1, label_y ="Relative abundance (Genus > 3%)", rel_widths = c(1,0.4))

Supp_Fig1_with_leg

