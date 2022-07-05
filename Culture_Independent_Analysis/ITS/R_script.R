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

### Set working directory
setwd("/Volumes/ShadeLab/WorkingSpace/MarcoMechan_WorkingSpace/sorghum_16S_its_R_analysis/")

### Import files

muc.its.otu=read.table("its_otu-table.tsv", header = T, row.names = 1)
muc.its.tax=read.delim("its_taxonomy.tsv",row.names = 1)
muc.its.map=read.table("metadata_ITS_20-21.csv", header = T, row.names = 1, sep=",")
muc.its.tree=read_tree("its_rooted_tree.nwk")

################################################################################
############################ Fungal community analysis  ########################
################################################################################

muc.its.tax=muc.its.tax[-2]
muc.its.tax.df <- colsplit(muc.its.tax$Taxon, ';', names =  c("Kingdom", "Phylum", "Class", 
                                                              "Order", "Family", "Genus", "Species"))
muc.its.tax.df[1:7] <- lapply(muc.its.tax.df[1:7], function(x) gsub(".*__", "", x))
rownames(muc.its.tax.df) <- rownames(muc.its.tax)

muc.its.OTU=otu_table(as.matrix(muc.its.otu), taxa_are_rows = T, errorIfNULL = TRUE)
muc.its.TAX=tax_table(as.matrix(muc.its.tax.df), errorIfNULL = TRUE)
muc.its.MAP=sample_data(data.frame(muc.its.map))

muc.its.otuPhyloseq=phyloseq(muc.its.OTU,muc.its.TAX,muc.its.MAP)
sample_sums(muc.its.otuPhyloseq)

summarize_phyloseq(muc.its.otuPhyloseq)
Muc.outphyloseq <- subset_samples(otuPhyloseq, Compartment%in%c("Mucilage")&Year%in%c("2020"))
Muc.outphyloseq <- prune_taxa(taxa_sums(Muc.outphyloseq) > 0, Muc.outphyloseq)
summarize_phyloseq(Muc.outphyloseq)
sample_sums(Muc.outphyloseq)

muc.its.phyloseq <- subset_samples(muc.its.otuPhyloseq, Compartment%in%c("Mucilage")&Year%in%c("2020"))
muc.its.phyloseq <- prune_taxa(taxa_sums(muc.its.phyloseq) > 0, muc.its.phyloseq)
sample_sums(muc.its.phyloseq)
summarize_phyloseq(muc.its.phyloseq)

#Filtering mitochondria, Chloroplast and Unclassified taxa
muc.its.otuPhyloseq.filt <- muc.its.otuPhyloseq %>%
  subset_taxa(Kingdom != "Unassigned")

muc.its.otuPhyloseq.filt <- muc.its.otuPhyloseq.filt %>%
  subset_taxa((Family != "Mitochondria") | is.na(Family))

muc.its.otuPhyloseq.filt <- muc.its.otuPhyloseq.filt %>%
  subset_taxa((Order != "Chloroplast") | is.na(Class))

sample_sums(muc.its.otuPhyloseq.filt) 
print(muc.its.otuPhyloseq.filt) 

summarize_phyloseq(muc.its.otuPhyloseq.filt)


muc.its.phyloseq.filt <- subset_samples(muc.its.otuPhyloseq.filt, Compartment%in%c("Mucilage")&Year%in%c("2021"))
muc.its.phyloseq.filt <- prune_taxa(taxa_sums(muc.its.phyloseq.filt) > 0, muc.its.phyloseq.filt)
sample_sums(muc.its.phyloseq.filt)
summarize_phyloseq(muc.its.phyloseq.filt)


muc.its.filtered.otus <- phyloseq_to_df(muc.its.otuPhyloseq.filt, addtax=T, addtot=T)
sample_sums(muc.its.otuPhyloseq.filt)

######### Rarefaction curves: Chloroplast/Mitochondria/Unassigned filtered data #############
muc.its <- subset_samples(muc.its.otuPhyloseq.filt, Compartment%in%c("Mucilage"))
muc.its <- prune_taxa(taxa_sums(muc.its) > 0, muc.its)

### remove singletons
muc.its <- prune_taxa(taxa_sums(muc.its) > 1, muc.its)
print(muc.its)### 5641taxa and 173 samples
sample_sums(muc.its)

### Samples with less than 30K reads were removed: 2 samples in total were removed
muc.its.1 <- prune_samples(sample_sums(muc.its) >= 30000, muc.its)
muc.its.1 <- prune_taxa(taxa_sums(muc.its.1) > 0, muc.its.1)
print(muc.its.1) ###### 5639 taxa and 171 samples 
sample_sums(muc.its.1)
summarize_phyloseq(muc.its.otuPhyloseq.filt)

### Samples were rarefied to 33975 reads/sample
set.seed(13)
muc.its.rare = rarefy_even_depth(muc.its.1, rngseed=1, sample.size=min(sample_sums(muc.its.1)), replace=F)
sample_sums(muc.its.rare)
print (muc.its.rare) ### 5413 taxa and 171 samples

## Fungal alpha diversity
## Fungal Rarefaction curves
rarecurve.muc.its = rarecurve(t(otu_table(muc.its)), step=1000, label=FALSE)
names(rarecurve.muc.its) <- sample_data(muc.its)$Sample.name

protox.its <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "rare_count"
  mydf$Sample.name <- y
  mydf$Subsample_size <- attr(x, "Subsample")
  mydf
}, x = rarecurve.muc.its, y = as.list(names(rarecurve.muc.its)), SIMPLIFY = FALSE)

rarecurve.muc.its.df = do.call(rbind, protox.its) %>%
  dplyr:::left_join(data.frame(sample_data(muc.its)), by = "Sample.name")

muc.its.rarecurves = ggplot(data=rarecurve.muc.its.df, aes(x=Subsample_size, y=rare_count)) +
  geom_vline(xintercept = 33975, color= "red",  linetype='dashed') +
  scale_color_manual(values=c("#325ddd")) +
  geom_line(aes(group=Sample.name, color=Compartment)) +
  labs(x="Rarefied read depth", y="Fungal ASVs count") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8,angle=0, hjust=0.5),
        axis.text.y = element_text(size = 8),
        axis.title=element_text(size=10,face="bold", vjust = 10),
        legend.position = c(.75, .1),
        legend.title = element_text(size=0),
        legend.text = element_text(colour="#325ddd", size = 12, face = "bold"))+
  ggtitle("C) Fungal alpha diversity")+
  theme(plot.title = element_text(size = 10, face = "bold"))
muc.its.rarecurves
  
### Diversity metrics
muc.its.alpha.df = rbind(data.frame(alpha_measure = specnumber(t(otu_table(muc.its.rare)))) %>%
                           mutate(Index = "Observed taxa", Compartment="Mucilage"),
                         data.frame(alpha_measure = picante::pd(t(otu_table(muc.its.rare)), muc.its.tree)$PD) %>%
                           mutate(Index = "Phylogenetic", Compartment="Mucilage"))

alpha.div.muc.its <- ggplot(muc.its.alpha.df, aes(x=Compartment, y=alpha_measure, color=Compartment)) +
  geom_boxplot(alpha=0.6) +
  scale_color_manual(values=c("#0088dc")) +
  geom_jitter(width=0.15, alpha=0.9) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme_bw() +
  theme(axis.text.x=element_text(size=0,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size=10,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=1, face = 'bold'),
        plot.title = element_text(size = rel(2)), 
        axis.title=element_text(size=12,face="bold", vjust = 10),
        legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~Index, scales = "free_y")

alpha.div.muc.its

its.rich = estimate_richness(muc.its.rare)
muc.its.DAE.obs = pairwise.wilcox.test(its.rich$Observed, sample_data(muc.its.rare)$Developmental.Stage)
muc.its.DAE.Chao1 = pairwise.wilcox.test(its.rich$Chao1, sample_data(muc.its.rare)$Developmental.Stage)
muc.its.DAE.Shan = pairwise.wilcox.test(its.rich$Shannon, sample_data(muc.its.rare)$Developmental.Stage)
muc.its.DAE.Simp = pairwise.wilcox.test(its.rich$Simpson, sample_data(muc.its.rare)$Developmental.Stage)

muc.its.DAE.alpha.wicox = data.frame(index = c("Observed", "Chao1", "Shannon", "Simpson"),
                                     pvalue = c(muc.its.DAE.obs$p.value, muc.its.DAE.Chao1$p.value,
                                                muc.its.DAE.Shan$p.value, muc.its.DAE.Simp$p.value))

muc.its.DAE.alpha.wicox 

muc.its.Y.obs = pairwise.wilcox.test(its.rich$Observed, sample_data(muc.its.rare)$Year)
muc.its.Y.Chao1 = pairwise.wilcox.test(its.rich$Chao1, sample_data(muc.its.rare)$Year)
muc.its.Y.Shan = pairwise.wilcox.test(its.rich$Shannon, sample_data(muc.its.rare)$Year)
muc.its.Y.Simp = pairwise.wilcox.test(its.rich$Simpson, sample_data(muc.its.rare)$Year)

muc.its.Y.alpha.wicox = data.frame(index = c("Observed", "Chao1", "Shannon", "Simpson"),
                                     pvalue = c(muc.its.Y.obs$p.value, muc.its.Y.Chao1$p.value,
                                                muc.its.Y.Shan$p.value, muc.its.Y.Simp$p.value))

muc.its.Y.alpha.wicox 


muc.its.FT.obs = pairwise.wilcox.test(its.rich$Observed, sample_data(muc.its.rare)$Fertilization.Status)
muc.its.FT.Chao1 = pairwise.wilcox.test(its.rich$Chao1, sample_data(muc.its.rare)$Fertilization.Status)
muc.its.FT.Shan = pairwise.wilcox.test(its.rich$Shannon, sample_data(muc.its.rare)$Fertilization.Status)
muc.its.FT.Simp = pairwise.wilcox.test(its.rich$Simpson, sample_data(muc.its.rare)$Fertilization.Status)

muc.its.FT.alpha.wicox = data.frame(index = c("Observed", "Chao1", "Shannon", "Simpson"),
                                    pvalue = c(muc.its.FT.obs$p.value, muc.its.FT.Chao1$p.value,
                                               muc.its.FT.Shan$p.value, muc.its.FT.Simp$p.value))

muc.its.FT.alpha.wicox 

muc.its.y.obs = pairwise.wilcox.test(its.rich$Observed, sample_data(muc.its.rare)$Year)
muc.its.y.Chao1 = pairwise.wilcox.test(its.rich$Chao1, sample_data(muc.its.rare)$Year)
muc.its.y.Shan = pairwise.wilcox.test(its.rich$Shannon, sample_data(muc.its.rare)$Year)
muc.its.y.Simp = pairwise.wilcox.test(its.rich$Simpson, sample_data(muc.its.rare)$Year)

muc.its.y.alpha.wicox = data.frame(index = c("Observed", "Chao1", "Shannon", "Simpson"),
                                   pvalue = c(muc.its.y.obs$p.value, muc.its.y.Chao1$p.value,
                                              muc.its.y.Shan$p.value, muc.its.y.Simp$p.value))

muc.its.y.alpha.wicox 


alpha_meas = c("Observed", "Shannon", "Simpson")
(p_muc <- plot_richness(muc.its.rare, x= "Year", measures=alpha_meas, color="Year"))
p_muc + geom_boxplot(alpha=0.6) +
  geom_jitter(aes(color=Year), width=0.15, alpha=0.9)+
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme(axis.text.x=element_text(size=14,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10))


### Merge Rarefaction curves and Div metrics
Figure1c <- cowplot::plot_grid(muc.its.rarecurves, alpha.div.muc.its, ncol=1, rel_heights = c(1, 1), align = "v",
                               rel_widths = c(1,1))

Figure1c

## Fungal Beta Diversity: Bray-Curtis dissimilarity 
bray.diss.muc.its = phyloseq::distance(muc.its.rare, method="bray")
ordination.muc.its = ordinate(muc.its.rare, method="PCoA", distance=bray.diss.muc.its)

#ADONIS test
muc.its.FT.adon<- adonis2(bray.diss.muc.its ~ sample_data(muc.its.rare)$Fertilization.Status)
its.ft.pVal=muc.its.FT.adon$aov.tab[6]$`Pr(>F)`[1]
its.ft.r2 = round(muc.its.FT.adon$R2[1], digits = 4)

muc.its.DAE.adon <- adonis2(bray.diss.muc.its ~ sample_data(muc.its.rare)$Developmental.Stage)
its.dae.pVal=muc.its.DAE.adon$`Pr(>F)`[1]
r2 = round(muc.its.DAE.adon$R2[1], digits = 4)

muc.its.year.adon <- adonis2(bray.diss.muc.its ~ sample_data(muc.its.rare)$Year)
year.pVal=muc.its.year.adon$`Pr(>F)`[1]
r2 = round(muc.its.year.adon$R2[1], digits = 4)

Muc.its.DAE.FS.adonOut<- adonis2(bray.diss.muc.its ~ sample_data(muc.its.rare)$Developmental.Stage*sample_data(muc.its.rare)$Fertilization.Status)
Muc.its.DAE.Year.adonOut<- adonis2(bray.diss.muc.its ~ sample_data(muc.its.rare)$Developmental.Stage*sample_data(muc.its.rare)$Year)
Muc.its.FS.year.adonOut<- adonis2(bray.diss.muc.its ~ sample_data(muc.its.rare)$Fertilization.Status*sample_data(muc.its.rare)$Year)


muc.its.PCoA = plot_ordination(muc.its.rare, ordination.muc.its, color = 'Year', shape = 'Developmental.Stage') + theme(aspect.ratio=1)+
  geom_point(na.rm = FALSE, size=4) + 
  scale_shape_manual(values = c(19, 1))+
  scale_color_manual(values=c("#314293", "#7692e9")) +
  labs(x="PCoA1 (51.6% var. explained)", y="PCoA1 (5.1% var. explained)") +
  theme_bw() +
  theme(axis.text.x=element_text(size=12,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size=12,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=12, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=12,face="bold", vjust = 10),
        legend.text = element_text(size=10),
        legend.title = element_text(size =12, face="bold"))+
  labs(colour="Year of sampling", shape="Developmental Stage")+
  annotate(geom = 'text', label=paste("PERMANOVA\nR2=", r2, "\nP-val=",year.pVal, sep=''), x=.05, y=-.2)+
  ggtitle("Fungal Root Mucilage")

muc.its.PCoA$layers <- muc.its.PCoA$layers[-1]
muc.its.PCoA
Figure3c <- muc.its.PCoA

#### Betadisper 
muc.its.dispr.year <- betadisper(bray.diss.muc.its, phyloseq::sample_data(muc.its.rare)$Year)
muc.its.dispr.year
## Tukey's Honest Significant Differences
(muc.its.dispr.year.HSD <- TukeyHSD(muc.its.dispr.year))
plot(muc.its.dispr.year.HSD)
## Plot the groups and distances to centroids on the first two PCoA axes
plot(muc.its.dispr.year)
## Draw a boxplot of the distances to centroid for each group
boxplot(muc.its.dispr.year)

#### Betadisper 
muc.its.fs.dispr <- betadisper(bray.diss.muc.its, phyloseq::sample_data(muc.its.rare)$Fertilization.Status)
muc.its.fs.dispr
## Tukey's Honest Significant Differences
(muc.its.fs.dispr.HSD <- TukeyHSD(muc.its.fs.dispr))
plot(muc.its.fs.dispr.HSD)
## Plot the groups and distances to centroids on the first two PCoA axes
plot(muc.its.fs.dispr)
## Draw a boxplot of the distances to centroid for each group
boxplot(muc.its.fs.dispr)

#### Betadisper 
muc.its.ds.dispr <- betadisper(bray.diss.muc.its, phyloseq::sample_data(muc.its.rare)$Developmental.Stage)
muc.its.ds.dispr
## Tukey's Honest Significant Differences
(muc.its.ds.dispr.HSD <- TukeyHSD(muc.its.ds.dispr))
plot(muc.its.ds.dispr.HSD)
## Plot the groups and distances to centroids on the first two PCoA axes
plot(muc.its.ds.dispr)
## Draw a boxplot of the distances to centroid for each group
boxplot(muc.its.ds.dispr)


## Fungal Relative abundance
muc.its.phylum <- muc.its.rare %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.001) %>%
  arrange(Phylum)

write.csv(muc.its.phylum, 'muc.its.phylum.RA.phylum.csv') 

n <- dim(muc.its.phylum)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

muc.its.RA.Phylum = ggplot(muc.its.phylum,aes(x=Sample,y=Abundance,fill=Phylum, order=as.factor(Phylum))) +
  facet_wrap(~Year, scales= "free_x", nrow=1) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Fungal Phylum")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Phylum) \n") +
  xlab("Sorghum aerial root mucilage") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=12),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Fungal Phylum Composition")
muc.its.RA.Phylum

muc.its.class <- muc.its.rare %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.001) %>%
  arrange(Class)
write.csv(muc.its.class, 'muc.its.RA.class.csv')

n <- dim(muc.its.class)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

muc.its.RA.Class = ggplot(muc.its.class,aes(x=Sample,y=Abundance,fill=Class)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1)+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Fungal Class")+
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
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
muc.its.RA.Class

muc.its.family <- muc.its.rare %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.02) %>%
  arrange(Family)
write.csv(muc.its.family, 'muc.its.RA.family.csv')

Devl.names = c(`60DAE` = "60 Days after emergence", `90DAE` = "90 Days after emergence")

n <- dim(muc.its.family)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

muc.its.RA.Family <- ggplot(muc.its.family,aes(x=Sample,y=Abundance,fill=Family)) +
  facet_wrap(~Developmental.Stage, scales= "free_x", nrow=1, labeller = as_labeller(Devl.names))+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Fungal Family")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Family > 3%) \n") +
  xlab("Aerial Root Mucilage - Fungal Microbiome") +
  theme(axis.text.x=element_text(size=0,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 16),
        strip.text.x = element_text(size=0,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title = element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=12),
        legend.title = element_text(size =12, face="bold"),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  #ggtitle("Family Composition")
muc.its.RA.Family

muc.its.genus <- muc.its.rare %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.03) %>%
  arrange(Genus)
#write.csv(muc.its.genus, 'muc.its.RA.genus.csv')

n <- dim(muc.its.genus)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

muc.its.RA.Genus = ggplot(muc.its.genus,aes(x=Sample,y=Abundance,fill=Genus)) +
  facet_wrap(~Year, scales= "free_x", nrow=1)+
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  labs(fill="Fungal Genus")+
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Genus > 2%) \n") +
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
        panel.grid.minor = element_blank()) 
  #ggtitle("Fungal Genus Composition")
muc.its.RA.Genus
