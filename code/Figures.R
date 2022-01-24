### Wombat GI tract 16S paper 2021
### Raphael Eisenhofer

library(phyloseq)
library(forcats)
library(qiime2R)
library(dplyr)
library(svglite)
library(cowplot)
library(grid)
library(gtable)
library(gplots)
library(ggplot2)
library(scales)
library(ggpubr)
library(ggrepel)
library(tidyr)
library(knitr)
library(microbiome)
library(microbiomeutilities)
library(VennDiagram)
library(ggVennDiagram)

## Import data

ps<-qza_to_phyloseq(
  features="data/QIIME2_output/Wombat-GIT-table.qza",
  tree="data/QIIME2_output/Wombat-GIT-sepp-tree.qza",
  taxonomy="data/QIIME2_output/Wombat-GIT-SILVA-138.qza",
  metadata = "data/QIIME2_output/Wombat-GIT-Metadata.txt"
)

###############
##Filter data##
###############

#Remove negative control for subsequent analyses:
psUF <- subset_samples(ps, sampletype != "negative control")

#What phyla contain the most ASVs?
##Code modified from Callahan et al. 2016: https://doi.org/10.12688/f1000research.8986.2
sort(table(tax_table(psUF)[,"Phylum"], exclude=NULL))

#Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(psUF),
             MARGIN = ifelse(taxa_are_rows(psUF), yes = 1, no = 2),
             FUN = function(x){sum(x > 0)})
#Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                   TotalAbundance = taxa_sums(psUF),
                   tax_table(psUF))

#OK, how many ASVs are only found in 1 replicate (we have replicates for each sample).
count(prevdf, vars = Prevalence <= 1)

#Looks like only 105/2524 ASVs. We'll still remove them. What % of samples is 1 replicate?
nsamples(psUF)
#30 samples
1/30
#1/30 = 0.0333

#Plot prevalance/abundance of each ASV, grouped by phylum, with the proposed prevalence
# filter threshold line
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(psUF,"Phylum"))

ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(psUF), color = Phylum)) +
  geom_hline(yintercept=0.0333, alpha=0.5, linetype=2) + geom_point(size=2, alpha=0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

ggsave("output/AbundancePrevalence_Per_Phylum.png",
       width = 19, height = 10, dpi = 300)

#Execute prevalence filter
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= 2)]
psUF.prevF <- prune_taxa(keepTaxa, psUF)


#####################
###Alpha diversity###
#####################

## Figure 1

#Set seed, check the sample with the lowest read count
set.seed(1337)

sort(sample_sums(psUF.prevF))

#Rarefy n number of reads where n = read count of samples with fewest reads
psUF.prevF.rar <- rarefy_even_depth(psUF.prevF, sample.size = 21710)

#Copmute alpha diversity metrics
psUF.prevF.alpha.div <- alpha(psUF.prevF.rar, index = "all")

#Append metadata to alpha diversity values
psUF.prevF.meta <- meta(psUF.prevF.rar)
psUF.prevF.meta$name <- rownames(psUF.prevF.meta)
psUF.prevF.alpha.div$name <- rownames(psUF.prevF.alpha.div)
psUF.prevF.alpha.div.df <- merge(psUF.prevF.alpha.div, psUF.prevF.meta, by = "name")

#Split dataframe by species
BNW.alpha.div.df <- filter(psUF.prevF.alpha.div.df, species == "Bare-nosed")
SHNW.alpha.div.df <- filter(psUF.prevF.alpha.div.df, species == "Hairy-nosed")

##GGPLOT!
#Set colours
colours <- c("#ae017e", "#f768a1", "#fbb4b9", "#feebe2")

#Reorder to match GI tract order
BNW.alpha.div.df.ordered <- BNW.alpha.div.df
BNW.alpha.div.df.ordered$sampleshort <- factor(BNW.alpha.div.df.ordered$sampleshort,
                                               levels = c('ST', 'SI', 'PC1', 'PC2',
                                                          'PC3', 'PC4', 'DC'))

SHNW.alpha.div.df.ordered <- SHNW.alpha.div.df
SHNW.alpha.div.df.ordered$sampleshort <- factor(SHNW.alpha.div.df.ordered$sampleshort,
                                                levels = c('ST', 'PSI', 'DSI', 'PC1',
                                                           'PC2', 'PC3', 'DC1', 'DC3'))

BNW.richness <- ggplot(BNW.alpha.div.df.ordered,
                       aes(x=sampleshort, y=observed, fill=sampletype))
SHNW.richness <- ggplot(SHNW.alpha.div.df.ordered,
                        aes(x=sampleshort, y=observed, fill=sampletype))

#Bare-nosed plot
gg.BNW.richness = BNW.richness +
  #Boxplot
  facet_wrap(~sampleshort, scales = "free_x", ncol = 7) +
  geom_line() +
  #Jitter, size, colour
  geom_jitter(position=position_dodge2(0.1), size=8, aes(colour=sampletype)) +
  #Custom manual colours
  scale_colour_manual(values=colours) +
  scale_fill_manual(values=colours) +
  #Tick labels
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold", size=20),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        #Legend
        legend.position = "none") +
  #Axis labels
  scale_y_continuous(breaks=seq(0,900,100)) +
  labs(x = "") +
  labs(y = "ASV Richness\n")

#Adjust sizes of facets to match proportions of the drawing
#Using 'grid' to adjust sizes. To figure out which element to select, use "gtable_show_layout(gt)"
#Also save as SVG file

svg(filename = "output/Figure1A.svg", width = 20, height = 6)

gtbnw = ggplot_gtable(ggplot_build(gg.BNW.richness))
gtbnw$widths[5] = 0.5*gtbnw$widths[5]
gtbnw$widths[9] = 5.1*gtbnw$widths[9]
gtbnw$widths[13] = 1.5*gtbnw$widths[13]
gtbnw$widths[17] = 1.75*gtbnw$widths[17]
gtbnw$widths[21] = 1.25*gtbnw$widths[21]
gtbnw$widths[25] = 1.25*gtbnw$widths[25]
gtbnw$widths[29] = 2.25*gtbnw$widths[29]
grid.draw(gtbnw)

dev.off()


#add jitterm, size, and colour
gg.SHNW.richness = SHNW.richness +
  #Boxplot
  facet_wrap(~sampleshort, scales = "free_x", ncol = 8) +
  geom_line() +
  #Jitter, size, colour
  geom_jitter(position=position_dodge2(0.1), size=8, aes(colour=sampletype)) +
  #Custom manual colours
  scale_colour_manual(values=colours) +
  scale_fill_manual(values=colours) +
  #Tick labels
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold", size=20),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        #Legend
        legend.position = "none") +
  #Axis labels
  scale_y_continuous(breaks=seq(0,900,100)) +
  labs(x = "") +
  labs(y = "ASV Richness\n")

svg(filename = "output/Figure1B.svg", width = 20, height = 6)

gtshnw = ggplot_gtable(ggplot_build(gg.SHNW.richness))
gtshnw$widths[5] = 0.5*gtshnw$widths[5]
gtshnw$widths[9] = 3.1*gtshnw$widths[9]
gtshnw$widths[13] = 3.1*gtshnw$widths[13]
gtshnw$widths[17] = 1*gtshnw$widths[17]
gtshnw$widths[21] = 0.35*gtshnw$widths[21]
gtshnw$widths[25] = 0.8*gtshnw$widths[25]
gtshnw$widths[27] = 10*gtshnw$widths[27]
gtshnw$widths[29] = 2.45*gtshnw$widths[29]
gtshnw$widths[33] = 2.25*gtshnw$widths[33]
grid.draw(gtshnw)

dev.off()


#####################
###Beta diversity####
#####################

## Figure 2


#Unweighted UniFrac
ord.unw.uni <- ordinate(psUF.prevF.rar, "PCoA", "unifrac", weighted=F)

#Weighted UniFrac
ord.w.uni <- ordinate(psUF.prevF.rar, "PCoA", "unifrac", weighted=T)

#Axes 1/2
unwt.unifrac.1.2 <- plot_ordination(psUF.prevF.rar,
                                    ord.unw.uni, color="sampletype", shape = "species",
                                    axes = c(1, 2))
w.unifrac.1.2 <- plot_ordination(psUF.prevF.rar,
                                 ord.w.uni, color="sampletype", shape = "species",
                                 axes = c(1, 2))
#Axes 1/3
unwt.unifrac.1.3 <- plot_ordination(psUF.prevF.rar,
                                    ord.unw.uni, color="sampletype", shape = "species",
                                    axes = c(1, 3))
w.unifrac.1.3 <- plot_ordination(psUF.prevF.rar,
                                 ord.w.uni, color="sampletype", shape = "species",
                                 axes = c(1, 3))

#Plot it
Fig2a <- unwt.unifrac.1.2 +
  scale_x_reverse() +
  geom_point(size=9) +
  geom_path(size = 1) +
  geom_text_repel(aes(label=samplepoint), color="black", fontface="bold", size=5) +
  scale_colour_manual(values=colours) +
  guides(color = guide_legend(reverse=TRUE)) +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

Fig2b <- w.unifrac.1.2 +
  scale_x_reverse() +
  geom_point(size=9) +
  geom_path(size = 1) +
  geom_text_repel(aes(label=samplepoint), color="black", fontface="bold", size=5) +
  geom_rect(mapping = aes(xmin=-0.06, xmax=-0.14, ymin=-0.13, ymax=-0.20),
            color = "black", alpha = 0, size = 2) +
  geom_text(aes(label="First proximal\n colon samples", x=-0, y=-0.163),
            size = 5, color = "black") +
  scale_colour_manual(values=colours) +
  guides(color = guide_legend(reverse=TRUE)) +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

#Create figure 2:

ggarrange(Fig2a, Fig2b, nrow = 2, common.legend = TRUE, legend="bottom",
          labels = c("A)", "B)"), font.label = list(size=20, face="bold", color="black"))

ggsave(filename = "output/Figure4.png", width = 10, height = 10, dpi = 300)
ggsave(filename = "output/Figure4.svg", width = 10, height = 10, dpi = 300)
ggsave(filename = "output/Figure4.pdf", width = 10, height = 10, dpi = 300)


## Supplementary figure 1
## Remove the stomach/SI samples to get a better look at what's going on in the colon!
ps.colon <- subset_samples(psUF.prevF, sampletype != "Stomach" &
                             sampletype != "Small intestine")

#Rarefy n number of reads where n = read count of samples with fewest reads
sort(sample_sums(ps.colon))
ps.colon.rar <- rarefy_even_depth(ps.colon, sample.size = 47301)

#Unweighted UniFrac
ord.unw.uni.colon <- ordinate(ps.colon.rar, "PCoA", "unifrac", weighted=F)

#Weighted UniFrac
ord.w.uni.colon <- ordinate(ps.colon.rar, "PCoA", "unifrac", weighted=T)

#Axes 1/2
unwt.unifrac.1.2.colon <- plot_ordination(ps.colon.rar, ord.unw.uni.colon,
                                          color="sampletype", shape = "species",
                                          axes = c(1, 2))
w.unifrac.1.2.colon <- plot_ordination(ps.colon.rar, ord.w.uni.colon,
                                       color="sampletype", shape = "species",
                                       axes = c(1, 2))
#Axes 1/3
unwt.unifrac.1.3.colon <- plot_ordination(ps.colon.rar, ord.unw.uni.colon,
                                          color="sampletype", shape = "species",
                                          axes = c(1, 3))
w.unifrac.1.3.colon <- plot_ordination(ps.colon.rar, ord.w.uni.colon,
                                       color="sampletype", shape = "species",
                                       axes = c(1, 3))

#Plot it
SI_figure_1A <- unwt.unifrac.1.2.colon +
  #stat_ellipse() +
  scale_x_reverse() +
  geom_point(size=9) +
  geom_text(aes(label=samplepoint), color="black", hjust=2.5, vjust=-1) +
  geom_path(size = 1) +
  scale_colour_manual(values=colours) +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

SI_figure_1B <- w.unifrac.1.2.colon +
  #stat_ellipse() +
  scale_x_reverse() +
  geom_point(size=9) +
  geom_text(aes(label=samplepoint), color="black", hjust=2.5, vjust=-1) +
  geom_path(size = 1) +
  scale_colour_manual(values=colours) +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

#Create supplmentary figure 1:

ggarrange(SI_figure_1A, SI_figure_1B, nrow = 2, common.legend = TRUE, legend="right",
          labels = c("A)", "B)"), font.label = list(size=20, face="bold", color="black"))

ggsave(filename = "output/SI_Figure1.png", width = 10, height = 10, dpi = 300)
ggsave(filename = "output/SI_Figure1.svg", width = 10, height = 10, dpi = 300)
ggsave(filename = "output/SI_Figure1.pdf", width = 10, height = 10, dpi = 300)


#####################
###Taxa bar plots####
#####################

## I want stacked bar plots for each wombat species, two separate plots using the same
## colours for each taxa

## Supplementary figure 2

#Filter out the top 10 most abundant phlya to reduce plot colour needs
psUF.phylum <- tax_glom(psUF, taxrank = "Phylum")
psUF.phylum.10 <- prune_taxa(names(sort(taxa_sums(psUF.phylum),TRUE)[1:10]), psUF.phylum)

#What % of reads do the top 10 phyla have?
print("mean:")
mean(sort(sample_sums(psUF.phylum.10)/sample_sums(psUF.phylum))) %>% print()
print("SD:")
sd(sort(sample_sums(psUF.phylum.10)/sample_sums(psUF.phylum))) %>% print()
#Looks like the top 10 phyla account for the vast majority of sequences

#Split into separate phyloseq objects by species
psUF.phylum.10.shnw <- subset_samples(psUF.phylum.10, species == "Hairy-nosed")
psUF.phylum.10.bnw <- subset_samples(psUF.phylum.10, species == "Bare-nosed")

#Merge replicate samples to look at the average of each site
psUF.phylum.10.shnw.merged <- merge_samples(psUF.phylum.10.shnw, "sampleshort2")
psUF.phylum.10.bnw.merged <- merge_samples(psUF.phylum.10.bnw, "sampleshort2")

#Convert to compositional (relative abundance)
psUF.phylum.10.shnw.merged.rel <- microbiome::transform(psUF.phylum.10.shnw.merged, "compositional")
psUF.phylum.10.bnw.merged.rel <- microbiome::transform(psUF.phylum.10.bnw.merged, "compositional")

#Melt into a dataframe
pd.shnw <- psmelt(psUF.phylum.10.shnw.merged.rel)
pd.bnw <- psmelt(psUF.phylum.10.bnw.merged.rel)

#Set a protanomaly-friendly colour palette
colour_phylum <- c("Actinobacteriota" = "#EE6677", "Bacteroidota" = "#228833",
                   "Cyanobacteria" = "#4477AA", "Desulfobacterota" = "#CCBB44",
                   "Fibrobacterota" = "#66CCEE", "Firmicutes" = "#AA3377",
                   "Fusobacteriota" = "#24ff24", "Proteobacteria" = "#db6d00",
                   "Spirochaetota" = "#b66dff", "Verrucomicrobiota" = "#000000")

#plot the suckers
tax.phylum.bnw <- ggplot(pd.bnw, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", size = 0) +
  facet_grid(~samplecollection, scales = "free_x") +
  labs(x = "", y = "Relative abundance") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.y = element_text(size=20, face = 'bold'),
        axis.title.y = element_text(size=20, face = 'bold'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 25, face = 'bold'),
        legend.position = "bottom",
        panel.spacing = unit(2, "lines"),
        panel.background = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = colour_phylum)

#Add spacing between bars to match the GI tract drawing
svg(filename = "output/SI_Figure2A.svg", width = 20, height = 6)

gtbnw.tax = ggplot_gtable(ggplot_build(tax.phylum.bnw))
#Bars
gtbnw.tax$widths[5] = 0.2*gtbnw.tax$widths[5]
gtbnw.tax$widths[7] = 2*gtbnw.tax$widths[7]
gtbnw.tax$widths[9] = 0.7*gtbnw.tax$widths[9]
gtbnw.tax$widths[11] = 0.7*gtbnw.tax$widths[11]
gtbnw.tax$widths[13] = 0.5*gtbnw.tax$widths[13]
gtbnw.tax$widths[15] = 0.5*gtbnw.tax$widths[15]
gtbnw.tax$widths[17] = 1*gtbnw.tax$widths[17]
#Bar spacing
gtbnw.tax$widths[6] = 2*gtbnw.tax$widths[6]
gtbnw.tax$widths[8] = 0.5*gtbnw.tax$widths[8]
gtbnw.tax$widths[10] = 0*gtbnw.tax$widths[10]
gtbnw.tax$widths[12] = 0.5*gtbnw.tax$widths[12]
gtbnw.tax$widths[14] = 0.5*gtbnw.tax$widths[14]
gtbnw.tax$widths[16] = 0.5*gtbnw.tax$widths[16]

#gtbnw.tax$heights[11,11] = 10*gtbnw.tax$heights[11,11]
grid.draw(gtbnw.tax)

dev.off()


tax.phylum.shnw <- ggplot(pd.shnw, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", size = 0) +
  facet_grid(~samplecollection, scales = "free_x") +
  labs(x = "", y = "Relative abundance") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.y = element_text(size=20, face = 'bold'),
      axis.title.y = element_text(size=20, face = 'bold'),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(size = 1),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 25, face = 'bold'),
      legend.position = "bottom",
      panel.spacing = unit(2, "lines"),
      panel.background = element_blank(),
      axis.text.x = element_blank()) +
  scale_fill_manual(values = colour_phylum)

#Add spacing between bars to match the GI tract drawing
svg(filename = "output/SI_Figure2B.svg", width = 20, height = 6)

gtshnw.tax = ggplot_gtable(ggplot_build(tax.phylum.shnw))
#Bars
gtshnw.tax$widths[5] = 0.5*gtshnw.tax$widths[5]
gtshnw.tax$widths[7] = 2.25*gtshnw.tax$widths[7]
gtshnw.tax$widths[9] = 2.25*gtshnw.tax$widths[9]
gtshnw.tax$widths[11] = 1*gtshnw.tax$widths[11]
gtshnw.tax$widths[13] = 0.35*gtshnw.tax$widths[13]
gtshnw.tax$widths[15] = 0.8*gtshnw.tax$widths[15]
gtshnw.tax$widths[17] = 2.25*gtshnw.tax$widths[17]
gtshnw.tax$widths[19] = 2.25*gtshnw.tax$widths[19]
#Bar spacing
gtshnw.tax$widths[6] = 0.75*gtshnw.tax$widths[6]
gtshnw.tax$widths[8] = 2*gtshnw.tax$widths[8]
gtshnw.tax$widths[10] = 1.5*gtshnw.tax$widths[10]
gtshnw.tax$widths[12] = 0.1*gtshnw.tax$widths[12]
gtshnw.tax$widths[14] = 0.1*gtshnw.tax$widths[14]
gtshnw.tax$widths[16] = 1.2*gtshnw.tax$widths[16]
gtshnw.tax$widths[18] = 0.001*gtshnw.tax$widths[18]
gtshnw.tax$widths[20] = 3*gtshnw.tax$widths[20]

#gtshnw.tax$heights[11,11] = 10*gtshnw.tax$heights[11,11]
grid.draw(gtshnw.tax)

dev.off()


## Collapse taxonomy to family level
## Print relative abundance values of the different families

psUF.family <- tax_glom(psUF, taxrank = "Family", NArm = FALSE)

psUF.family.merged <- merge_samples(psUF.family, "sampleshort2")

sam_data_fixed <- psUF.family@sam_data[- grep("B", psUF.family@sam_data$replicate),]
rownames(sam_data_fixed) <- gsub("-A", "", rownames(sam_data_fixed))
sam_data_fixed
psUF.family.merged@sam_data = sam_data_fixed
#Had to replace sample_data, as the merge_sample function wrecks it...

psUF.family.merged.relab <- transform_sample_counts(psUF.family.merged, function(x) x / sum(x))
psUF.family.merged.relab@otu_table <- phyloseq::t(psUF.family.merged.relab@otu_table)
#Had to transpose otu table, as merge_samples flips it for some reason...

psUF.family.merged.relab.tax <- as.data.frame(tax_table(psUF.family.merged.relab))
psUF.family.merged.relab.otu.and.tax <- cbind(psUF.family.merged.relab.tax,
                                              psUF.family.merged.relab@otu_table)

DT::datatable(psUF.family.merged.relab.otu.and.tax)

write.csv(psUF.family.merged.relab.otu.and.tax, file = "output/Taxonomic_classifications_rel_abu.csv")


## Figure 2
## Taxa bar plots at family level. I want to colour only the top 20-most abundant.

#First, figure out what the 20-most abundant families are, then create a column to
# separate them.
top20families = names(sort(taxa_sums(psUF.family.merged), TRUE)[1:20])
taxtab20 = cbind(tax_table(psUF.family.merged), family_20 = NA)
taxtab20[top20families, "family_20"] <- as(tax_table(psUF.family.merged)
                                           [top20families, "Family"],
                                           "character")
tax_table(psUF.family.merged) <- tax_table(taxtab20)

#Transform to relative abundance
psUF.family.merged <- transform_sample_counts(psUF.family.merged, function(x) 100 * x/sum(x))

#What is the average total relative abundance of the top 20-most abundant ASVs?
mean(sample_sums(prune_taxa(top20families, psUF.family.merged)))

#Split into separate phyloseq objects by species (for each panel of figure)
psUF.family.merged.shnw <- subset_samples(psUF.family.merged, species == "Hairy-nosed")
psUF.family.merged.bnw <- subset_samples(psUF.family.merged, species == "Bare-nosed")

#Melt into a dataframe
pd.family.shnw <- psmelt(psUF.family.merged.shnw)
pd.family.bnw <- psmelt(psUF.family.merged.bnw)

#Replace NA with 'other', to help assign colours
pd.family.shnw.sorted <- arrange(pd.family.shnw, family_20)
pd.family.bnw.sorted <- arrange(pd.family.bnw, family_20)
pd.family.shnw.sorted$family_20[is.na(pd.family.shnw.sorted$family_20)] <- c("Other")
pd.family.bnw.sorted$family_20[is.na(pd.family.bnw.sorted$family_20)] <- c("Other")

#Set a protanomaly-friendly 20-colour palette

colour_family <- c("Bacteroidaceae" = "#228B22", "Bacteroidales_BS11_gut_group" = "#3CB371",
                   "Bacteroidales_RF16_group" = "#2E8B57", "Bacteroidales_UCG-001" = "#8FBC8F",
                   "Christensenellaceae" = "#191970", "Clostridiaceae" = "#db6d00",
                   "Fibrobacteraceae" = "#66CCEE", "Fusobacteriaceae" = "#FF1493",
                   "Izemoplasmatales" = "#0000FF", "Lachnospiraceae" = "#B8860B",
                   "Oscillospiraceae" = "#FF6347", "p-251-o5" = "#F0E68C",
                   "Pasteurellaceae" = "#FF69B4", "Peptostreptococcaceae" = "#FFC0CB",
                   "Prevotellaceae" = "#FFFF00", "Rikenellaceae" = "#00FF7F",
                   "Spirochaetaceae" = "#b66dff", "Streptococcaceae" = "#FF00FF",
                   "Veillonellaceae" = "##8B0000", "WCHB1-41" = "#000000", "Other" = "#808080")

#Plot the suckers
pd.family.bnw.plot <- ggplot(pd.family.bnw.sorted,
                             aes(x = Sample, y = Abundance,
                                 fill = fct_reorder(family_20, -Abundance))) +
  geom_bar(stat = "identity", size = 0) +
  facet_grid(~samplecollection, scales = "free_x") +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.y = element_text(size=20, face = 'bold'),
        axis.title.y = element_text(size=20, face = 'bold'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"),
        panel.background = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = colour_family)

#Add spacing between bars to match the GI tract drawing
svg(filename = "output/Figure3A.svg", width = 20, height = 6)

gtbnw.tax.family = ggplot_gtable(ggplot_build(pd.family.bnw.plot))
#Bars
gtbnw.tax.family$widths[5] = 0.2*gtbnw.tax.family$widths[5]
gtbnw.tax.family$widths[7] = 2*gtbnw.tax.family$widths[7]
gtbnw.tax.family$widths[9] = 0.7*gtbnw.tax.family$widths[9]
gtbnw.tax.family$widths[11] = 0.7*gtbnw.tax.family$widths[11]
gtbnw.tax.family$widths[13] = 0.5*gtbnw.tax.family$widths[13]
gtbnw.tax.family$widths[15] = 0.5*gtbnw.tax.family$widths[15]
gtbnw.tax.family$widths[17] = 1*gtbnw.tax.family$widths[17]
#Bar spacing
gtbnw.tax.family$widths[6] = 2*gtbnw.tax.family$widths[6]
gtbnw.tax.family$widths[8] = 0.5*gtbnw.tax.family$widths[8]
gtbnw.tax.family$widths[10] = 0*gtbnw.tax.family$widths[10]
gtbnw.tax.family$widths[12] = 0.5*gtbnw.tax.family$widths[12]
gtbnw.tax.family$widths[14] = 0.5*gtbnw.tax.family$widths[14]
gtbnw.tax.family$widths[16] = 0.5*gtbnw.tax.family$widths[16]

#gtbnw.tax$heights[11,11] = 10*gtbnw.tax$heights[11,11]
grid.draw(gtbnw.tax.family)

dev.off()


pd.family.shnw.plot <- ggplot(pd.family.shnw.sorted,
                              aes(x = Sample, y = Abundance,
                                  fill = fct_reorder(family_20, -Abundance))) +
  geom_bar(stat = "identity", size = 0) +
  facet_grid(~samplecollection, scales = "free_x") +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.y = element_text(size=20, face = 'bold'),
        axis.title.y = element_text(size=20, face = 'bold'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"),
        panel.background = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = colour_family)

#Add spacing between bars to match the GI tract drawing
svg(filename = "output/Figure3B.svg", width = 20, height = 6)

gtshnw.tax.family = ggplot_gtable(ggplot_build(pd.family.shnw.plot))
#Bars
gtshnw.tax.family$widths[5] = 0.5*gtshnw.tax.family$widths[5]
gtshnw.tax.family$widths[7] = 2.25*gtshnw.tax.family$widths[7]
gtshnw.tax.family$widths[9] = 2.25*gtshnw.tax.family$widths[9]
gtshnw.tax.family$widths[11] = 1*gtshnw.tax.family$widths[11]
gtshnw.tax.family$widths[13] = 0.35*gtshnw.tax.family$widths[13]
gtshnw.tax.family$widths[15] = 0.8*gtshnw.tax.family$widths[15]
gtshnw.tax.family$widths[17] = 2.25*gtshnw.tax.family$widths[17]
gtshnw.tax.family$widths[19] = 2.25*gtshnw.tax.family$widths[19]
#Bar spacing
gtshnw.tax.family$widths[6] = 0.75*gtshnw.tax.family$widths[6]
gtshnw.tax.family$widths[8] = 2*gtshnw.tax.family$widths[8]
gtshnw.tax.family$widths[10] = 1.5*gtshnw.tax.family$widths[10]
gtshnw.tax.family$widths[12] = 0.1*gtshnw.tax.family$widths[12]
gtshnw.tax.family$widths[14] = 0.1*gtshnw.tax.family$widths[14]
gtshnw.tax.family$widths[16] = 1.2*gtshnw.tax.family$widths[16]
gtshnw.tax.family$widths[18] = 0.001*gtshnw.tax.family$widths[18]
gtshnw.tax.family$widths[20] = 3*gtshnw.tax.family$widths[20]

#gtshnw.tax.family$heights[11,11] = 10*gtshnw.tax.family$heights[11,11]
grid.draw(gtshnw.tax.family)

dev.off()


#Now for SI figure 4, heatmap of all families
## Collapse taxonomy to family level
## Print relative abundance values of the different families

psUF.family <- tax_glom(psUF, taxrank = "Family", NArm = FALSE)

psUF.family.merged <- merge_samples(psUF.family, "sampleshort2")

#Fix sample data
sam_data_fixed <- psUF.family@sam_data[- grep("B", psUF.family@sam_data$replicate),]
rownames(sam_data_fixed) <- gsub("-A", "", rownames(sam_data_fixed))
sam_data_fixed
psUF.family.merged@sam_data = sam_data_fixed

#Get best hit:
psUF.family.merged.besthit <- microbiomeutilities::format_to_besthit(
  psUF.family.merged, prefix = NULL)

#Create new column with Phylum+Family for ambiguous IDs:
new_tax_table_fam <- as.data.frame(tax_table(psUF.family.merged.besthit)) %>%
  unite("Phylum_Family", c("Phylum", "Family"), sep = " -- ", remove = FALSE)

new_tax_table_fam_matrix <- as.matrix(new_tax_table_fam)

psUF.family.merged.besthit@tax_table <- tax_table(new_tax_table_fam_matrix)

#Set sample order
sampleorder <- c("BNW-ST", "BNW-SI", "BNW-PC1", "BNW-PC2", "BNW-PC3", "BNW-PC4", "BNW-DC",
                 "SHNW-ST", "SHNW-PSI", "SHNW-DSI", "SHNW-PC1", "SHNW-PC2", "SHNW-PC3", "SHNW-DC1", "SHNW-DC3")

plot_heatmap(psUF.family.merged.besthit, taxa.label = "Phylum_Family", method = "RDA",
             sample.label = "sampleshort", low = "white", high = "red", na.value = "black",
             trans = log_trans(10), sample.order = sampleorder) +
  theme(axis.text.y = element_text(size=10, face = 'italic'),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, face = "bold"),
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white")) +
  facet_grid(.~species, scale = "free_x", space = "free_x")

ggsave(filename = "output/FigureSS4.png", width = 10, height = 20, dpi = 300)
ggsave(filename = "output/FigureSS4.pdf", width = 10, height = 20, dpi = 300)


###############
###HEAT MAP####
###############

#Make a colon-only table glommed to the genus level to send to Erin so she can run ANCOM2
#Note, I'm using the 'NArm = FALSE' option, such that ASVs without genus-level classifications
# are not discarded.

#How many ASVs in the colon only table?
prune_taxa(taxa_sums(ps.colon.rar) > 0, ps.colon.rar)

#Agglomerate classifications to the genus-level
ps.colon.rar.genus.NAs <- tax_glom(ps.colon.rar, taxrank = "Genus", NArm = FALSE)

#Sanity check to see how tax_glom agglomerates ASVs
Methano.test <- subset_taxa(ps.colon.rar, Genus=="Methanocorpusculum")
taxa_sums(Methano.test)
Methano.glommed <- subset_taxa(ps.colon.rar.genus.NAs, Genus=="Methanocorpusculum")
taxa_sums(Methano.glommed)
##Looks like it collapses all 5 ASVs of that rank to the ASV ID that has the highest abundance

#Export to a csv file for Erin
otu <- as.data.frame(ps.colon.rar.genus.NAs@otu_table)
tax <- as.data.frame(tax_table(ps.colon.rar.genus.NAs))
otu.and.tax <- cbind(tax, otu)
DT::datatable(otu.and.tax)
write.csv(otu.and.tax, file = "output/ps.colon.rarefied.genus.NAs.csv")

### See 'code/Wombat-GIT-ANCOM-Rcode.rmd' for the ANCOM2 analysis.

##Now, make a heatmap from Erin's ANCOM2 analysis

#Import Erin's ANCOM2 results
ANCOM.sampleposition <- read.csv("output/colon_rarefied_sampleposition_ANCOM_sums.csv")

#Pull out taxa that are significantly differentially abundant at atleast the 60% level
ANCOM.detected.0.7.taxa <- ANCOM.sampleposition %>% filter(detected_0.7=="TRUE") %>% pull(1)

#Merge replicate samples in genus globbed table (only want one cell per sample/taxon in figure)
ps.colon.rar.genus.NAs.merged <- merge_samples(ps.colon.rar.genus.NAs, "sampleshort2")
sam_data_fixed2 <- ps.colon.rar.genus.NAs@sam_data[- grep("B", ps.colon.rar.genus.NAs@sam_data$replicate),]
rownames(sam_data_fixed2) <- gsub("-A", "", rownames(sam_data_fixed2))
ps.colon.rar.genus.NAs.merged@sam_data = sam_data_fixed2
ps.colon.rar.genus.NAs.merged@otu_table <- phyloseq::t(ps.colon.rar.genus.NAs.merged@otu_table)

#Filter our phyloseq object to only keep significantly differentially abundant taxa
all.taxa = taxa_names(ps.colon.rar.genus.NAs.merged)
all.taxa.non.sig <- all.taxa[!(all.taxa %in% ANCOM.detected.0.7.taxa)]

#Sanity checking numbers
length(all.taxa)-length(all.taxa.non.sig)

#Math checks out
ps.colon.rar.genus.NAs.merged.ANCOM <- prune_taxa(ps.colon.rar.genus.NAs.merged,
                                                  taxa = ANCOM.detected.0.7.taxa)
#Format taxonomy to best hit
ps.colon.rar.genus.NAs.merged.ANCOM.besthit <- microbiomeutilities::format_to_besthit(
  ps.colon.rar.genus.NAs.merged.ANCOM, prefix = NULL)

#Create new rank which is phylum + best hit (for non-informative best hits -- e.g. 'uncultured')
new_tax_table <- data.frame(tax_table(ps.colon.rar.genus.NAs.merged.ANCOM.besthit)) %>%
  unite("Phylum_BestHit", c("Phylum", "Genus"), sep = " -- ", remove = FALSE)

new_tax_table_matrix <- as.matrix(new_tax_table)

ps.colon.rar.genus.NAs.merged.ANCOM.besthit@tax_table <- tax_table(new_tax_table_matrix)

#Add phylogenetic tree
ps.colon.rar.genus.NAs.merged.ANCOM.besthit@phy_tree = ps.colon.rar.genus.NAs.merged.ANCOM@phy_tree

#Set sample order
sampleorder <- c("BNW-PC1", "BNW-PC2", "BNW-PC3", "BNW-PC4", "BNW-DC",
                 "SHNW-PC1", "SHNW-PC2", "SHNW-PC3", "SHNW-DC1", "SHNW-DC3")

plot_heatmap(ps.colon.rar.genus.NAs.merged.ANCOM.besthit, taxa.label = "Phylum_BestHit", method = "RDA",
             sample.label = "sampleshort", low = "white", high = "red", na.value = "black",
             trans = log_trans(10), sample.order = sampleorder) +
  theme(axis.text.y = element_text(size=10, face = 'italic'),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, face = "bold"),
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white")) +
  facet_grid(.~species, scale = "free_x", space = "free_x")

ggsave(filename = "output/Figure4.png", width = 10, height = 10, dpi = 300)
ggsave(filename = "output/Figure4.svg", width = 10, height = 10, dpi = 300)
ggsave(filename = "output/Figure4.pdf", width = 10, height = 10, dpi = 300)


#####################
########Venn#########
#####################

## Figure 6 and SI figure 3 & 5
## Here, we're asking the question: "How many ASVs are shared between PC1 and the last
## distal colon sample for each wombat?"

#Pull out the first proximal colon and last distal colon sample for each species
ps.colon.rar.BNW <- subset_samples(ps.colon.rar, species == "Bare-nosed")
ps.colon.rar.SHNW <- subset_samples(ps.colon.rar, species == "Hairy-nosed")

ps.PC.BNW <- subset_samples(ps.colon.rar.BNW, sampleshort == "PC1")
ps.PC.BNW.table <- otu_table(ps.PC.BNW)

ps.DC.BNW <- subset_samples(ps.colon.rar.BNW, sampleshort == "DC")
ps.DC.BNW.table <- otu_table(ps.DC.BNW)

ps.PC.SHNW <- subset_samples(ps.colon.rar.SHNW, sampleshort == "PC1")
ps.PC.SHNW.table <- otu_table(ps.PC.SHNW)

ps.DC.SHNW <- subset_samples(ps.colon.rar.SHNW, sampleshort == "DC3")
ps.DC.SHNW.table <- otu_table(ps.DC.SHNW)

#Merge tables
ps.colon.all <- merge_phyloseq(ps.PC.BNW, ps.DC.BNW, ps.PC.SHNW, ps.DC.SHNW)

#Remove ASVs with counts of 0.
ps.colon.all <- prune_taxa(taxa_sums(ps.colon.all) > 0, ps.colon.all)

##I've set a detection threshold of >2 read counts per sample for this.
## I.e., an ASV needs at least 3 reads assigned to a sample type to be considered present.
## My rationale for this is to remove noise -- if it's only detected with 1 read, is it really there?

#For each ASV (row), if abundance > 2, print ASV (rowname) to a vector
venn.pc1.bnw.0 <- rownames(ps.PC.BNW.table[ apply(ps.PC.BNW.table, MARGIN = 1,
                                                  function(x) any(x > 2))])

venn.dc1.bnw.0 <- rownames(ps.DC.BNW.table[ apply(ps.DC.BNW.table, MARGIN = 1,
                                                  function(x) any(x > 2))])

venn.pc1.shnw.0 <- rownames(ps.PC.SHNW.table[ apply(ps.PC.SHNW.table, MARGIN = 1,
                                                    function(x) any(x > 2))])

venn.dc1.shnw.0 <- rownames(ps.DC.SHNW.table[ apply(ps.DC.SHNW.table, MARGIN = 1,
                                                    function(x) any(x > 2))])

#Create lists for each species and both species
bnw <- list(PC1 = venn.pc1.bnw.0, DC = venn.dc1.bnw.0)

shnw <- list(PC1 = venn.pc1.shnw.0, DC = venn.dc1.shnw.0)

bothspecies <- list(BNW.DC = venn.dc1.bnw.0, BNW.PC1 = venn.pc1.bnw.0,
                    SHNW.PC1 = venn.pc1.shnw.0, SHNW.DC = venn.dc1.shnw.0)

#Plot em using venn.diagram (figure 5)
venn.diagram(bnw, cex = 1.8, cat.cex = 2.5, print.mode = c("raw","percent"), fill = c("#f768a1", "#ae017e"), inverted = TRUE,
             imagetype = "tiff", filename = "output/Figure6A.tiff", cat.pos = c(0,0))

venn.diagram(shnw, cex = 1.8, cat.cex = 2.5, print.mode = c("raw","percent"), fill = c("#f768a1", "#ae017e"), inverted = TRUE,
             imagetype = "tiff", filename = "output/Figure6B.tiff", cat.pos = c(0,0))


##Same again, but at genus level
ps.colon.rar.genus <- tax_glom(ps.colon.rar, taxrank = "Genus", NArm = FALSE)


#Pull out the first proximal colon and last distal colon sample for each species
ps.colon.rar.BNW <- subset_samples(ps.colon.rar.genus, species == "Bare-nosed")
ps.colon.rar.SHNW <- subset_samples(ps.colon.rar.genus, species == "Hairy-nosed")

ps.PC.BNW <- subset_samples(ps.colon.rar.BNW, sampleshort == "PC1")
ps.PC.BNW.table <- otu_table(ps.PC.BNW)

ps.DC.BNW <- subset_samples(ps.colon.rar.BNW, sampleshort == "DC")
ps.DC.BNW.table <- otu_table(ps.DC.BNW)

ps.PC.SHNW <- subset_samples(ps.colon.rar.SHNW, sampleshort == "PC1")
ps.PC.SHNW.table <- otu_table(ps.PC.SHNW)

ps.DC.SHNW <- subset_samples(ps.colon.rar.SHNW, sampleshort == "DC3")
ps.DC.SHNW.table <- otu_table(ps.DC.SHNW)

#Merge tables
ps.colon.all <- merge_phyloseq(ps.PC.BNW, ps.DC.BNW, ps.PC.SHNW, ps.DC.SHNW)

#Remove ASVs with counts of 0.
ps.colon.all <- prune_taxa(taxa_sums(ps.colon.all) > 0, ps.colon.all)

##I've set a detection threshold of >2 read counts per sample for this.
## I.e., an ASV needs at least 3 reads assigned to a sample type to be considered present.
## My rationale for this is to remove noise -- if it's only detected with 1 read, is it really there?

#For each ASV (row), if abundance > 2, print ASV (rowname) to a vector
venn.pc1.bnw.0 <- rownames(ps.PC.BNW.table[ apply(ps.PC.BNW.table, MARGIN = 1,
                                                  function(x) any(x > 2))])

venn.dc1.bnw.0 <- rownames(ps.DC.BNW.table[ apply(ps.DC.BNW.table, MARGIN = 1,
                                                  function(x) any(x > 2))])

venn.pc1.shnw.0 <- rownames(ps.PC.SHNW.table[ apply(ps.PC.SHNW.table, MARGIN = 1,
                                                    function(x) any(x > 2))])

venn.dc1.shnw.0 <- rownames(ps.DC.SHNW.table[ apply(ps.DC.SHNW.table, MARGIN = 1,
                                                    function(x) any(x > 2))])

#Create lists for each species and both species
bnw <- list(PC1 = venn.pc1.bnw.0, DC = venn.dc1.bnw.0)

shnw <- list(PC1 = venn.pc1.shnw.0, DC = venn.dc1.shnw.0)

bothspecies <- list(BNW.DC = venn.dc1.bnw.0, BNW.PC1 = venn.pc1.bnw.0,
                    SHNW.PC1 = venn.pc1.shnw.0, SHNW.DC = venn.dc1.shnw.0)

#Plot em using venn.diagram (figure 5)
venn.diagram(bnw, cex = 1.4, cat.cex = 1.7, print.mode = c("raw","percent"), fill = c("#f768a1", "#ae017e"), inverted = FALSE,
             imagetype = "png", filename = "output/FigureSS5A_GENUS.png", cat.pos = c(320, 30))

venn.diagram(shnw, cex = 1.4, cat.cex = 1.7, print.mode = c("raw","percent"), fill = c("#f768a1", "#ae017e"), inverted = TRUE,
             imagetype = "png", filename = "output/FigureSS5B_GENUS.png", cat.pos = c(30,320))



#Alternative venn plot using ggVennDiagram https://github.com/gaospecial/ggVennDiagram
#Latest github version lets you pull out list of features in each region
#items <- get_region_items(x)
#openxlsx::write.xlsx(get_region_items(x), filename))

#Both species (SI figure 3)
ggVennDiagram(bothspecies, label_alpha = 0.5, label = "count") +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)
        ) +
  labs(fill = "ASVs")

ggsave(filename = "output/SI_Figure_3.png", width = 10, height = 10, dpi = 300)


##Now, let's pull out stats about the venn diagrams (i.e. what % of ASVs and rel ab.)
#Pull out ASVs belonging to specific regions of the cross-species Venn diagram
bothspecies.asv.lists <- get_region_items(bothspecies)
bothspecies.bothregions.shared <- bothspecies.asv.lists$ABCD
bothspecies.pc.shared <- bothspecies.asv.lists$BC
bothspecies.dc.shared <- bothspecies.asv.lists$AD

#BNW only Venn diagram
bnw.ASV.lists <- get_region_items(bnw)
bnw.pc.only <- bnw.ASV.lists$A
bnw.dc.only <- bnw.ASV.lists$B
bnw.shared <- bnw.ASV.lists$AB

#SHNW only Venn diagram
shnw.ASV.lists <- get_region_items(shnw)
shnw.pc.only <- shnw.ASV.lists$A
shnw.dc.only <- shnw.ASV.lists$B
shnw.shared <- shnw.ASV.lists$AB

##What's the relative abundance of site-specific ASVs in each region?
#Calculated as: number of reads for site-specific ASVs / number of reads for all ASVs at site

#SHNW PC
percent(
  mean(sample_sums(prune_taxa(shnw.pc.only, ps.PC.SHNW.table)))/
  mean(sample_sums(ps.PC.SHNW.table))
)
#SHNW DC
percent(
  mean(sample_sums(prune_taxa(shnw.dc.only, ps.DC.SHNW.table)))/
  mean(sample_sums(ps.DC.SHNW.table))
)
#BNW PC
percent(
  mean(sample_sums(prune_taxa(bnw.pc.only, ps.PC.BNW.table)))/
  mean(sample_sums(ps.PC.BNW.table))
)
#BNW DC
percent(
  mean(sample_sums(prune_taxa(bnw.dc.only, ps.DC.BNW.table)))/
  mean(sample_sums(ps.DC.BNW.table))
)
