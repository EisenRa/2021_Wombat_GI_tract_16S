### Wombat GI tract 16S paper 2021
### Raphael Eisenhofer, December 2021
### Supplementary analysis for reviewer 1
### The reviewer wanted to see what would happen if we didn't filter ASVs found in less
### than 1 biological replicate.

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

psUF <- subset_samples(ps, sampletype != "negative control")

#Rename to get object to work with subsequent scripts
psUF.prevF <- psUF

##Sanity check to make sure we have not filtered object:
#Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(psUF.prevF),
                MARGIN = ifelse(taxa_are_rows(psUF.prevF), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
#Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(psUF.prevF),
                     tax_table(psUF.prevF))

count(prevdf, vars = Prevalence <= 1)

#What's the relative abundance of these ASVs?
OneReplicateASVs <- rownames(prevdf)[(prevdf$Prevalence <= 1)]
OneReplicateASVs_PS <- prune_taxa(OneReplicateASVs, psUF)

percent(mean(sample_sums(OneReplicateASVs_PS)/sample_sums(psUF.prevF)), accuracy = 0.04)
write.csv(percent(sample_sums(OneReplicateASVs_PS)/sample_sums(psUF.prevF), accuracy = 0.04), file = "OneRepASV.csv")


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
colours <- c("red", "orange", "yellow", "forestgreen")

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
  geom_boxplot() +
  #Jitter, size, colour
  geom_jitter(position=position_dodge2(0.3), size=3.5, aes(colour=sampletype)) +
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
ggsave(filename = "output/FigureSS1A.png", width = 20, height = 6)


svg(filename = "output/FigureSS1A.png", width = 20, height = 6)

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
  geom_boxplot() +
  #Jitter, size, colour
  geom_jitter(position=position_dodge2(0.3), size=3.5, aes(colour=sampletype)) +
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

ggsave(filename = "output/FigureSS1B.png", width = 20, height = 6)

svg(filename = "output/FigureSS1B.svg", width = 20, height = 6)

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

ggsave(filename = "output/FigureSS2.png", width = 10, height = 10, dpi = 300)
ggsave(filename = "output/FigureSS2.svg", width = 10, height = 10, dpi = 300)
ggsave(filename = "output/FigureSS2.pdf", width = 10, height = 10, dpi = 300)


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

ggsave(filename = "output/SI_FigureSS1.png", width = 10, height = 10, dpi = 300)
ggsave(filename = "output/SI_FigureSS1.svg", width = 10, height = 10, dpi = 300)
ggsave(filename = "output/SI_FigureSS1.pdf", width = 10, height = 10, dpi = 300)


###############
###HEAT MAP####
###############

## Collapse taxonomy to family level
## Print relative abundance values of the different families

psUF.family <- tax_glom(psUF, taxrank = "Family", NArm = FALSE)

psUF.family.merged <- merge_samples(psUF.family, "sampleshort2")


#Set sample order
sampleorder <- c("BNW-ST", "BNW-SI", "BNW-PC1", "BNW-PC2", "BNW-PC3", "BNW-PC4", "BNW-DC",
                 "SHNW-ST", "SHNW-PSI", "SHNW-DSI", "SHNW-PC1", "SHNW-PC2", "SHNW-PC3", "SHNW-DC1", "SHNW-DC3")

plot_heatmap(psUF.family.merged, taxa.label = "Family", method = "RDA",
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




#####################
########Venn#########
#####################

## Figure 5 and SI figure 3
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
bnw <- list(BNW.PC1 = venn.pc1.bnw.0, BNW.DC = venn.dc1.bnw.0)

shnw <- list(SHNW.PC1 = venn.pc1.shnw.0, SHNW.DC = venn.dc1.shnw.0)

bothspecies <- list(BNW.DC = venn.dc1.bnw.0, BNW.PC1 = venn.pc1.bnw.0,
                    SHNW.PC1 = venn.pc1.shnw.0, SHNW.DC = venn.dc1.shnw.0)

#Plot em using venn.diagram (figure 5)
venn.diagram(bnw, print.mode = c("raw","percent"), fill = c("orange", "red"), inverted = TRUE,
             imagetype = "png", filename = "output/FigureSS5A.png", cat.pos = c(0,0))

venn.diagram(shnw, print.mode = c("raw","percent"), fill = c("orange", "red"), inverted = TRUE,
             imagetype = "png", filename = "output/FigureSS5B.png", cat.pos = c(0,0))


## SAME, BUT WITHOUT 2> reads per ASV per site filter:

#For each ASV (row), if abundance > 2, print ASV (rowname) to a vector
venn.pc1.bnw.0 <- rownames(ps.PC.BNW.table[ apply(ps.PC.BNW.table, MARGIN = 1,
                                                  function(x) any(x > 0))])

venn.dc1.bnw.0 <- rownames(ps.DC.BNW.table[ apply(ps.DC.BNW.table, MARGIN = 1,
                                                  function(x) any(x > 0))])

venn.pc1.shnw.0 <- rownames(ps.PC.SHNW.table[ apply(ps.PC.SHNW.table, MARGIN = 1,
                                                    function(x) any(x > 0))])

venn.dc1.shnw.0 <- rownames(ps.DC.SHNW.table[ apply(ps.DC.SHNW.table, MARGIN = 1,
                                                    function(x) any(x > 0))])

#Create lists for each species and both species
bnw <- list(BNW.PC1 = venn.pc1.bnw.0, BNW.DC = venn.dc1.bnw.0)

shnw <- list(SHNW.PC1 = venn.pc1.shnw.0, SHNW.DC = venn.dc1.shnw.0)

bothspecies <- list(BNW.DC = venn.dc1.bnw.0, BNW.PC1 = venn.pc1.bnw.0,
                    SHNW.PC1 = venn.pc1.shnw.0, SHNW.DC = venn.dc1.shnw.0)

#Plot em using venn.diagram (figure 5)
venn.diagram(bnw, print.mode = c("raw","percent"), fill = c("orange", "red"), inverted = TRUE,
             imagetype = "tiff", filename = "output/FigureSS_no_site_filter_5A.tiff", cat.pos = c(0,0))

venn.diagram(shnw, print.mode = c("raw","percent"), fill = c("orange", "red"), inverted = TRUE,
             imagetype = "tiff", filename = "output/FigureSS5B_no_site_filter_.tiff", cat.pos = c(0,0))


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

ggsave(filename = "output/SI_Figure_SS3.png", width = 10, height = 10, dpi = 300)

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



#####################
########Venn#########
#####################

## This time, glommed at the genus level

#Glom to genus
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
bnw <- list(BNW.PC1 = venn.pc1.bnw.0, BNW.DC = venn.dc1.bnw.0)

shnw <- list(SHNW.PC1 = venn.pc1.shnw.0, SHNW.DC = venn.dc1.shnw.0)

bothspecies <- list(BNW.DC = venn.dc1.bnw.0, BNW.PC1 = venn.pc1.bnw.0,
                    SHNW.PC1 = venn.pc1.shnw.0, SHNW.DC = venn.dc1.shnw.0)

#Plot em using venn.diagram (figure 5)
venn.diagram(bnw, print.mode = c("raw","percent"), fill = c("orange", "red"), inverted = TRUE,
             imagetype = "png", filename = "output/FigureSS5A_GENUS.png", cat.pos = c(0,0))

venn.diagram(shnw, print.mode = c("raw","percent"), fill = c("orange", "red"), inverted = TRUE,
             imagetype = "png", filename = "output/FigureSS5B_GENUS.png", cat.pos = c(0,0))


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

