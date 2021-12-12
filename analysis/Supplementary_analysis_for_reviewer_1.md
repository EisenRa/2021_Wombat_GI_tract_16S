# Supplementary analyses for reviewer 1
#### Raphael Eisenhofer, 12_2021

# Analysis 1:
**Reviewer 1 comment**: "My main concern is related to the preliminary data manipulation prior to the analyses. Specifically, the exclusion of ASVs found in a single sample. While these ASVs will not give you any information about sharing among samples they may provide important information on the individuality of samples and on alpha diversity. The removal of singltons (ASVs for which only a single read is returned) is common and reasonable as they often represent sequencing error but the removal of a sequence because it is only seen in a single sample could be misleading. If my interpretation of how the analyses were performed is correct then this will artificially reduce the variation between samples, making them appear more homogeneous. Reducing variation in this way could make differences between the GI sites appear more convincing. Therefore, I think it is important for all ASVs to be included in the analyses or at the very least for the authors to show (in an online supplement) that it does not change the findings."

**Response**: See (Supp_analyses_for_reviewers.R)[https://github.com/EisenRa/2021_Wombat_GI_tract_16S/blob/master/code/Supp_Analysis_Reviewer1.R] in the 'code' folder for the code used for these. The number of ASVs found in only 1 sample (i.e. replicate, as we have two samples per region -- e.g. PC1a & PC1b) is **105/2419 = ~4.3%**

These ASVs account for only **0.16%** of the average relative abundance for the dataset. 

### Figure 1 (from original submission -- now figure 2 with new figure addition as recommended by reviewer 2)
#### Bare-nosed wombat ASV richness
![Figure 1a](https://github.com/EisenRa/2021_Wombat_GI_tract_16S/blob/master/analysis/figures/FigureSS1A.png)
#### Southern hairy-nosed wombat ASV richness
![Figure 1b](https://github.com/EisenRa/2021_Wombat_GI_tract_16S/blob/master/analysis/figures/FigureSS1B.png)

### Figure 2
![Figure 2](https://github.com/EisenRa/2021_Wombat_GI_tract_16S/blob/master/analysis/figures/FigureSS2.png)
#### A = unweighted unifrac, B = weighted unifrac

### Supplementary figure 1 (colon-only ordination)
![Figure SI1](https://github.com/EisenRa/2021_Wombat_GI_tract_16S/blob/master/analysis/figures/SI_FigureSS1.png)
#### A = unweighted unifrac, B = weighted unifrac

Note, figure 3 would not be influnced by this analysis, as it displayed only the top 20 most abundant families. Likewise for figure 4, as these ultra-low abundance features would be filtered by the ANCOM-II analysis. Hence they are not displayed here.

### Figure 5
#### Bare-nosed wombat
![Figure 5a](https://github.com/EisenRa/2021_Wombat_GI_tract_16S/blob/master/analysis/figures/FigureSS5A.png)
#### Southern hairy-nosed wombat
![Figure 5b](https://github.com/EisenRa/2021_Wombat_GI_tract_16S/blob/master/analysis/figures/FigureSS5B.png)

### Conclusions:
As expected, this did not have a discenable impact on the figures/interpretations of the initial analyses. 


# Analysis 2:
**Reviewer 1 comment**:
"Figure 3: I would like to see a larger number of families represented in this plot as for the DC in the SHNW a large proportion of the microbiome is assigned to “other” making it difficult to see the detail of this community."

Unfortunately, it is too difficult to display 161 colours in a stacked bar chart. As a compromise, I have created a heatmap containing all the families. Note that blank family names correspond to ASVs that did not have a family classification (i.e. they were classified higher up).

![Figure 4](https://github.com/EisenRa/2021_Wombat_GI_tract_16S/blob/master/analysis/figures/FigureSS4.png)


# Analysis 3:
**Reviewer 1 comment**:
"L232-240: I would be really interested to see what the overlap is between the sites if you repeat this analysis at the genera or family level. Given the known issues with intragenomic heterogeneity in 16S genes and also in assigning ASVs to lower taxonomic levels the differences seen in composition may not be biologically meaningful."

Here is the venn results at the genus-level:

#### BNW
![Figure 5a_gen](https://github.com/EisenRa/2021_Wombat_GI_tract_16S/blob/master/analysis/figures/FigureSS5A_GENUS.png)

#### SHNW
![Figure 5b_gen](https://github.com/EisenRa/2021_Wombat_GI_tract_16S/blob/master/analysis/figures/FigureSS5B_GENUS.png)

Now, the relative abundance of site-specific (i.e. DC only) genera are:
BNW DC = 5%
BNW PC = 20%
SHNW DC = 3%
SHNW PC = 12%

