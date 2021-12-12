# Supplementary analyses for reviewer 1
#### Raphael Eisenhofer, 12_2021

# Analysis 1:
**Reviewer 1 comment**: "My main concern is related to the preliminary data manipulation prior to the analyses. Specifically, the exclusion of ASVs found in a single sample. While these ASVs will not give you any information about sharing among samples they may provide important information on the individuality of samples and on alpha diversity. The removal of singltons (ASVs for which only a single read is returned) is common and reasonable as they often represent sequencing error but the removal of a sequence because it is only seen in a single sample could be misleading. If my interpretation of how the analyses were performed is correct then this will artificially reduce the variation between samples, making them appear more homogeneous. Reducing variation in this way could make differences between the GI sites appear more convincing. Therefore, I think it is important for all ASVs to be included in the analyses or at the very least for the authors to show (in an online supplement) that it does not change the findings."

**Response**: See Supp_analyses_for_reviewers.R in the 'code' folder for the code used for these. The number of ASVs found in only 1 sample (i.e. replicate, as we have two samples per region -- e.g. PC1a & PC1b) is **105/2419 = ~4.3%**

These ASVs account for only **0.16%** of the average relative abundance for the dataset. 

### Figure 1 (from original submission -- now figure 2 with new figure addition as recommended by reviewer 2)
![Figure 1](figure/FigureSS1A.png)
