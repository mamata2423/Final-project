# Final Project

## DOI
[![DOI](https://zenodo.org/badge/969865982.svg)](https://doi.org/10.5281/zenodo.15278537)

# Growth and Biofilm Formation of *Xylella fastidiosa* Almond Strains in Grape Sap Amended Media
This repository contains data, code, and results for the project analyzing the growth and biofilm formation of 11 different strains of *Xylella fastidiosa* isolated from Almonds under four different media. The primary goal of this study is to compare the phenotypic variation in growth and biofilm formation across strains and media, including media supplemented with grape sap, to gain insights into host adaptation and potential virulence of almond strains in grapevines.
## DataSets
### Raw Data: 
The folder contains Excel files with 96-well plate biofilm data measured on Day 7 and growth data from Day 0 to Day 7 in the Cytation 3 plate reader.
Biofilm.xlsx: Biofilm quantification raw data (OD600).
Day0.xlsx - Day7.xlsx: Growth raw data (OD600).
- [Raw Data](Data/RawData)
### Clean Data: 
The folder contains cleaned and processed .csv files of biofilm and growth data for statistical tests.
Biofilm.csv: Final cleaned and processed biofilm data for analysis.
GrowthData.csv: Final cleaned and processed growth data from Day 0 to Day 7 for analysis.
- [Clean Data](Data/CleanData)
## Media Types
The study used four types of media to investigate the effect of grape sap on different almond-derived strains of *Xylella fastidiosa*.
-	PD3: Artificial media
-	25SAP: 25% grape sap in PD3
-	50SAP: 50% grape sap in PD3
-	75SAP: 75% grape sap in PD3
## Experimental Summary
*	**Strains:** 11 almond-derived strains of *Xylella fastidiosa*.
*	**Experimental Format:**
    -	96-well plate assays
    -	7-day growth monitoring and biofilm quantification on the 7th day.
    -	The first column (Column 1) was control, and Columns 2 to 12 contained 11 different strains of *Xylella fastidiosa*.
    -	Each well was inoculated with respective 0.01 OD600 of the bacterial strains, except the control wells.
*	**Phenotypes Measured:**
    -	Bacterial growth over time every 24 hours for 7 days (Absorbance reading at OD600)
    -	Area Under Curve (AUC) from growth curves
    -	Biofilm formation using the crystal violet assay.
## Results
-	[Figures](Figures) contain manuscript-ready figures, which include Growth curves, Area under Growth curves, and biofilm plots for all the strains combined within a single plot based on different media. 
- [figure-gfm](FinalProject_files/figure-gfm) includes figures for each strain during the analysis.
-	[FinalProject.Rmd](FinalProject.Rmd) contains the complete code for the analysis.
-	[FinalProject.md](FinalProject.md) is a markdown file that is favored by GitHub.
-	In addition, [FinalProject.pdf](FinalProject.pdf) and [FinalProject.html](FinalProject.html) include all the codes with results of statistical analysis, figures with a table of contents, and headers as well.
## Statistical Analysis
-	**Biofilm Analysis:**
Raw OD600 data were cleaned by removing irrelevant rows, combining sheets in Excel files, subtracting control values, and making a long-format dataset. OD600 values were log-transformed to meet linear model assumptions. A linear mixed-effects model was fitted for each strain to assess the effect of different media on log-transformed biofilm values using a replicate plate as a random effect. Estimated marginal means (emmeans) were calculated to compare media effects within each strain, and compact letter display (cld) was used to assign significance groupings (p<0.05) among treatments in boxplots generated using the ggplot2 package.

-	**Growth Data Analysis:** 
Raw data were processed similarly to above in the biofilm data, and OD600 values were log-transformed to visualize relative growth more effectively, and mean and standard errors were calculated by strain, media, and day. Line plots were used to visualize the growth curves. Furthermore, the area under the growth curve was also calculated, and a linear mixed effects model was fitted for each strain with media as a fixed effect and replicate as a random effect. Barplots were created comparing mean AUC values across media for each strain.

## File Tree

```r
fs::dir_tree()
```

```
├── Data                         #Includes all data
│   ├── CleanData                #cleaned and processed data
│   │   ├── Biofilm.csv          #cleaned and processed biofilm data
│   │   └── GrowthData.csv       #cleaned and processed growth data
│   └── RawData                  #raw data obtained directly from the plate reader
│       ├── Biofilm.xlsx         #raw biofilm data
│       ├── Day0.xlsx            #raw data from day 0 to day 7
│       ├── Day1.xlsx
│       ├── Day2.xlsx
│       ├── Day3.xlsx
│       ├── Day4.xlsx
│       ├── Day5.xlsx
│       ├── Day6.xlsx
│       └── Day7.xlsx
├── Figures                     #manuscript ready figures 
│   ├── GrowthCurves.png        #growth curve figure for all strains combined
│   ├── combined_AUC_plot.png   #Area under curve figure for all strains combined
│   └── combined_biofilm_plot.png  #biofilm plots for all strains combined
├── FinalProject.Rmd            #R Markdown file for analysis
├── FinalProject.Rproj          #top level directory .Rproj file = working directory
├── FinalProject.html           #html rendered file for analysis
├── FinalProject.md             #GitHub formatted markdown
├── FinalProject.pdf            #pdf document rendered file for analysis
├── FinalProject_files  
│   └── figure-gfm              #includes all the intermediate unnamed figures generated during analysis
│       ├── unnamed-chunk-2-1.png
│       ├── unnamed-chunk-2-10.png
│       ├── unnamed-chunk-2-11.png
│       ├── unnamed-chunk-2-12.png
│       ├── unnamed-chunk-2-2.png
│       ├── unnamed-chunk-2-3.png
│       ├── unnamed-chunk-2-4.png
│       ├── unnamed-chunk-2-5.png
│       ├── unnamed-chunk-2-6.png
│       ├── unnamed-chunk-2-7.png
│       ├── unnamed-chunk-2-8.png
│       ├── unnamed-chunk-2-9.png
│       ├── unnamed-chunk-4-1.png
│       ├── unnamed-chunk-5-1.png
│       ├── unnamed-chunk-5-10.png
│       ├── unnamed-chunk-5-11.png
│       ├── unnamed-chunk-5-12.png
│       ├── unnamed-chunk-5-2.png
│       ├── unnamed-chunk-5-3.png
│       ├── unnamed-chunk-5-4.png
│       ├── unnamed-chunk-5-5.png
│       ├── unnamed-chunk-5-6.png
│       ├── unnamed-chunk-5-7.png
│       ├── unnamed-chunk-5-8.png
│       └── unnamed-chunk-5-9.png
└── README.md                     #top level directory README
```

