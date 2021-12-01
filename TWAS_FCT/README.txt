
This folder contains all scripts for transcriptome wide association study (TWAS) and Fisher's combined test (FCT) analysis, including the RNAseq data processing such as gene filtering. These scripts cover the “Expression data analysis”, “Transcriptome-wide association study” and “Fisher’s combined tests of GWAS and TWAS” sections in the manuscript.

0.1-MergeGenes_B73.R merges overlapped gene regions from the GFF file. Result of this script will be used only for the FCT analysis.

1.1-MakeDataset.R merges multiple files with a few manual corrections of the wrong record in raw data. This code does not do any calculation. Output files will be used in next steps.

1.2-AttachGDD.R calculates GDD based on the pollination and harvest date by using weather (minimum and maximum temperature) data.

1.3-GeneFiltering.R removes genes if less than half of the samples showed zero rlog2 value.

1.4-OutlierRemoval_and_Imputation.R removes outliers based on median absolute deviation and impute those values with the median expression profile. If more than 10% of the samples were detected as outliers, that gene was removed from the data in this script.

1.5.1-BLUE_each.R calculates BLUE values based on the mixed linear model. This R code was executed by using the bash file 1.5-RunMultipleR_v1.1_B73.sh for all genes.

1.5.2-BLUE_merge.R merges all output files of BLUE values (previous R script results in one file for one gene).

1.5.3-BLUE_RmErr.R removes genes for which BLUE calculation did not converge if there is. BLUE calculation converged for all genes and therefore this code does nothing in our pipeline. However, we kept using this code so that we can diagnose the BLUE calculation.  

1.5.4-BLUE_SampleFiltering.R excludes 104 lines classified as sweet corn, popcorn, or with an endosperm mutation. This also removes 15 lines not analyzed in GWAS.

2.1-Peer_Use25Fact.R executes PEER with 25 factors.

2.2-MakeScreePlot.R draws the diagnostic plot of the factor relevance. The figure was visually assessed, and the “elbow” was identified.


2.3-Peer_UseOptFact.R executes PEER with the optimal number of factors.

3.1-OutlierRemoval.R calculates Studentized deleted residuals by fitting a simple linear model that has only the grand mean and removes outliers by using a Bonferroni adjusted significance threshold of α = 0.05.

4.1-BIC_for_TWAS.R calculates BIC to select the optimal model for TWAS.

4.2-TWAS.R executes TWAS by using the optimal model selected in the previous step.

5.1-CombinedTest.R executes Fisher’s combined test. Prior to this analysis

999.999.001-AnalyseVte7_01.R executes the same data processing pipeline on the CPM data of vte7, from gene filtering to BLUE calculation.

999.999.001-AnalyseVte7_02.R executes the same data processing pipeline on the CPM data of vte7, from sample filtering (before PEER) to TWAS.
  
All the scripts were exectued as written in ExecuteScripts.sh file.