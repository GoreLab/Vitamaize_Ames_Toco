# count numbers etc

# package
library(data.table)

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"

###################### 
ver <- "v1.1" # "v1" or "v1.1"
ref <- "PH207" # "B73" or "PH207"
######################


# gene filtering
f <- paste0("RESULT/1.3-GeneFiltering/Summary_", ver, "_", ref, ".csv")
summary.filt <- read.csv(f)
table(summary.filt$flag)

# outlier removal & imputation
f <- paste0("RESULT/1.4-OutlierRemoval_and_Imputation/number_of_outliers_", ver, "_", ref, ".csv")
summary.outl <- read.csv(f)
sum(summary.outl$n.outlier)
table(summary.outl$flag)
sum(summary.outl$n.outlier[summary.outl$flag == "keep"])
table(summary.outl$n.outlier == 0)
table((summary.outl$n.outlier == 0)[summary.outl$flag == "keep"])


# BLUE with all accessions (include checks)
f <- paste0("RESULT/1.5-BLUE/expression_BLUE_", ver, "_", ref, "_raw.csv")
BLUE.all <- fread(file = f, data.table = F)
BLUE.all.mat <- t(as.matrix(BLUE.all[, -1]))
colnames(BLUE.all.mat) <- BLUE.all[, 1]
BLUE.all.mat[1:5, 1:5]
dim(BLUE.all.mat)

# BLEU for the experimental & phenotyped 545 (in v1.1) accessions 
f <- paste0("RESULT/1.5-BLUE/expression_BLUE_final_", ver, "_", ref, ".csv")
BLUE.fin <- fread(file = f, data.table = F)
BLUE.fin.mat <- as.matrix(BLUE.fin[, -1])
rownames(BLUE.fin.mat) <- BLUE.fin[, 1]
BLUE.fin.mat[1:5, 1:5]





