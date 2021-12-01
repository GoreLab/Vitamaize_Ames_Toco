# outlier removal via MAD

# library
library(data.table)

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"
cutoff <- as.numeric(args[3]) # numeric value for the MAD threshold

# mkdir
dir.save <- "RESULT/1.4-OutlierRemoval_and_Imputation"
dir.create(dir.save)

# file to load
filename.data <- paste0("RESULT/1.3-GeneFiltering/ExpressionData_Filtered_", ver, "_", ref, ".csv")

# load
ExprDat <- fread(file = filename.data, data.table = F)

# data frame to matrix
ExprMat <- as.matrix(ExprDat[, -1])
rownames(ExprMat) <- ExprDat[, 1]

# MAD calculation
MAD.matrix <- matrix(NA, nr = nrow(ExprMat), nc = ncol(ExprMat),
                     dimnames = dimnames(ExprMat))
for ( i in 1:ncol(ExprMat) ) {
   x <- ExprMat[, i]
   thres.i <- mad(x, constant = 1)
   mad.i <- abs(x - median(x)) / thres.i
   MAD.matrix[, i] <- mad.i
}

# save MAD matrix
MAD.Dat <- data.frame("Sample_ID" = ExprDat[, 1], MAD.matrix)
fwrite(MAD.Dat, file = paste0(dir.save, "/MAD_matrix_", ver, "_", ref, ".csv"))

# remove observations when MAD > cutoff
bool.mat <- cutoff < MAD.matrix
ExprMat.RmOut <- ExprMat
ExprMat.RmOut[bool.mat] <- NA
print(sum(bool.mat))

# write the number of outliers
n.rm <- apply(bool.mat, 2, sum)
tf <- nrow(ExprMat.RmOut) * 0.1 < n.rm
tab.n.rm <- table(n.rm)
df.save.01 <- data.frame("GeneID" = names(n.rm),
                         "n.outlier" = n.rm,
                         "flag" = c("keep", "remove")[as.numeric(tf)+1],
                         row.names = NULL)
fwrite(df.save.01, file = paste0(dir.save, "/number_of_outliers_", ver, "_", ref, ".csv"))

# write summary table of the outlier removal
df.save.02 <- data.frame(tab.n.rm)
colnames(df.save.02) <- c("Number_of_removed_outliers", "Freq")
fwrite(df.save.02, file = paste0(dir.save, "/summary_table_", ver, "_", ref, ".csv"))

# when more than than 10% of the samples were identified as outliers, remove the gene
tf <- nrow(ExprMat.RmOut) * 0.1 < n.rm
genes.remove <- colnames(ExprMat.RmOut)[tf]
ExprMat.RmOut.RmGene <- ExprMat.RmOut[, !tf]

# imputation by using median
tf <- apply(is.na(ExprMat.RmOut.RmGene), 2, sum) != 0
genes.with.na <- colnames(ExprMat.RmOut.RmGene)[tf]
ExprMat.Imp <- ExprMat.RmOut.RmGene
for ( j in 1:length(genes.with.na) ) {
   gene.j <- genes.with.na[j]
   expr.vec <- ExprMat.RmOut.RmGene[, gene.j]
   med.j <- median(expr.vec, na.rm = T)
   ExprMat.Imp[is.na(expr.vec), gene.j] <- med.j
}

# write the outlier-removed & median-imputed matrix
ExprDat.Imp <- data.frame("Sample_ID" = ExprDat[, 1], ExprMat.Imp)
fwrite(ExprDat.Imp, file = paste0(dir.save, "/ExpressionData_For_BLUE_", ver, "_", ref, ".csv"))
