# outlier removal via MAD

# library
library(data.table)
library(e1071)
library(reshape2)
library(ggplot2)

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"
cutoff.all <- c(5, 10, 25, 50, 100)

# ---------------------------------------------------------------------------- #
# mkdir
dir.save.up <- paste0("RESULT/999.003-Diagnosis_", ref)
dir.create(dir.save.up, recursive = TRUE)

# file to load
filename.data <- paste0("RESULT/1.3-GeneFiltering/ExpressionData_Filtered_",
                        ver, "_", ref, ".csv")

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


# ---------------------------------------------------------------------------- #
mad.max <- abs(apply(MAD.matrix, 2, max))
genes.mkfig <- names(mad.max)[10 < mad.max]
for (i in 1:length(genes.mkfig)) {
   gene.id <- genes.mkfig[i]
   expr.vec <- ExprMat[, gene.id]
   mad.vec <- MAD.matrix[, gene.id]
   df.fig <- data.frame("rlog2" = expr.vec,
                        "Filter" = "before_outlier_removal",
                        row.names = NULL)
   for ( k in 1:length(cutoff.all) ) {
      cutoff <- cutoff.all[k]
      df.k <- data.frame("rlog2" = expr.vec[mad.vec <= cutoff],
                         "Filter" = paste0("MAD_", cutoff),
                         row.names = NULL)
      df.fig <- rbind(df.fig, df.k)
   }
   df.fig$Filter <- factor(df.fig$Filter, levels = c("before_outlier_removal",
                                                     "MAD_100", "MAD_50", "MAD_25",
                                                     "MAD_10", "MAD_5"))
   p <- ggplot(df.fig, aes(x = rlog2))
   p <- p + geom_histogram()
   p <- p + facet_wrap(~ Filter, ncol = 1)
   p <- p + ggtitle(gene.id)
   ggsave(filename = paste0(dir.save.up, "/Fig-", gene.id, ".png"), p,
          width = 6, height = 8)
}