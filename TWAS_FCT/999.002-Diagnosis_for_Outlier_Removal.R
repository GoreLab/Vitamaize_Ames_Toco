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
dir.save.up <- paste0("RESULT/999.002-Diagnosis_", ref)
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

# make summary
n.all <- prod(dim(MAD.matrix))
df.summary <- data.frame("MAD_threshold" = cutoff.all,
                         "n.removed.obs" = sapply(cutoff.all, 
                                                  function(x){sum(MAD.matrix > x)}))
df.summary$prcnt.removed.obs <- df.summary$n.removed.obs / n.all
write.csv(df.summary, file = paste0(dir.save.up, "/summary.csv"), row.names = F)

# ---------------------------------------------------------------------------- #
# for all cutoff, run the diagnosis
for ( c in 1:length(cutoff.all) ) {
   # cutoff
   cutoff <- cutoff.all[c]
   
   # mkdir
   dir.save <- paste0(dir.save.up, "/MAD_", cutoff)
   dir.create(dir.save, recursive = T)
   
   # Remove extreme values
   bool.mat <- cutoff < MAD.matrix
   ExprMat.RmOut <- ExprMat
   ExprMat.RmOut[bool.mat] <- NA
   
   # overall distribution
   file <- paste0(dir.save, "/Fig01.01-Hist_all_before.png")
   png(file,  width = 600, height = 500)
   hist(ExprMat, breaks = 100, 
        main = "Histogram BEFORE outlier removal",
        xlab = "rlog2 value")
   dev.off()
   
   file <- paste0(dir.save, "/Fig01.02-Hist_all_after.png")
   png(file,  width = 600, height = 500)
   hist(ExprMat.RmOut, breaks = 100, 
        main = "Histogram AFTER outlier removal",
        xlab = "rlog2 value")
   dev.off()
   
   # number of removed observations for each gene
   n.rm.per.gene <- apply(bool.mat, 2, sum)
   file <- paste0(dir.save, "/Fig02-Hist_of_the_number_of_removed_genes.png")
   png(file, width = 600, height = 500)
   hist(n.rm.per.gene[n.rm.per.gene != 0], 
        xlab = "Number of removed observations for each gene",
        main = paste0("Histogram of ", sum(n.rm.per.gene != 0), " genes"),
        breaks = 100)
   dev.off()
   
   # make figures
   genes.with.outliers <- names(n.rm.per.gene)[10 < n.rm.per.gene]
   ExprMat.with.outl <- ExprMat[, genes.with.outliers]
   ExprMat.RmOut.with.outl <- ExprMat.RmOut[, genes.with.outliers]
   o <- order(apply(is.na(ExprMat.RmOut.with.outl), 2, sum), decreasing = T)
   ExprMat.with.outl.o <- ExprMat.with.outl[, o]
   ExprMat.RmOut.with.outl.o <- ExprMat.RmOut.with.outl[, o]
   for ( k in 1:ncol(ExprMat.with.outl.o) ) {
      num <- formatC(k, width = 4, flag = "0")
      gene.id <- colnames(ExprMat.with.outl.o)[k]
      n.rm <- sum(is.na(ExprMat.RmOut.with.outl.o[, gene.id]))
      df.fig.01 <- melt(ExprMat.with.outl.o[, gene.id])
      df.fig.02 <- melt(ExprMat.RmOut.with.outl.o[, gene.id])
      df.fig.01$Data <- "all rlog2 values"
      df.fig.02$Data <- "after the outlier removal"
      df.fig <- rbind(df.fig.01, df.fig.02)
      p <- ggplot(df.fig, aes(x = value))
      p <- p + geom_histogram(bins = 100)
      p <- p + facet_wrap(~ Data, ncol = 1)
      p <- p + ggtitle(paste0(gene.id, ": we remove ", n.rm, " observations"))
      file <- paste0(dir.save, "/Fig03.", num, "-Hist_of_", gene.id, ".png")
      ggsave(filename = file, p, width = 6, height = 4)
   }
   
   # write the number of observations removed per gene
   df.save <- data.frame("Gene_ID" = names(n.rm.per.gene),
                         "n.obs.removed" = n.rm.per.gene)
   write.csv(df.save, file = paste0(dir.save, "/Table_number_of_observations_removed.csv"), row.names = F)
}













