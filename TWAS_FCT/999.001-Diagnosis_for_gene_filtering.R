# packages
library(gdata)
library(data.table)
library(ggplot2)

# mkdir
dir.save <- "RESULT/999.001-Diagnosis_for_gene_filtering"
dir.create(dir.save, recursive = T)

# params
args <- commandArgs(trailingOnly = T)
ver <- "v1.1" # this script is "v1.1" only!!
ref <- args[1] # "B73" or "PH207"

# load files
ExprDat <- fread(paste0("RESULT/1.1-MakeDataset/ExpressionData_v1.1_", ref, ".csv"),
                 data.table = F)
Info <- read.csv(paste0("RESULT/1.1-MakeDataset/SampleInfo_v1.1_", ref, ".csv"))
expr.mat.all <- as.matrix(ExprDat[, -1])
rownames(expr.mat.all) <- ExprDat[, 1]

# load file
KinData <- read.csv("RAWDATA/Kinship_NEW/GAPIT.Kin.VanRaden.csv", header = F)
K <- as.matrix(KinData[, -1])
rownames(K) <- colnames(K) <- KinData[, 1]

#
info.geno.vec <- gsub("_", "", Info$GENOTYPE)
tf <- info.geno.vec %in% rownames(K)
expr.mat.sub <- expr.mat.all[tf, ]
dim(expr.mat.sub)

# filter on 741 samples
pr.zero.UseAll <- apply(expr.mat.all == 0, 2, sum) / nrow(expr.mat.all)
df.n.rm <- data.frame("pr.zero" = seq(0, 0.9, by = 0.1),
                      "n.gene.rm" = NA)
for ( i in 1:nrow(df.n.rm) ) {
   cutoff <- df.n.rm$pr.zero[i]
   df.n.rm$n.gene.rm[i] <- sum(cutoff < pr.zero.UseAll)
}
df.n.rm$pr.zero.prcnt <- paste0(df.n.rm$pr.zero * 100, "%")
p <- ggplot(df.n.rm, aes(x = pr.zero.prcnt, y = n.gene.rm))
p <- p + geom_bar(stat = "identity")
p <- p + geom_text(aes(label = n.gene.rm), vjust = -1)
p <- p + xlab("threshold: % of zero samples")
p <- p + ylab("# of removed genes")
p <- p + ylim(c(0, 3500))
p <- p + ggtitle("Gene filtering on the all 741 samples")
ggsave(filename = paste0(dir.save, "/Fig01-Num.Rm.Genes.UseAll_", ref, ".png"),
       p, width = 6, height = 4)

# filter on 546 samples
pr.zero.UseSub <- apply(expr.mat.sub == 0, 2, sum) / nrow(expr.mat.sub)
df.n.rm <- data.frame("pr.zero" = seq(0, 0.9, by = 0.1),
                      "n.gene.rm" = NA)
for ( i in 1:nrow(df.n.rm) ) {
   cutoff <- df.n.rm$pr.zero[i]
   df.n.rm$n.gene.rm[i] <- sum(cutoff < pr.zero.UseSub)
}
df.n.rm$pr.zero.prcnt <- paste0(df.n.rm$pr.zero * 100, "%")
p <- ggplot(df.n.rm, aes(x = pr.zero.prcnt, y = n.gene.rm))
p <- p + geom_bar(stat = "identity")
p <- p + geom_text(aes(label = n.gene.rm), vjust = -1)
p <- p + xlab("threshold: % of zero samples")
p <- p + ylab("# of removed genes")
p <- p + ggtitle("Gene filtering on the all 546 samples")
p <- p + ylim(c(0, 3500))
ggsave(filename = paste0(dir.save, "/Fig02-Num.Rm.Genes.UseSub_", ref, ".png"),
       p, width = 6, height = 4)


# # ---------------------------------------------------------------------------- #
# genes.rm.UseAll <- names(pr.zero.UseAll)[0.5 < pr.zero.UseAll]
# genes.rm.UseSub <- names(pr.zero.UseSub)[0.5 < pr.zero.UseSub]
# length(intersect(genes.rm.UseAll, genes.rm.UseSub))
# setdiff(genes.rm.UseAll, genes.rm.UseSub)
# setdiff(genes.rm.UseSub, genes.rm.UseAll)








