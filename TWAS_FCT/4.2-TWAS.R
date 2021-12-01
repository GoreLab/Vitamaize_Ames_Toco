# source
library(data.table)
library(rrBLUP)

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"

# mkdir
dir.create("RESULT/4.2-TWAS", recursive = T)

# file I/O
file.PeerResid <- paste0("RESULT/3.1-OutlierRemoval/PeerResiduals_RmOut_", ver, "_", ref, ".txt")
filename.summary <- paste0("RESULT/4.1-TWAS/TwasSummary_", ver, "_", ref, ".txt")

# file I; depends on reference genome
if ( ref == "B73" ) { file.gff <- "RAWDATA/Annotation/Zea_mays.B73_RefGen_v4.59_anno.csv" }
if ( ref == "PH207" ) { file.gff <- "RAWDATA/Annotation/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0.v03eksc_anno_merged.csv" } 

# load data
PeerResid <- fread(file.PeerResid, data.table = F)
GeneInfo <- fread(file.gff, data.table = F)
KinMat <- read.csv("RAWDATA/Kinship_NEW/GAPIT.Kin.VanRaden.csv", header = FALSE, stringsAsFactors = F)
Pheno.All <- read.delim("RAWDATA/Phenotype/ames_vita_blue_trans_v6version_locked_24_traits.txt")

# Peer residual matrix
PeerResidMat <- as.matrix(PeerResid[, -1])
myFun.01 <- function(VEC) { MIN <- min(VEC, na.rm = T); MAX <- max(VEC, na.rm = T); VEC.NEW <- (VEC - MIN) / (MAX - MIN); return(VEC.NEW)}
PeerResidMat.sc <- apply(X = PeerResidMat, 2, FUN = myFun.01)
rownames(PeerResidMat.sc) <- gsub("_", "", as.character(PeerResid$ID))

# chr & pos
m <- match(colnames(PeerResidMat.sc), GeneInfo$ID)
chr <- GeneInfo$chr[m]
pos <- GeneInfo$start[m]
snp <- GeneInfo$ID[m]

# kinship matrix
K <- as.matrix(KinMat[, -1])
rownames(K) <- colnames(K) <- KinMat$V1

# LOOP FOR ALL TRAITS
trait.all <- colnames(Pheno.All)[2:ncol(Pheno.All)]
for ( i in 1:length(trait.all) ) {
   # trait
   trait <- trait.all[i]
   
   # print 
   print(paste0("Start TWAS for ", trait, ": ", i, "-th trait out of ", length(trait.all), " traits."))
   
   # filename to save result
   file.TwasRes <- paste0("RESULT/4.2-TWAS/TwasResult_", ver, "_", ref, "_", trait, ".csv")
   
   # target phenotype data
   Pheno <- Pheno.All[, c("Taxa", trait)]
   
   # use intersect accessions
   common.accessions <- intersect(rownames(PeerResidMat.sc), Pheno$Taxa); length(common.accessions)
   my.Pheno <- Pheno[match(common.accessions, Pheno$Taxa), ]
   my.Geno <- data.frame("SNP" = snp, "Chr" = chr, "Pos" = pos, 
                         t(PeerResidMat.sc[match(common.accessions, rownames(PeerResidMat.sc)), ]),
                         row.names = 1:length(snp))
   my.K <- K[match(common.accessions, rownames(K)), match(common.accessions, rownames(K))]
   
   # Run GWAS
   b <- Sys.time()
   res.GWAS <- GWAS(pheno = my.Pheno, geno = my.Geno, K = my.K, n.PC = 0,
                    min.MAF = -Inf, P3D = FALSE, plot = FALSE, n.core = 40)
   colnames(res.GWAS)[4] <- "neg.log.P"
   colnames(res.GWAS)[1] <- "Gene"
   a <- Sys.time()
   print(a-b)
   
   # save result
   fwrite(res.GWAS, file = file.TwasRes)
}

# write summary
text.a <- paste0("TWAS: use ", ver, " ", ref, " dataset")
text.b <- paste0("Phenotype data has ", nrow(Pheno.All), " accessions")
text.c <- paste0("Peer residual matrix has ", nrow(PeerResidMat.sc), " accessions and ", ncol(PeerResidMat.sc), " genes")
text.d <- paste0("Number of overlapped accessions =  ", length(common.accessions))
SummaryText <- c(text.a, text.b, text.c, text.d)
write(x = SummaryText, file = filename.summary, sep = "/n")



