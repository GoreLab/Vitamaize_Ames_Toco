# package
library(data.table)
library(MASS)
library(gdata)
library(rrBLUP)

# library
dir.save <- "RESULT/999-VTE7"
dir.create(dir.save, recursive = T)

# ---------------------------------------------------------------------------- #
# load data & BLUE diagnosis
# ---------------------------------------------------------------------------- #
# load 
blue.dat <- read.csv("RESULT/999-VTE7/BLUE.csv")
expr.dat <- read.csv("RESULT/999-VTE7/ExprData.csv")
info <- read.csv("RESULT/999-VTE7/SampleInfo.csv")

# average expression
avg.expr <- tapply(X = expr.dat$CPM, INDEX = info$GENOTYPE, mean)

# experimental lines
expr.gen.id <- unique(info$GENOTYPE[info$CHECK != "Check"])
avg.expr.664 <- avg.expr[expr.gen.id]
blue.664 <- blue.dat$BLUE[match(expr.gen.id, blue.dat$Genotype)]
names(blue.664) <- blue.dat$Genotype[match(expr.gen.id, blue.dat$Genotype)]

# check -> OK
all.equal(names(avg.expr.664), names(blue.664))

# avg and BLUE
f <- paste0(dir.save, "/Fig004-CPM_and_BLUE.png")
png(filename = f, width = 600, height = 600)
lim <- range(c(avg.expr.664, blue.664))
plot(x = avg.expr.664, y = blue.664,
     xlim = lim, ylim = lim,
     xlab = "Average CPM", ylab = "BLUE", 
     main = "664 experimental lines",
     pch = 20)
abline(0, 1, lty = 2)
dev.off()



# ---------------------------------------------------------------------------- #
# Sample Filtering
# ---------------------------------------------------------------------------- #
# load accession info data
AccInfo <- read.xls("RAWDATA/Transcript/vitamaize_accession_removal_20191113.xlsx")
geno.names.info <- as.character(AccInfo$Accession.Number_ion)

# convert name
names(blue.664) <- gsub("_", "", names(blue.664))

# common ones etc
geno.names.common <- intersect(names(blue.664), geno.names.info)


# keep/remove flag -> accessions to be kept
m <- match(geno.names.common, geno.names.info)
tmp <- AccInfo[m, ]
tf.tmp.remove <- tmp$Removal_decision == "remove"
accessions.remove <- as.character(tmp$Accession.Number_ion)[tf.tmp.remove]
accessions.remove <- c(accessions.remove, "PI644099") # manually add PI644099
accessions.keep <- setdiff(names(blue.664), accessions.remove)

# further filter accessions: if there is no phenotype, remove them
pheno.data <- read.table("RAWDATA/Phenotype/ames_vita_blue_trans_v6version_locked_24_traits.txt",
                         header = T)
accessions.pheno <- pheno.data$Taxa
accessions.non.pheno <- setdiff(accessions.keep, accessions.pheno)
accessions.keep <- intersect(accessions.keep, accessions.pheno)

# 545 accessions
blue.545 <- blue.664[accessions.keep]



# ---------------------------------------------------------------------------- #
# outlier removal
# ---------------------------------------------------------------------------- #
# Studentized residual
lmod <- lm(blue.545 ~ 1)
stud.res <- studres(lmod)

# numbers
N <- length(blue.545)
threshold <- qt(p =  1 - 0.05 / (2 * N), df = (N - 2))

# remove two accessions
blue.545.RmOut <- blue.545
blue.545.RmOut[stud.res > threshold] <- NA

# avg and BLUE
f <- paste0(dir.save, "/Fig005-Stud_resid.png")
png(filename = f, width = 600, height = 600)
plot(x = blue.545, y = stud.res,
		 xlab = "BLUE", ylab = "Studentized residual", 
		 main = "664 experimental lines",
		 pch = 20)
abline(h = threshold, lty = 2)
dev.off()



# ---------------------------------------------------------------------------- #
# TWAS
# ---------------------------------------------------------------------------- #
# load data
KinMat <- read.csv("RAWDATA/Kinship_NEW/GAPIT.Kin.VanRaden.csv", header = FALSE, stringsAsFactors = F)
Pheno.All <- read.delim("RAWDATA/Phenotype/ames_vita_blue_trans_v6version_locked_24_traits.txt")

# Peer residual matrix
myFun.01 <- function(VEC) {
	MIN <- min(VEC, na.rm = T)
	MAX <- max(VEC, na.rm = T)
	VEC.NEW <- (VEC - MIN) / (MAX - MIN)
	return(VEC.NEW)
}
blue.545.RmOut.sc <- myFun.01(blue.545.RmOut)

# kinship matrix
K <- as.matrix(KinMat[, -1])
rownames(K) <- colnames(K) <- KinMat$V1

# add pesuedo SNPs
pesuedo.SNPs <- matrix(runif(n = length(blue.545.RmOut.sc) * 3, min = 0, max = 1),
											 nr = length(blue.545.RmOut.sc),
											 nc = 3)
PeerResidMat.sc <- cbind(blue.545.RmOut.sc, pesuedo.SNPs)
colnames(PeerResidMat.sc) <- c("VTE7", "pesuedo1", "pesuedo2", "pesuedo3")

# LOOP FOR ALL TRAITS
trait.all <- c("a.T", "a.T3", "d.T", "d.T3", "g.T", "g.T3",
							 "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
pval.all <- rep(NA, length(trait.all))
names(pval.all) <- trait.all
for ( i in 1:length(trait.all) ) {
	# trait
	trait <- trait.all[i]
	
	# target phenotype data
	Pheno <- Pheno.All[, c("Taxa", trait)]
	
	# use intersect accessions
	common.accessions <- intersect(rownames(PeerResidMat.sc), Pheno$Taxa); length(common.accessions)
	my.Pheno <- Pheno[match(common.accessions, Pheno$Taxa), ]
	my.Geno <- data.frame("SNP" = NA, "Chr" = NA, "Pos" = NA, 
												t(PeerResidMat.sc[match(common.accessions, rownames(PeerResidMat.sc)), ]),
												row.names = 1:4)
	my.K <- K[match(common.accessions, rownames(K)), match(common.accessions, rownames(K))]
	
	# Run GWAS
	res.GWAS <- GWAS(pheno = my.Pheno, geno = my.Geno, K = my.K, n.PC = 0,
									 min.MAF = -Inf, P3D = FALSE, plot = FALSE, n.core = 1)
	colnames(res.GWAS)[4] <- "neg.log.P"
	colnames(res.GWAS)[1] <- "Gene"
	
	# save P-value
	pval <- res.GWAS$neg.log.P[1]
	pval.all[i] <- pval
}

# make figures
df.blue <- data.frame("Taxa" = names(blue.545.RmOut.sc),
											"vte7" = blue.545.RmOut.sc,
											row.names = NULL)
df.pheno <- Pheno.All[, c("Taxa", trait.all)]
df.fig <- merge(x = df.blue, y = df.pheno, by.x = "Taxa", by.y = "Taxa")
for ( i in 1:length(trait.all) ) {
	trait <- trait.all[i]
	f <- paste0(dir.save, "/Fig006-Stud_resid_", trait, ".png")
	png(filename = f, width = 600, height = 600)
	plot(x = df.fig$vte7, y = df.fig[, trait], pch = 20,
			 xlab = "vte7 (expression BLUE)", ylab = trait,
			 main = paste0("TWAS -log(P) = ", round(pval.all[i], 3)))
	dev.off()
}

# 
rank.prcnt.all <- rep(NA, length(pval.all))
for ( i in 1:length(pval.all) ) {
	trait <- names(pval.all)[i]
	f <- paste0("RESULT/4.2-TWAS/TwasResult_v1.1_B73_", trait, ".csv")
	TwasRes <- fread(f, data.table = F)
	pval.vte7 <- pval.all[trait]
	pval.all.snp <- TwasRes$neg.log.P
	rank.prcnt.all[i] <- 100 * sum(pval.all.snp > pval.vte7) / length(pval.all.snp)
}

# make file
df.save <- data.frame("neg.log.P" = round(pval.all, 3), 
											"rank.prcnt" = round(rank.prcnt.all, 2))
write.csv(df.save, file = paste0(dir.save, "/TwasSummary.csv"))






