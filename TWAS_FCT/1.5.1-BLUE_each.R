# calculate BLUE via asreml

# get argument from shell
args <- commandArgs(trailingOnly = T)
s <- as.numeric(args[1])
ver <- args[2] # "v1" or "v1.1"
ref <- args[3] # "B73" or "PH207"
rand.seed <- 2018

# source
library(asreml)
library(data.table)
library(doParallel)

# file I/O
filename.data <- paste0("RESULT/1.4-OutlierRemoval_and_Imputation/ExpressionData_For_BLUE_", ver, "_", ref, ".csv")
filename.info <- paste0("RESULT/1.2-AttachGDD/SampleInfo_withGDD_", ver, "_", ref, ".csv")
foldername.save <- paste0("RESULT/1.5-BLUE/ResEach_", ver, "_", ref)
filename.save.1 <- paste0("RESULT/1.5-BLUE/ResEach_", ver, "_", ref, "/exp_BLUE_all_", s, ".Rdata")
filename.save.2 <- paste0("RESULT/1.5-BLUE/ResEach_", ver, "_", ref, "/exp_BLUE_pred_", s, ".Rdata")
filename.save.3 <- paste0("RESULT/1.5-BLUE/ResEach_", ver, "_", ref, "/exp_BLUE_lastm_", s, ".Rdata")
filename.save.4 <- paste0("RESULT/1.5-BLUE/ResEach_", ver, "_", ref, "/exp_BLUE_seed_", s, ".txt")

# make folder to save result
dir.create(foldername.save, recursive = TRUE)

# load data
info <- read.csv(filename.info, row.names = 1)
exp.data <- fread(filename.data)

# re-shape: make matrix
sample.name <- exp.data[, 1]
exp.mat <- as.matrix(exp.data[, -1])
rownames(exp.mat) <- unlist(sample.name)

# make data
y <- exp.mat[, s]
dat <- cbind(info, "y" = y)

# add "GENO" column to the data
dat$GENO <- as.character(dat$GENOTYPE)
dat$GENO[dat$CHECK == "Check"] <- as.character(dat$GRP[dat$CHECK == "Check"])
dat$GENO <- as.factor(dat$GENO)

# int(and char, just in case) -> factor
dat$TIER <- as.factor(dat$TIER)
dat$BLOCK <- as.factor(dat$BLOCK)

# run asreml
set.seed(rand.seed)
asr.mod <- asreml(fixed = y ~ -1 + GENO + GDD, random = ~ TIER + BLOCK + PLATE + PLATE:LANE, maxiter = 100, data = dat, trace = F)
last.m <- asr.mod$last.message
pred <- predict(asr.mod, classify = "GENO", data = dat)
pred.val <- pred$predictions$pvals$predicted.value
names(pred.val) <- pred$predictions$pvals$GENO

# make object for the output
coef.all <- coef(asr.mod)
coef.fixed <- setNames(coef.all$fixed, rownames(coef.all$fixed))
coef.random <- setNames(coef.all$random, rownames(coef.all$random))
gamma.rand <- asr.mod$gammas
asr.summary.all <- c(coef.fixed, coef.random, gamma.rand)

# write
saveRDS(asr.summary.all, file = filename.save.1)
saveRDS(pred.val, file = filename.save.2)
saveRDS(last.m, file = filename.save.3)
write(rand.seed, file = filename.save.4)

# print
if ( asr.mod$converge == FALSE ) { print(s); print("Not converged") }
