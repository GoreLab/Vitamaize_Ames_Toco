# source
library(data.table)
library(rrBLUP)

# function for BIC calculation
myfun.calc.BIC <- function(y, PC, K) {
   # number of max. PC
   n.max.pc <- ncol(PC)
   
   # remove NA
   tf.na <- !is.na(y)
   y.common <- y[tf.na]
   K.common <- K[tf.na, tf.na]
   PC.common <- PC[tf.na, ]
   
   # ----- BIC with K ----- #
   # fit model by using mixed.solve function
   LL.vec <- rep(NA, n.max.pc + 1)
   names(LL.vec) <- paste0("PC", 0:n.max.pc)
   for ( k in 0:n.max.pc ) {
      if ( k == 0 ) { X.matrix <- NULL }
      if ( k >= 1 ) { X.matrix <- cbind(1, PC.common[, 1:k]) }
      ms <- mixed.solve(y = y.common, Z = NULL, K = K.common, X = X.matrix, method = "REML")
      n <- length(y.common)
      q <- length(ms$beta)
      H <- K.common + diag(nrow(K.common)) * (ms$Ve / ms$Vu)
      H.inv <- MASS::ginv(H)
      a <- (-1) * n * log(2 * pi * ms$Vu)
      b <- (-1) * as.numeric(determinant(H, logarithm = TRUE)$modulus)
      if ( k == 0 ) {
         X.matrix <- matrix(1, nr = n, nc = 1)
         c <- as.numeric((-1) * (1 / ms$Vu) * t(y.common - X.matrix %*% ms$beta) %*% H.inv %*% (y.common - X.matrix %*% ms$beta))
      } else {
         c <- as.numeric((-1) * (1 / ms$Vu) * t(y.common - X.matrix %*% ms$beta) %*% H.inv %*% (y.common - X.matrix %*% ms$beta))
      }
      LogLik.ML <- 0.5 * (a + b + c) # full likelihood based on REML 
      LL.vec[k+1] <- LogLik.ML
   }
   n <- length(y.common)
   n.par.vec <- seq(from = 3, by = 1, length.out = n.max.pc + 1) # n.par starts from 3 & increases by 1
   BIC.vec <- -2 * LL.vec + n.par.vec * log(n)
   
   # ----- BIC without K ----- #
   LL.vec.woK <- rep(NA, n.max.pc + 1)
   names(LL.vec.woK) <- paste0("PC", 0:n.max.pc)
   for ( k in 0:n.max.pc ) {
      if ( k == 0 ) {
         LL.vec.woK[k+1] <- as.numeric(logLik(lm(y.common ~ 1)))
      } else {
         df.lm <- data.frame("y" = y.common, PC.common[, 1:k])
         colnames(df.lm)[2:ncol(df.lm)] <- paste0("PC", 1:k)
         LL.vec.woK[k+1] <- as.numeric(logLik(lm(y ~., df.lm)))
      }
   }
   n.par.vec.woK <- seq(from = 2, by = 1, length.out = n.max.pc + 1) # n.par starts from 2
   BIC.vec.woK <- n.par.vec.woK * log(n) - 2 * LL.vec.woK
   
   # return
   ResList <- list("n.input" = length(y),
                   "n.sample" = length(y.common),
                   "BIC" = list("withK" = data.frame("n.PC" = names(LL.vec),
                                                     "LL" = LL.vec,
                                                     "n.par" = n.par.vec,
                                                     "BIC" = BIC.vec),
                                "withoutK" = data.frame("n.PC" = names(LL.vec.woK),
                                                        "LL" = LL.vec.woK,
                                                        "n.par" = n.par.vec.woK,
                                                        "BIC" = BIC.vec.woK)))
   return(ResList)
}

# mkdir
dir.save <- "RESULT/4.1-BIC_for_TWAS"
dir.create(dir.save, recursive = T)

# load data
KinMat <- read.csv("RAWDATA/Kinship_NEW/GAPIT.Kin.VanRaden.csv", header = FALSE, stringsAsFactors = F)
Pheno.All <- read.delim("RAWDATA/Phenotype/ames_vita_blue_trans_v6version_locked_24_traits.txt")

# kinship matrix
K <- as.matrix(KinMat[, -1])
rownames(K) <- colnames(K) <- KinMat$V1

# calculate PCs
pr <- prcomp(K)
PC <- pr$x[, 1:10]

# use intersect accessions
m <- match(rownames(K), Pheno.All$Taxa)
pheno.data <- Pheno.All[m, ]

# traits
trait.all <- colnames(pheno.data)[2:ncol(pheno.data)]
df.all <- NULL
for ( i in 1:length(trait.all) ) {
   # use a trait
   trait <- trait.all[i]
   
   # target phenotype data
   pheno <- pheno.data[, c("Taxa", trait)]
   
   # calculate
   ResList <- myfun.calc.BIC(y = pheno[, 2], PC = PC, K = K)
   
   # make df to merge all results
   df.i <- data.frame("Trait" = trait,
                      "Kinship" = rep(c("With.K", "Without.K"), each = nrow(ResList$BIC$withK)),
                      "PC" = c(rownames(ResList$BIC$withK), 
                               rownames(ResList$BIC$withoutK)),
                      "BIC" = c(ResList$BIC$withK$BIC, ResList$BIC$withoutK$BIC))
   df.all <- rbind(df.all, df.i)
}

# save
write.csv(df.all, file = paste0(dir.save, "/result_BIC_myFun.csv"), row.names = F)

# make figure
library(ggplot2)
df.all$n.PC <- as.integer(gsub("PC", "", df.all$PC))
p <- ggplot(data = df.all, aes(x = n.PC, y = BIC))
p <- p + facet_grid(Trait ~ Kinship, scales = "free_y")
p <- p + geom_point() + geom_line()
ggsave(filename = paste0(dir.save, "/Figure.png"), p, width = 5, height = 20)
