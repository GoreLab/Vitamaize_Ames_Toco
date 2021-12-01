# Combined test

# library
library(data.table)
library(metap)
library(dplyr)

# mkdir
dir.create("RESULT/5.1-CombinedTest")

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"

# filename I/O
if ( ref == "B73" ) { file.gff <- paste0("RAWDATA/Annotation/Zea_mays.B73_RefGen_v4.59_anno_merged.csv") }
if ( ref == "PH207" ) { file.gff <- paste0("RAWDATA/Annotation/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0.v03eksc_anno_merged.csv") }
file.gwas.all <- "RAWDATA/GWAS/vita_v6_locked_hmp_LD0.1_K_GAPIT_MLM_all_top0.1_24_traits.csv"
file.pheno.all <- "RAWDATA/Phenotype/ames_vita_blue_trans_v6version_locked_24_traits.txt"

# load phenotype file (to get names of traits)
Pheno.All <- read.delim(file.pheno.all)
trait.all <- colnames(Pheno.All)[2:ncol(Pheno.All)]

# load GFF & GWAS result
GeneInfo <- fread(file.gff, data.table = F)
gwas.res.all <- fread(file.gwas.all, data.table = F)

# LOOP FOR ALL TRAITS
for ( trait in trait.all ) {
   # print
   print(paste0("Run FC test for ", trait))
   
   # filename I/O
   file.twas <- paste0("RESULT/4.2-TWAS/TwasResult_", ver, "_", ref, "_", trait, ".csv")
   file.save <- paste0("RESULT/5.1-CombinedTest/Result_FCtest_", ver, "_", ref, "_", trait, ".csv")
   
   # load TWAS result
   twas.res <- read.csv(file.twas, stringsAsFactors = F)
   twas.res <- twas.res[twas.res$Chr %in% 1:10, ] # retain chr = 1:10
   
   # GWAS result of the target trait
   gwas.res <- gwas.res.all[gwas.res.all$trait == trait, ]
   
   # 1. remake twas data (assign P = 1 for non-avairable genes)
   dat.twas.all <- data.frame("ID" = GeneInfo$ID,
                              "CHR" = GeneInfo$chr,
                              "START" = GeneInfo$start,
                              "END" = GeneInfo$end,
                              "GROUP.NEW" = GeneInfo$group,
                              "START.NEW" = GeneInfo$new.start,
                              "END.NEW" = GeneInfo$new.end,
                              "pvalue" = 1,
                              stringsAsFactors = F)
   m <- match(twas.res$Gene, GeneInfo$ID)
   dat.twas.all$pvalue[m] <- 10 ^ -(twas.res$neg.log.P) # p-values in TWAS
   
   # 2. MERGE TWAS + GWAS
   snp.border <- c()
   df.pair.all <- NULL
   # for all chromosomes
   for ( i in 1:10 ) {
      # print
      b <- Sys.time()
      print(paste0("Merge GWAS and TWAS result for the ", i, "-th chromosome: start"))
      
      # result on i-th chr
      dat.twas.chr <- dat.twas.all[dat.twas.all$CHR == i, ]
      gwas.res.chr <- gwas.res[gwas.res$Chromosome == i, ]
      gwas.res.chr <- gwas.res.chr[order(gwas.res.chr$Position), ]
      P <- nrow(dat.twas.chr)
      N <- nrow(gwas.res.chr)
      
      # make intervals
      dat.twas.chr.unique.grp <- dat.twas.chr[!duplicated(dat.twas.chr$GROUP.NEW), ]
      P2 <- nrow(dat.twas.chr.unique.grp)
      mid.pos <- (dat.twas.chr.unique.grp$START.NEW[2:P2] + dat.twas.chr.unique.grp$END.NEW[1:(P2-1)]) / 2
      tmp.vec <- c(0, rep(mid.pos, each = 2), Inf)
      interval.matrix <- t(matrix(tmp.vec, nr = 2))
      rownames(interval.matrix) <- dat.twas.chr.unique.grp$GROUP.NEW
      colnames(interval.matrix) <- c("START", "END")
      
      # assign SNPs
      SNP.all.a <- group.all.a <- rep(NA, N)
      SNP.all.b <- group.all.b <- c() # for the border of two intervals
      for ( j in 1:N ) {
         SNP.j <- gwas.res.chr$SNP[j]
         pos.j <- gwas.res.chr$Position[j]
         
         # find which bin?
         if ( j == 1 ) { count <- 0 } else { count <- num - 1 }
         repeat {
            count <- count + 1
            l <- interval.matrix[count, 1]
            r <- interval.matrix[count, 2]
            tf <- between(x = pos.j, left = l, right = r)
            if ( tf ) {
               num <- count
               break
            }
         }
         
         # make vec of result pairs 
         SNP.all.a[j] <- SNP.j
         group.all.a[j] <- rownames(interval.matrix)[num]
         
         # check border
         tf.border <- pos.j == interval.matrix[num, 2]
         if ( tf.border ) {
            SNP.all.b <- c(SNP.all.b, SNP.j)
            group.all.b <- c(group.all.b, rownames(interval.matrix)[num+1])
            print(paste0("The ", j, "-th marker (pos = ", pos.j,") was on the border."))
            snp.border <- c(snp.border, SNP.j)
         }
      }
      
      # make data frame
      df.tmp <- data.frame("SNP" = c(SNP.all.a, SNP.all.b), "group" = c(group.all.a, group.all.b))
      df.tmp.2 <- merge(x = df.tmp, y = dat.twas.chr, by.x = "group", by.y = "GROUP.NEW")
      df.pair <- merge(x = df.tmp.2, y = gwas.res.chr, by.x = "SNP", by.y = "SNP")
      df.pair <- df.pair[order(df.pair$Position), ] # re-order
      df.pair <- df.pair[, c("SNP", "Chromosome", "Position", "P.value", "FDR_Adjusted_P-values",
                             "ID", "START", "END", "pvalue",
                             "group", "START.NEW", "END.NEW", "trait")]
      colnames(df.pair) <- c("snp", "chr", "pos", "gwas.pval", "gwas.adj.pval",
                             "gene", "start", "end", "twas.pval",
                             "gene.group.id", "start.merged", "end.merged", "trait")
      
      # merge data.frame
      df.pair.all <- rbind(df.pair.all, df.pair)
      
      # print
      a <- Sys.time()
      print(paste0("Merge GWAS and TWAS result for the ", i, "-th chromosome: done"))
      print(a-b)
   }
   
   # 3. FC.TEST
   pval.mat <- cbind(df.pair.all$gwas.pval, df.pair.all$twas.pval)
   combined.p <- apply(pval.mat, 1, FUN = function(x){sumlog(x)$p})
   df.pair.all$FCT.pval <- combined.p
   
   # 4. re-format result
   df.save <- df.pair.all[, c("snp", "chr", "pos", "gene", "start", "start.merged", "end", "end.merged", "gwas.pval", "twas.pval", "FCT.pval", "trait")]
   
   # save result
   fwrite(x = df.save, file = file.save)
}

