# Merge Genes

# library
library(data.table)
library(IRanges)

# load data
GeneInfo <- fread("Zea_mays.B73_RefGen_v4.59_anno.csv", data.table = F)

# remove genes not on chromosome 1-10
tf <- GeneInfo$chr %in% 1:10
GeneInfo.OnChr <- GeneInfo[tf, ]

# check & remove duplication
df.DupCheck <- GeneInfo.OnChr[, c("chr", "start", "end", "strand")]
tf <- duplicated(df.DupCheck) | duplicated(df.DupCheck, fromLast = TRUE)
GeneInfo.Dup <- GeneInfo.OnChr[tf, ]

# Remove "EN..." from duplicated ones
tmp.vec <- paste(GeneInfo.Dup$chr, GeneInfo.Dup$start, GeneInfo.Dup$end, GeneInfo.Dup$strand, sep = "_")
tmp.unique <- unique(tmp.vec)
GeneID.remove.all <- c()
text.all <- c()
for (tmp in tmp.unique ) {
   tf <- tmp.vec == tmp
   GeneInfo.Dup.tmp <- GeneInfo.Dup[tf, ]
   GeneID.dup <- GeneInfo.Dup.tmp$ID
   if ( all(substr(GeneID.dup, 1, 2) == "Zm") ) {
      GeneID.remove <- GeneID.dup[GeneInfo.Dup.tmp$type != "gene"] # This is an exception!! BE CAREFULL!!
   } else {
      GeneID.remove <- GeneID.dup[substr(GeneID.dup, 1, 2) == "EN"]
   }
   GeneID.remove.all <- c(GeneID.remove.all, GeneID.remove)
   if ( length(GeneID.dup) == 2 ) {
      text.x <- paste(GeneID.dup[1], "and", GeneID.dup[2], "were duplicated, and ", GeneID.remove, "was removed")
      text.all <-c(text.all, text.x)
      print(text.x)
   } else {
      print("ERROR!!")
   }
}
GeneInfo.OnChr.retain <- GeneInfo.OnChr[!(GeneInfo.OnChr$ID %in% GeneID.remove.all), ]

# ----------- merge intervals
df.all <- NULL
for ( i in 1:10 ) {
   df.i <- GeneInfo.OnChr.retain[GeneInfo.OnChr.retain$chr == i, ]
   ir <- IRanges(df.i$start + 0.01, df.i$end - 0.01) # +/- 0.01 is an ad-hoc method to get a correct result
   vec <- subjectHits(findOverlaps(ir, reduce(ir)))
   df.i$group <- paste0(i, "-", vec)
   df.i$new.start <- tapply(df.i$start, vec, min)[vec]
   df.i$new.end <- tapply(df.i$end, vec, max)[vec]
   if ( i == 1 ) {
      df.all <- df.i
   } else {
      df.all <- rbind(df.all, df.i)
   }
}

# save
fwrite(x = df.all, file = "Zea_mays.B73_RefGen_v4.59_anno_merged.csv")

