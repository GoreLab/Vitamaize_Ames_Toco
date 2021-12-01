# package
library(data.table)

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"

# filenames Input
filename.info <- paste0("RESULT/1.2-AttachGDD/SampleInfo_withGDD_", ver, "_", ref, ".csv")
filename.blue <- paste0("RESULT/1.5-BLUE/expression_BLUE_", ver, "_", ref, "_raw.csv")
filename.lastm <- paste0("RESULT/1.5-BLUE/expression_BLUE_", ver, "_", ref, "_message.csv")
filename.out <- paste0("RESULT/1.5-BLUE/expression_BLUE_RmErr_", ver, "_", ref, ".csv")
filename.summary <- paste0("RESULT/1.5-BLUE/Summary_RmErr_", ver, "_", ref, ".txt")

# load data
Info <- read.csv(filename.info, stringsAsFactors = F)
BLUE <- fread(filename.blue)
LastM <- fread(filename.lastm)

# add column
Info$GENO <- Info$GENOTYPE
Info$GENO[Info$CHECK == "Check"] <- Info$GRP[Info$CHECK == "Check"]

# BLUE
BLUE.mat <- as.matrix(BLUE[, -1])
rownames(BLUE.mat) <- unlist(BLUE[, 1])

# remove genes which we get error message
tf <- LastM$last.m == "LogLikelihood Converged"
BLUE.RmOut <- BLUE.mat[tf, ]

# save final BLUE data
df.save <- data.frame("GeneID" = rownames(BLUE.RmOut), BLUE.RmOut)
fwrite(df.save, filename.out)

# make summary
n.sample.all <- ncol(BLUE.mat)
n.gene.all <- nrow(LastM)
n.gene.conv <- sum(LastM$last.m == "LogLikelihood Converged")
lastm.non.conv <- LastM$last.m[LastM$last.m != "LogLikelihood Converged"]
table.lastm.non.conv <- table(lastm.non.conv)
text.err <- c()
if ( length(table.lastm.non.conv) != 0 ) {
   for ( i in 1:length(table.lastm.non.conv) ) {
      text.i <- paste0(names(table.lastm.non.conv[i]), ": ",  table.lastm.non.conv[i])
      text.err <- c(text.err, text.i)
   }
}

# write summary
text.summary <- c(paste0("Summary numbers for the BLUE calcualtion of ", ref, "-mapped, ", ver, " dataset"), 
                  "",
                  paste0("The raw BLUE matrix has ", n.sample.all, " samples and ", n.gene.all, " genes"),
                  "",
                  paste0("Sumamry of the error message:"),
                  paste0("LogLikelihood Converged: ", n.gene.conv),
                  text.err,
                  "",
                  paste0("The BLUE matrix has ", ncol(BLUE.RmOut), " samples and ", nrow(BLUE.RmOut),
                         " genes, after removing errors")
                  )
write(x = text.summary, file = filename.summary, sep = "/n")
