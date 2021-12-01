# packages
library(gdata)
library(data.table)

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"

# file I/O
filename.blue <- paste0("RESULT/1.5-BLUE/expression_BLUE_RmErr_", ver, "_", ref, ".csv")
filename.out <- paste0("RESULT/1.5-BLUE/expression_BLUE_final_", ver, "_", ref, ".csv")
filename.summary <- paste0("RESULT/1.5-BLUE/Summary_Filtering_", ver, "_", ref, ".txt")

# load BLUE
blue <- fread(filename.blue, data.table = F)
blue.mat <-  as.matrix(blue[, -1])
rownames(blue.mat) <- blue[, 1]
colnames(blue.mat) <- gsub("_", "", colnames(blue.mat))

# expr.line vs check line
geno.names.all <- colnames(blue.mat)
geno.names.exp <- geno.names.all[substr(geno.names.all, 1, 5) != "Check"] # remove Check lines
n.check <- sum(substr(geno.names.all, 1, 5) == "Check")
text.b0 <- paste0("There are ", length(geno.names.all), " accessions in the BLUE matrix")
text.b1 <- paste0("There are ", n.check, " check lines in the BLUE matrix")
text.b2 <- paste0("There are ", length(geno.names.exp), " experimental lines in the BLUE matrix")

# remove check lines
blue.mat.rm.check <- t(blue.mat[, geno.names.exp])


# ---------------------------------------------------------------------------- #
# load accession info data
AccInfo <- read.xls("RAWDATA/Transcript/vitamaize_accession_removal_20191113.xlsx")
geno.names.info <- as.character(AccInfo$Accession.Number_ion)
text.a <- paste0("There are ", length(geno.names.info), " accessions in vitamaize_accession_removal_20191113.xlsx")

# common ones etc
geno.names.common <- intersect(geno.names.exp, geno.names.info)
text.c <- paste0(length(geno.names.common), " accessions out of ", 
                 length(geno.names.exp), " accessions are included in the vitamaize_accession_removal_20191113.xlsx")
text.d <- paste0("Accessions without remove/keep flag: ", setdiff(geno.names.exp, geno.names.common))
text.e <- "Ryokei and Xiaowei looked at GRIN and found that PI644099 is sweet corn."

# keep/remove flag -> accessions to be kept
m <- match(geno.names.common, geno.names.info)
tmp <- AccInfo[m, ]
tf.tmp.remove <- tmp$Removal_decision == "remove"
accessions.remove <- as.character(tmp$Accession.Number_ion)[tf.tmp.remove]
accessions.remove <- c(accessions.remove, "PI644099") # manually add PI644099
accessions.keep <- setdiff(geno.names.exp, accessions.remove)
text.f <- paste0("By using vitamaize_accession_removal_20191113.xlsx, we remove ", length(accessions.remove), " accessions")
text.g <- paste0("Therefore, in total of ", length(accessions.keep), " accessions were retained in this step")

# further filter accessions: if there is no phenotype, remove them
pheno.data <- read.table("RAWDATA/Phenotype/ames_vita_blue_trans_v6version_locked_24_traits.txt",
                         header = T)
accessions.pheno <- pheno.data$Taxa
accessions.non.pheno <- setdiff(accessions.keep, accessions.pheno)
accessions.keep <- intersect(accessions.keep, accessions.pheno)

# text to write
text.h <- paste0("There were ", length(accessions.non.pheno), " accessions without phenotype")
text.i <- paste0("There were ", length(accessions.keep), " accessions retaind in the end!")

# save filtered BLUE data
M <- t(as.matrix(blue[, -1]))
colnames(M) <- blue$GeneID
blue.save <- data.frame("Accession_ID" = rownames(M), M, stringsAsFactors = F)
m <- match(accessions.keep, gsub("_", "", blue.save$Accession_ID))
blue.save <- blue.save[m, ]
fwrite(blue.save, filename.out)
text.j <- paste0("The final BLUE matrix has ", nrow(blue.save), " samples and ", ncol(blue.save) - 1, " genes")

# write summary
text.summary <- c(text.a, text.b0, text.b1, text.b2, "", text.c, text.d, "",
                  text.e, text.f, text.g, "", text.h, text.i, "",
                  text.j)
write(x = text.summary, file = filename.summary, sep = "/n")
