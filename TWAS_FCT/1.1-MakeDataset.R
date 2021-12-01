# Re-format raw data for the next steps
# Also write out a summary numbers (# of samples, # of positive controls etc.)

# packages
library(data.table)
library(e1071)
library(gdata)
library(tibble)

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"

# TagSeq filename to load
if ( ver == "v1.1" & ref == "B73" ) { filename.TagSeq <- "RAWDATA/Transcript/B73/RLOG_curated_count_matrix_all_info_v1.1.txt" }
if ( ver == "v1" & ref == "B73" ) { filename.TagSeq <- "RAWDATA/Transcript/B73/RLOG_curated_count_matrix_all_info_v1.txt" }
if ( ver == "v1.1" & ref == "PH207" ) { filename.TagSeq <- "RAWDATA/Transcript/PH207/PH207_RLOG_curated_count_matrix_all_info_v1.1.txt" }
if ( ver == "v1" & ref == "PH207" ) { filename.TagSeq <- "RAWDATA/Transcript/PH207/PH207_RLOG_curated_count_matrix_all_info_v1.txt" }

# other filenames to load
filename.field <- "RAWDATA/Transcript/Source_to_Block_key_18_REVISED.csv"
filename.harvest.all <- "RAWDATA/HarvestDate/HarvestInfo_all.xlsx"
filename.harvest.check <- "RAWDATA/HarvestDate/check_harvest_and_pollination_dates_9apr20.xlsx"
filename.harvest.double <- "RAWDATA/HarvestDate/double_harvested_samples_correct_harvest_and_pollination_date_9Apr20.xlsx"

# filename to save
filename.info <- paste0("RESULT/1.1-MakeDataset/SampleInfo_", ver, "_", ref, ".csv")
filename.data <- paste0("RESULT/1.1-MakeDataset/ExpressionData_", ver, "_", ref, ".csv")
figname <- paste0("RESULT/1.1-MakeDataset/GeneFiltering_", ver, "_", ref, ".png")
filename.summary <- paste0("RESULT/1.1-MakeDataset/Summary_", ver, "_", ref, ".txt")

# mkidr
folder.save <- "RESULT/1.1-MakeDataset"
dir.create(folder.save, recursive = TRUE)

# load data
TagSegDat.raw <- fread(filename.TagSeq)
field.map <- read.csv(filename.field, stringsAsFactors = FALSE)
harvest.all.raw <- read.xls(filename.harvest.all)
harvest.check.raw <- read.xls(filename.harvest.check)
harvest.double.raw <- read.xls(filename.harvest.double)

# re-format
if ( ref == "B73" ) { num.info.col <- 5 } else { num.info.col <- 6 } # 5 cols in B73, 6 cols in PH207
gene.id <- TagSegDat.raw$gene_id
gene.info <- as.data.frame(TagSegDat.raw[, 2:num.info.col])
sample.names <- colnames(TagSegDat.raw)[(num.info.col + 1):ncol(TagSegDat.raw)]
dat.matrix <- t(as.matrix(TagSegDat.raw[, (num.info.col + 1):ncol(TagSegDat.raw)]))
colnames(dat.matrix) <- gene.id
harvest.all <- data.frame("Source" = harvest.all.raw$Source,
                          "RecId" = harvest.all.raw$RecId,
                          "Entry" = harvest.all.raw$Entry,
                          "Range" = harvest.all.raw$Range,
                          "Pass" = harvest.all.raw$Pass,
                          "Pedigree" = harvest.all.raw$Pedigree,
                          "poll.date" = harvest.all.raw$poll.date..of.fresh.harvest.ear.,
                          "harv.date" = harvest.all.raw$harv.,
                          stringsAsFactors = F)
harvest.check <- harvest.check.raw[, 2:13]

# remove controls
tf <- sapply(sample.names, FUN = function(x){substr(x, nchar(x)-15, nchar(x)) == "Positive_Control"})
dat.matrix.samples <- dat.matrix[!tf, ]; dim(dat.matrix.samples)

# summary number
n.all.samples <- length(sample.names)
n.positive.control <- sum(tf)
n.real.sample <- sum(!tf)
n.all.genes <- ncol(dat.matrix)

# --- decompose sample name: easy for BLUE/BLUP etc
MSU.ID <- PLATE <- WELL <- ACC.ID <- ACC.NAME <- TIER <- ROW <- PASS <- c()
sample.names.new <- rownames(dat.matrix.samples)
for ( i in 1:length(sample.names.new) ) {
   # target name
   x <- sample.names.new[i]
   
   # 01. get field info
   x.sep01 <- strsplit(x, "_Field:")[[1]]
   field.info <- strsplit(x.sep01[2], "_")[[1]]
   tier.num <- as.numeric(gsub("T", "", field.info[1]))
   row.num <- as.numeric(gsub("R", "", field.info[2]))
   pass.num <- as.numeric(gsub("P", "", field.info[3]))
   
   # 02. get accession.id
   x.sep02 <- strsplit(x.sep01[1], "_Acc:")[[1]]
   acc.id <- x.sep02[2]
   
   # 03. get accession.name
   x.sep03 <- strsplit(x.sep02[1], "_Ped:")[[1]]
   acc.name <- x.sep03[2]
   
   # 04. get MSU-ID
   x.sep04 <- strsplit(x.sep03[1], "_GT_")[[1]]
   msu.id <- paste0("GT_", x.sep04[2])
   
   # 04. get MSU-ID
   x.info <- x.sep04[1]
   plate <- strsplit(x.info, "_")[[1]][1]
   well <- strsplit(x.info, "_")[[1]][2]
   
   # save
   MSU.ID <- c(MSU.ID, msu.id)
   PLATE <- c(PLATE, plate)
   WELL <- c(WELL, well)
   ACC.ID <- c(ACC.ID, acc.id)
   ACC.NAME <- c(ACC.NAME, acc.name)
   TIER <- c(TIER, tier.num)
   ROW <- c(ROW, row.num)
   PASS <- c(PASS, pass.num)
}
# here is the information from the raw data
df.info <- data.frame("MSU.ID" = MSU.ID,
                      "PLATE" = PLATE,
                      "WELL" = WELL,
                      "ACCESSION.ID" = ACC.ID,
                      "ACCESSION.NAME" = ACC.NAME,
                      "TIER" = TIER,
                      "ROW" = ROW,
                      "PASS" = PASS,
                      stringsAsFactors = F)

# assign field info
info.trp <- paste(df.info$TIER, df.info$ROW, df.info$PASS, sep = "-")
field.trp <- paste(field.map$Tier, field.map$Range, field.map$Pass, sep = "-")
length(info.trp) == length(unique(info.trp)) # not unique (ok)
length(field.trp) == length(unique(field.trp)) # unique
m <- match(info.trp, field.trp)
df.info <- cbind(df.info, field.map[m, c("Range", "Block", "Grp", "Pedigree.GRIN", "Source.GRIN")])
all(df.info$ROW == df.info$Range) # check -> OK

# remove redundant column
df.info <- data.frame("MSU.ID" = df.info$MSU.ID,
                      "TIER" = df.info$TIER,
                      "BLOCK" = df.info$Block,
                      "RANGE" = df.info$Range,
                      "PASS" = df.info$PASS,
                      "PLATE" = df.info$PLATE,
                      "WELL" = df.info$WELL,
                      "GENOTYPE" = df.info$ACCESSION.ID,
                      "PEDIGREE" = df.info$ACCESSION.NAME,
                      "GRP" = df.info$Grp,
                      stringsAsFactors = F)

table(df.info$GRP)
df.info[df.info$GRP == "Check1", ]


# add column CHECK (not necessary)
df.info$CHECK <- "New_Lines"
df.info$CHECK[df.info$GRP != "New_Lines"] <- "Check"
head(df.info)

# add column LANE
df.info <- add_column(df.info, LANE = NA, .after = "WELL") # assign L1 for all (temporary)
for ( i in 1:nrow(df.info) ) {
   plate.i <- df.info$PLATE[i]
   well.i <- df.info$WELL[i]
   well.num.i <- substr(as.character(well.i), 2, 3)
   if ( plate.i %in% c("P1", "P15") ) {
      df.info$LANE[i] <- "L1" # Plate 1 was not divided into two lanes. So, all samples should be lane 1. Plate 15 as well.
   } else {
      if ( well.num.i %in% c("07", "08", "09", "10", "11", "12") ) {
         df.info$LANE[i] <- "L2"
      } else {
         df.info$LANE[i] <- "L1"
      }
   }
}

# assign pollination & harvest date
rp.info <- paste(df.info$RANGE, df.info$PASS, sep = "-")
rp.harvest <- paste(harvest.all$Range, harvest.all$Pass, sep = "-")
m <- match(rp.info, rp.harvest)
harvest.all.sorted <- harvest.all[m, c("Source", "RecId", "poll.date", "harv.date")]
harvest.all.sorted$poll.date <- as.character(harvest.all.sorted$poll.date)
harvest.all.sorted$harv.date <- as.character(harvest.all.sorted$harv.date)
head(df.info); head(harvest.all[m, ])
df.info <- cbind(df.info, harvest.all.sorted[, c("poll.date", "harv.date")])

# add pollination & harvest date for those in Tier 4
num <- which(df.info$poll.date == "")
for ( i in num ) {
   msu.id.i <- df.info[i, "MSU.ID"]
   harvest.check.i <- harvest.check[harvest.check$MSU_Sample_Name %in% msu.id.i, ]
   df.info[i, "harv.date"] <- as.character(harvest.check.i$Harvest.Date)
   df.info[i, "poll.date"] <- as.character(harvest.check.i$Pollination.Date)
}

# add pollination & harvest date for those harvested twice
m <- match(harvest.double.raw$MSU_Sample_Name, df.info$MSU.ID)
df.info[m, "poll.date"] <- as.character(harvest.double.raw$Correct.Pollination.Date)
df.info[m, "harv.date"] <- as.character(harvest.double.raw$Correct.Harvest.Date)

# visual check
unique(df.info$poll.date)
unique(df.info$harv.date) # one wrong encoding "8/21"
df.info[df.info$harv.date == "8/21", "harv.date"] <- "2018-08-21" # manual correction

# info: visual check
head(df.info)

# summary number
n.exp <- sum(df.info$CHECK == "New_Lines")
n.check <- sum(df.info$CHECK == "Check")

# save
write.csv(dat.matrix.samples, file = filename.data)
write.csv(df.info, file = filename.info, row.names = F)

# write summary
text.summary <- c(paste0("Summary numbers for the ", ref, "-mapped, ", ver, "dataset"), 
                  "",
                  "Summary of samples",
                  paste0("number of all samples in the raw data: ", n.all.samples),
                  paste0("number of positive controls: ", n.positive.control),
                  paste0("number of plant samples: ", n.real.sample),
                  paste0("number of experimental-line samples: ", n.exp),
                  paste0("number of check samples: ", n.check),
                  "",
                  "Summary of genes",
                  paste0("number of all genes: ", n.all.genes))
write(x = text.summary, file = filename.summary, sep = "/n")
