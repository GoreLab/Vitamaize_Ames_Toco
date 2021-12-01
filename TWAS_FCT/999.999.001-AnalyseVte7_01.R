# package
library(data.table)
library(tibble)
library(gdata)
library(asreml)

# library
dir.save <- "RESULT/999-VTE7"
dir.create(dir.save, recursive = T)

# ---------------------------------------------------------------------------- #
# Re-format raw data for the next steps
# ---------------------------------------------------------------------------- #
# TagSeq filename to load
filename.TagSeq <- "RAWDATA/VTE7/Zm00001d006778_Zm00001d006779_merged_gene_CPM_matrix_v1.1.txt"
filename.field <- "RAWDATA/Transcript/Source_to_Block_key_18_REVISED.csv"
filename.harvest.all <- "RAWDATA/HarvestDate/HarvestInfo_all.xlsx"
filename.harvest.check <- "RAWDATA/HarvestDate/check_harvest_and_pollination_dates_9apr20.xlsx"
filename.harvest.double <- "RAWDATA/HarvestDate/double_harvested_samples_correct_harvest_and_pollination_date_9Apr20.xlsx"

# load data
TagSegDat.raw <- read.table(filename.TagSeq, skip = 2)
field.map <- read.csv(filename.field, stringsAsFactors = FALSE)
harvest.all.raw <- read.xls(filename.harvest.all)
harvest.check.raw <- read.xls(filename.harvest.check)
harvest.double.raw <- read.xls(filename.harvest.double)

# re-format
colnames(TagSegDat.raw) <- c("Sample.ID", "Count", "CPM")
TagSegDat.raw$Sample.ID <- as.character(TagSegDat.raw$Sample.ID)
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
sample.names <- as.character(TagSegDat.raw$Sample.ID)
tf <- sapply(sample.names, FUN = function(x){substr(x, nchar(x)-15, nchar(x)) == "Positive_Control"})
TagSegDat.samples <- TagSegDat.raw[!tf, ]

# decompose sample name: easy for BLUE/BLUP etc
MSU.ID <- PLATE <- WELL <- ACC.ID <- ACC.NAME <- TIER <- ROW <- PASS <- c()
sample.names.new <- as.character(TagSegDat.samples$Sample.ID)
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
   
   # 05. get plate/well
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

# MANUAL CURATION
df.info[df.info$MSU.ID == "GT_425", "PASS"] <- 84

# assign field info
info.trp <- paste(df.info$TIER, df.info$ROW, df.info$PASS, sep = "-")
field.trp <- paste(field.map$Tier, field.map$Range, field.map$Pass, sep = "-")
m <- match(info.trp, field.trp)
df.info <- cbind(df.info, field.map[m, c("Range", "Block", "Grp", "Pedigree.GRIN", "Source.GRIN")])

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

# add column CHECK (not necessary)
df.info$CHECK <- "New_Lines"
df.info$CHECK[df.info$GRP != "New_Lines"] <- "Check"

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
df.info[df.info$harv.date == "8/21", "harv.date"] <- "2018-08-21" # manual correction



# ---------------------------------------------------------------------------- #
# attach GDD
# ---------------------------------------------------------------------------- #
# read temp data
weather.dat <- read.csv("RAWDATA/HarvestDate/weather_2018.csv")

# convert to date object
df.info$poll.date <- as.Date(df.info$poll.date)
df.info$harv.date <- as.Date(df.info$harv.date)

# get min/max temp of each day
date.vec <- substr(weather.dat$valid, 1, 10)
T.min <- tapply(X = weather.dat$tmpf, INDEX = date.vec, FUN = min)
T.max <- tapply(X = weather.dat$tmpf, INDEX = date.vec, FUN = max)

# AGDD
AGDD.all <- rep(NA, times = nrow(df.info))
for ( i in 1:nrow(df.info) ) {
   date.vec.i <- seq(from = df.info$poll.date[i] + 1, to = df.info$harv.date[i] - 1, by = "days") # remove boder date!!
   T.min.vec.i <- T.min[as.character(date.vec.i)]
   T.max.vec.i <- T.max[as.character(date.vec.i)]
   T.min.vec.i
   T.min.vec.i[T.min.vec.i < 50] <- 50 # should be greater than 50
   T.max.vec.i[T.max.vec.i > 86] <- 86 # should be smaller than 86
   GDD.vec.i <- (T.min.vec.i + T.max.vec.i) / 2 - 50
   AGDD.vec.i <- cumsum(GDD.vec.i)
   AGDD.at.harvest.i <- as.numeric(tail(AGDD.vec.i, 1)) # last one is the final value
   AGDD.all[i] <- AGDD.at.harvest.i
}
df.info$GDD <- AGDD.all # attach

# save
fwrite(x = df.info, file = paste0(dir.save, "/SampleInfo.csv"))
fwrite(x = TagSegDat.samples, file = paste0(dir.save, "/ExprData.csv"))

# figure 001
f <- paste0(dir.save, "/Fig001-Hist_Count.png")
png(filename = f, width = 600, height = 500)
hist(x = TagSegDat.samples$Count, 
     xlab = "Total Counts in Regions of Interest (n = 2)",
     main = paste0("Histogram of ", nrow(TagSegDat.samples), " samples"),
     breaks = 20)
dev.off()

# figure 002
f <- paste0(dir.save, "/Fig002-Hist_CPM.png")
png(filename = f, width = 600, height = 500)
hist(x = TagSegDat.samples$CPM, 
     xlab = "CPM",
     main = paste0("Histogram of ", nrow(TagSegDat.samples), " samples"),
     breaks = 20)
dev.off()

# figure 003
f <- paste0(dir.save, "/Fig003-Count_and_CPM.png")
png(filename = f, width = 600, height = 600)
plot(x = TagSegDat.samples$Count, 
     y = TagSegDat.samples$CPM,
     xlab = "Total Counts in Regions of Interest (n = 2)",
     ylab = "CPM",
     main = paste0(nrow(TagSegDat.samples), " samples in v1.1 data"),
     pch = 20)
dev.off()



# ---------------------------------------------------------------------------- #
# MAD-based outlier removal (there is no outlier!)
# ---------------------------------------------------------------------------- #
# MAD calculation
x <- TagSegDat.samples$CPM
thres.x <- mad(x, constant = 1)
mad.x <- abs(x - median(x)) / thres.x

# data frame to save
df.save <- data.frame("Sample.ID" = TagSegDat.samples$Sample.ID,
                      "CPM" = TagSegDat.samples$CPM,
                      "MAD" = mad.x)

# save
fwrite(df.save, file = paste0(dir.save, "/MAD.csv"))

# ---------------------------------------------------------------------------- #
# BLUE calculation
# ---------------------------------------------------------------------------- #
# set seed
rand.seed <- 2021

# make data
y <- TagSegDat.samples$CPM
dat <- cbind(df.info, "y" = y)

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

# data frame to save
df.save <- data.frame("Genotype" = names(pred.val),
                      "BLUE" = pred.val)

# save
fwrite(df.save, file = paste0(dir.save, "/BLUE.csv"))