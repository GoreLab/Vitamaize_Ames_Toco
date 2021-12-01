# calculate and add GDD to the data

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"

# file I/O
filename.load <- paste0("RESULT/1.1-MakeDataset/SampleInfo_", ver, "_", ref, ".csv")
filename.save <- paste0("RESULT/1.2-AttachGDD/SampleInfo_withGDD_", ver, "_", ref, ".csv")

# mkdir
dir.create("RESULT/1.2-AttachGDD")

# read info data
info <- read.csv(filename.load)
info$poll.date <- as.Date(info$poll.date)
info$harv.date <- as.Date(info$harv.date)
head(info)

# read temp data
weather.dat <- read.csv("RAWDATA/HarvestDate/weather_2018.csv")

# get min/max temp of each day
date.vec <- substr(weather.dat$valid, 1, 10)
T.min <- tapply(X = weather.dat$tmpf, INDEX = date.vec, FUN = min)
T.max <- tapply(X = weather.dat$tmpf, INDEX = date.vec, FUN = max)

# AGDD
AGDD.all <- rep(NA, times = nrow(info))
for ( i in 1:nrow(info) ) {
   date.vec.i <- seq(from = info$poll.date[i] + 1, to = info$harv.date[i] - 1, by = "days") # remove boder date!!
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
info$GDD <- AGDD.all # attach

# save
write.csv(info, file = filename.save, row.names = F)

