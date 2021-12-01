# make scree plot

# mkdir
dir.create("RESULT/2.2-MakeScreePlot")

# 
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"

# make figure
file.in <- paste0("RESULT/2.1-Peer_Use25Fact/PeerResult_Use25Fact_", ver, "_", ref, "_precision.txt")
file.out <- paste0("RESULT/2.2-MakeScreePlot/ScreePlot_Use25Fact_", ver, "_", ref, ".jpeg")
precision <- read.delim(file.in)
variacne.from.second <- 1 / precision$precision[-1]
jpeg(filename = file.out)
plot(x = 2:25, y = variacne.from.second,
     type = "l",
     lwd = 2,
     col = "red",
     pch = 16,
     xlab = "Factors",
     ylab = "Variance",
     main = paste0("Scree plot of PEER on BLUE, ", ver, ", ", ref))
points(x = 2:25, y = variacne.from.second, pch = 16)
dev.off()


