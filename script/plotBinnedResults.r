args <- commandArgs(TRUE)

binResult <- read.delim("args[1]", header=FALSE)

colors <- c("red")	# default color for fpkm-binned / replication_timing-binned / 1Mb-binned results is red

plot(profile$V1, profile$V2, main=type, ylab="Translocation/Mb/sample ", pch=19, col=colors)