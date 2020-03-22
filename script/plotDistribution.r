args <- commandArgs(TRUE)

distribution <- read.delim("args[1]", row.names=1)

par(mar = c(4, 2, 4, 0))

colors <- c("steelblue4", "skyblue2", "powderblue")

barplot(as.matrix(distribution), ylab = "Translocation/Mb/sample", beside=TRUE, col=colors)

legend("topleft", c("silent", "nonsilent", "intergenic"), bty="n", fill=colors)
