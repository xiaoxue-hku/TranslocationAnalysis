args <- commandArgs(TRUE)

profile <- read.delim("args[1]", header=FALSE)

colors <- c("red")	# default color for profiles is red

x <- seq(1, 90, 1)

plot(x, profile, main=cohort, ylab="Translocation/Mb/sample ", type="l", col=colors)	