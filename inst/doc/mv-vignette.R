## ---- fig.width=7, fig.height=7------------------------------------------
library(multiview)
fourier  <- read.table("../digits_data/mfeat-fou.txt", header=FALSE, sep="")
profcorr <- read.table("../digits_data/mfeat-fac.txt", header=FALSE, sep="")
pixels   <- read.table("../digits_data/mfeat-pix.txt", header=FALSE, sep="")
morpho   <- read.table("../digits_data/mfeat-mor.txt", header=FALSE, sep="")

classes <- as.vector(sapply(0:9, function(x) rep(x,200)))

projection <- mvmds(list(fourier, profcorr, pixels, morpho), k=2)

mypalette <- c("chartreuse", "blueviolet", "deeppink1", "cyan2", "black", "blue3", 
              "gold1", "seagreen1", "gray60", "red")
plot(projection[,1:2], col = mypalette[classes+1], 
    pch = as.character(classes), axes = FALSE,  xlab="", ylab="")

## ---- fig.width=7, fig.height=7------------------------------------------
library(multiview)
fourier  <- read.table("../digits_data/mfeat-fou.txt", header=FALSE, sep="")
profcorr <- read.table("../digits_data/mfeat-fac.txt", header=FALSE, sep="")
pixels   <- read.table("../digits_data/mfeat-pix.txt", header=FALSE, sep="")
morpho   <- read.table("../digits_data/mfeat-mor.txt", header=FALSE, sep="")

classes <- as.vector(sapply(0:9, function(x) rep(x,200)))

projection <- mvtsne(list(fourier, profcorr, pixels, morpho), k=2)

mypalette <- c("chartreuse", "blueviolet", "deeppink1", "cyan2", "black", "blue3", 
              "gold1", "seagreen1", "gray60", "red")

plot(projection$embedding, col = mypalette[classes+1], 
    pch = as.character(classes), axes = FALSE,  xlab="", ylab="")

## ---- fig.width=7, fig.height=7------------------------------------------
library(multiview)
fourier  <- read.table("../digits_data/mfeat-fou.txt", header=FALSE, sep="")
profcorr <- read.table("../digits_data/mfeat-fac.txt", header=FALSE, sep="")
pixels   <- read.table("../digits_data/mfeat-pix.txt", header=FALSE, sep="")
morpho   <- read.table("../digits_data/mfeat-mor.txt", header=FALSE, sep="")

classes <- as.vector(sapply(0:9, function(x) rep(x,200)))

clust   <- mvsc(list(fourier, profcorr, pixels, morpho), k=10)

# $clustering member has the clustering assignment vector
knitr::kable(table(classes, clust$clustering))

## ---- fig.width=7, fig.height=7------------------------------------------
clust   <- mvsc(list(fourier, profcorr, pixels, morpho), k=10, neighbours=2)

# $clustering member has the clustering assignment vector
knitr::kable(table(classes, clust$clustering))

