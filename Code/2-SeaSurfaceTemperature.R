# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  LIBRARIES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(magrittr)
library(tidyverse)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('./data/rawData/sst.RData')
ys <- 2013:2017
nC <- colnames(sst[[1]])[-ncol(sst[[1]])]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   ANOMALIES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Transform values between -.5 and +.5 to NA
for(i in 1:length(ys)) {
  sst[[i]][, nC] <- apply(sst[[i]][, nC, drop = T], 2, function(x) ifelse(x > -.5 & x < .5, NA, x))
}

# Positive & negative anomalies
pos <- neg <- sst
for(i in 1:length(ys)) {
  pos[[i]][, nC] <- apply(pos[[i]][, nC, drop = T], 2, function(x) ifelse(x < 0, NA, x))
  neg[[i]][, nC] <- apply(neg[[i]][, nC, drop = T], 2, function(x) ifelse(x > 0, NA, x))
}

# Transform negative anomalies as positive values
for(i in 1:length(ys)) neg[[i]][, nC] <- apply(neg[[i]][, nC, drop = T], 2, abs)

# Sum rows to generate annual anomaly index
pos <- lapply(pos, function(x) x %>% mutate(posSST = rowSums(x[,nC,drop=T], na.rm = T)))
neg <- lapply(neg, function(x) x %>% mutate(negSST = rowSums(x[,nC,drop=T], na.rm = T)))

# Single dataset for all annual anomaly indices
posSST <- cbind(pos[[1]][, 'posSST'], pos[[2]]$posSST, pos[[3]]$posSST, pos[[4]]$posSST, pos[[5]]$posSST)
negSST <- cbind(neg[[1]][, 'negSST'], neg[[2]]$negSST, neg[[3]]$negSST, neg[[4]]$negSST, neg[[5]]$negSST)
colnames(posSST)[1:5] <- colnames(negSST)[1:5] <- as.character(2013:2017)

# Mean of all annual indices as the anomaly index
posSST$PositiveSST <-rowMeans(posSST[, 1:5, drop = T])
negSST$NegativeSST <-rowMeans(negSST[, 1:5, drop = T])

# Remove 0 values
posSST <- posSST[posSST$PositiveSST > 0, ]
negSST <- negSST[negSST$NegativeSST > 0, ]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  EXPORT DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Change names for uniformity
PositiveSST <- posSST
NegativeSST <- negSST

# Export object as .RData
save(PositiveSST, file = './Data/Driver/PositiveSST.RData')
save(NegativeSST, file = './Data/Driver/NegativeSST.RData')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 VISUALIZE DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png('./Figures/PositiveSST.png', width = 1280, height = 1000, res = 200, pointsize = 6)
plot(posSST[, 'PositiveSST'], pch = 20, cex = .5)
dev.off()

png('./Figures/NegativeSST.png', width = 1280, height = 1000, res = 200, pointsize = 6)
plot(negSST[, 'NegativeSST'], pch = 20, cex = .5)
dev.off()
