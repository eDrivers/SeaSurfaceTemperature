# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                     LIBRARIES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(magrittr)
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   DOWNLOAD DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The data used to characterize comes from DFO and cannot be shared
# For more information read the repo's README.md document.

# Output location for downloaded data
output <- './Data/RawData'

# Data will need to be archived to Zenodo with restricted access and downloaded
# using an access token.
# Eventually it would ideally be part of the SLGO web portal

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       IMPORT AND FORMAT ANOMALIES 2013-2017
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------------------------------------- #
#      Sea surface temperature anomalies 2013-2017      #
# ----------------------------------------------------- #
# File name
fileName <- 'SSTanomaly2013-2017'

# Unzip
unzip(zipfile = paste0(output, '/', fileName, '.zip'),
      exdir = paste0(output, '/', fileName))

# Get folder names
ys <- dir(paste0(output, '/', fileName)) %>%
      as.numeric() %>%
      na.omit()

ms <- dir(paste0(output,  '/', fileName, '/', ys[1]))

# Create array to store monthly data per year
sst <- vector('list',12)

# Replicate for the number of years
sst <- replicate(length(ys), sst, simplify = F)

# Load files
for(i in 1:length(ys)) {
  for(j in 1:length(ms)) {
    filePath <- dir(paste0(output,  '/', fileName, '/', ys[i], '/', ms[j], '/'),
                    pattern = 'anomaly.dat',
                    full.names = T)
    sst[[i]][[j]] <- read.table(filePath)
  }
}

# ----------------------------------------------------- #
#                     CHECK OUTLIERS                    #
# ----------------------------------------------------- #
# Visualize outliers
par(mfrow = c(3,2))
for(i in 1:5) {
  graphicsutils::plot0(xlim = c(0,12), ylim = c(-20, 20))
  mtext(ys[i], 3, font = 2)
  for(j in 1:12) {
    boxplot(sst[[i]][[j]]$V3, add = T, at = j, frame = F, range = 3)
  }
  axis(1, at = 1:12, labels = 1:12)
}

# Identify and modify outliers in the dataset
for(i in 1:length(ys)) {
  for(j in 1:length(ms)) {
    out <- boxplot.stats(sst[[i]][[j]]$V3, coef = 3)$stat[c(1,5)]
    cap <- quantile(sst[[i]][[j]]$V3, probs=c(.05, .95), na.rm = T)
    sst[[i]][[j]]$V3[sst[[i]][[j]]$V3 < out[1]] <- cap[1]
    sst[[i]][[j]]$V3[sst[[i]][[j]]$V3 > out[2]] <- cap[2]
  }
}

# Visualize again
par(mfrow = c(3,2))
for(i in 1:5) {
  graphicsutils::plot0(xlim = c(0,12), ylim = c(-20, 20))
  mtext(ys[i], 3, font = 2)
  for(j in 1:12) {
    boxplot(sst[[i]][[j]]$V3, add = T, at = j, frame = F, range = 3)
  }
  axis(1, at = 1:12, labels = 1:12)
}


# ----------------------------------------------------- #
#                    SPATIAL OBJECTS                    #
# ----------------------------------------------------- #
# Data.frame with all unique coordinates
coords <- data.frame(V1 = numeric(), V2 = numeric())
for(i in 1:2) {
  for(j in 1:length(ms)) {
   coords <- rbind(coords, sst[[i]][[j]][, 1:2])
  }
}
coords <- unique(coords)

# Join datasets per year
sst2 <- vector('list', length(ys))
for(i in 1:length(ys)) {
  sst2[[i]] <- coords
  for(j in 1:12) {
    sst2[[i]] <- dplyr::left_join(sst2[[i]], sst[[i]][[j]],
                                  by = c('V1','V2'))
  }
  # Change column names
  colnames(sst2[[i]]) <- c('lat','lon', ms)
}

# Replace and remove objects to save memory
sst <- sst2
rm(sst2)

# Select only months from May to November due to biases in anomalies measurements
for(i in 1:length(ys)) {
  sst[[i]] <- sst[[i]][, -c(3,4,5,6,14)]
}

# Transform as spatial object and keep coordinats
sst <- lapply(sst, function(x) st_as_sf(x, coords = c("lon", "lat"), crs = 4326, remove = FALSE))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        IMPORT AND FORMAT REFERENCE 1980-2010
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------------------------------------- #
# Reference data (1980-2010) for normalization #
# ----------------------------------------------------- #
# Unzip file
# File name
fileName <- 'Clim1980-2010-Var'

# Unzip
unzip(zipfile = paste0(output, '/', fileName, '.zip'),
      exdir = paste0(output, '/', fileName))

# Get data names and import
ref <- dir(paste0(output, '/', fileName), full.names = T) %>%
       lapply(read.table)


# ----------------------------------------------------- #
#                    SPATIAL OBJECTS                    #
# ----------------------------------------------------- #
# Get all unique coordinates
xy <- matrix(ncol = 2)
for(i in 1:length(ref)) xy <- rbind(xy, ref[[i]][, 1:2])
xy <- unique(xy) %>% na.omit()

# Join all data in a single dataset
for(i in 1:length(ref)) xy <- dplyr::left_join(xy, ref[[i]], by = c('V1' = 'V1', 'V2' = 'V2'))
colnames(xy) <- c('lat','lon','05','06','07','08','09','10','11')

# Create spatial object
ref <- st_as_sf(xy, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
rm(xy)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        EXTENT REFERENCE vs CURRENT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify observations in current datasets (2013-2017) that overlap with reference (1980-2010)
ext <- st_bbox(ref)
id <- cbind(c(ext$xmin, ext$xmin, ext$xmax, ext$xmax, ext$xmin),
      c(ext$ymin, ext$ymax, ext$ymax, ext$ymin, ext$ymin)) %>%
      list() %>%
      st_polygon() %>%
      st_intersects(., sst[[1]]) %>%
      unlist()

# Subset current dataset
for(i in 1:length(ys)) sst[[i]] <- sst[[i]][id, ]

# Join to locations in sst
ref <- sst[[1]][, 1:2, drop = T] %>%
       dplyr::left_join(ref, by = c('lon' = 'lon', 'lat' = 'lat')) %>%
       st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          GAP-FILLING MISSING DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gap-fill missing data with mean of neighbouring values for each month
# Takes > 10000s to run, i.e. ~ 3h
for(i in 3:(ncol(ref)-1)) {
  idNA <- which(is.na(ref[, i]))
  for(j in idNA) {
    ref[j, i] <- st_distance(ref[j, ], ref) %>% # Distance between point j and all other points
                 t() %>%
                 as.data.frame() %>%
                 rename(dist = '.') %>%
                 mutate(id = 1:nrow(ref), # add ID
                        val = ref[, i, drop = T]) %>% # get normalizing values
                 arrange(dist) %>% # sort table min to max distance
                 .[-1, ] %>% # Remove point j, i.e. distance with itself = 0
                 na.omit() %>% # remove NAs
                 .[1:4, ] %>% # Select 4 closest points for gap-filling
                 mutate(wt = (min(dist) / dist) / sum(min(dist) / dist)) %>% # weighted to use to measure weighted average as a function of distance, with closest points having the highest impact on the average
                 mutate(wmean = wt * val) %>% # Calculate weighted average contribution from each point
                 summarise(wmean = sum(wmean)) # Calculate weighted average
    }
}
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                             NORMALIZE ANOMALIES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Normalize anomalies
for(i in 1:length(ys)) {
  for(j in 3:(ncol(ref)-1)) {
    sst[[i]][, j] <- sst[[i]][, j, drop = T] / ref[, j, drop = T]
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   SPATIAL
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Transform projection
sst <- lapply(sst, function(x) st_transform(x, crs = 32198))

# Remove coordinates
sst <- lapply(sst, function(x) x[, 3:ncol(x)])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  EXPORT DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Export object as .RData
save(sst, file = './data/rawData/sst.RData')
