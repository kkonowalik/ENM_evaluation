# ENM evaluation
Some of the R code that was used in the publication: 
Konowalik, K., Nosol, A. Evaluation metrics and validation of presence-only species distribution models based on distributional maps with varying coverage. Sci Rep 11, 1482 (2021). https://doi.org/10.1038/s41598-020-80062-1

Hi,
I was going to publish my R code that was used in this publication. I hope to do so but in the meantime I will also publish one that is prbably the most useful. One of the findings of that paper was that you can use AUC and MAE to select which models have the best scores for your data (see figure 6 in the publication). It is one of the ways that may be used - I'm not an R expert but it works quite easily! Let's start.

You need three input files:

1) raster
>library(raster)
>#"mm_TESTING01" is your raster file (model), you can use a result from maxent or any other algorithm, even if you do it outside R
>#for example
>mm_TESTING01 <- raster("M:/GIS/models/maxent_model.tif")

2) presence points
'points_pres <- read.csv("D:/rs_pres.csv",header=TRUE)
points_pres1 <- points_pres[,2:3]
#column 2 and 3 are latitude and longitude'

3) absence points
points_abs <- read.csv("D:/rot1/dysk_google/rs_abs.csv",header=TRUE)
points_abs1 <- points_abs[,2:3]

#if your absence or presence points are in wgs84 you can assign a crs, sometimes it is necessary to have a crs assigned
points_abs1 <- SpatialPointsDataFrame(coords = points_abs1, data = points_abs1, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

Then it goes like this:
library("Metrics")
#you extract the values of raster at you presence and absence points
points_pres2 <- raster::extract(mm_TESTING01, points_pres1, cellnumbers=TRUE)
points_abs2 <- raster::extract(mm_TESTING01, points_abs1, cellnumbers=TRUE)

testECE <- c(rep(1, nrow(points_pres2)), rep(0, nrow(points_abs2)))
testECE <- data.frame(cbind(testECE, rbind(points_pres2, points_abs2)))
testECE <- na.omit(testECE)
colnames(testECE) <- c('testECE', 'cells', 'layer')
testECE <- subset(testECE, select = c(testECE, layer))
bias <- bias(testECE$testECE, testECE$layer)
MAE <- mae(testECE$testECE, testECE$layer)


In case you will use it I will appreciate if you can provide a citation to the paper - the details are in the beginning and the linkis here: https://rdcu.be/cjB5b
Thank's!
Best regards,
Kamil
