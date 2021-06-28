# ENM evaluation and chosing the algorithms
Some of the R code that was used in the publication: 
>Konowalik, K., Nosol, A. Evaluation metrics and validation of presence-only species distribution models based on distributional maps with varying coverage. Sci Rep 11, 1482 (2021). https://doi.org/10.1038/s41598-020-80062-1

I hope to publish the whole R code and the data that was used for this publication but in the meantime, I will also publish one piece that is probably the most useful. One of the findings of that paper was that you can use simultaneously AUC and MAE to compare the modeling results and select which models have the best scores for your data (see figure 6 in the publication). It is one of the ways how to apply it - probably there is a nicer way to do it, but like this, it works quite easily. 
You need three input files:

1) raster
```
library(raster)
#"mm_TESTING01" is your raster file (model), you can use a result from maxent or any other algorithm, even if you do it outside R
#for example
mm_TESTING01 <- raster("M:/GIS/models/maxent_model.tif")
```
2) presence points
```
points_pres <- read.csv("D:/rs_pres.csv",header=TRUE)
points_pres1 <- points_pres[,2:3]
#column 2 and 3 are latitude and longitude
```
3) absence points
(if you don't have true absences you can also use your background points)
```
points_abs <- read.csv("D:/rot1/dysk_google/rs_abs.csv",header=TRUE)
points_abs1 <- points_abs[,2:3]
#if your absence or presence points are in wgs84 you can assign a crs, sometimes it is necessary to have a crs assigned
points_abs1 <- SpatialPointsDataFrame(coords = points_abs1, data = points_abs1, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
```
We are going to use functions from package "Metrics". First we load the package, assign raster values to points, make a table (it is easier to keep also the cell numbers), and finally compute metrics AUC, MAE, and bias. It goes like this:
```
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
AUC <- auc(testECE$testECE, testECE$layer)
```

In case you will use this method I will appreciate if you can provide a citation to the paper - the details are in the beginning and the link is here: https://rdcu.be/cjB5b

Thank's and good luck with modeling!
