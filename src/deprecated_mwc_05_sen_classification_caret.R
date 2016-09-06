### environmental settings -----

## set working directory
if(Sys.info()["sysname"] == "Windows") {
  filepath_base <- "F:/GFO/BushEncroachment/"
} else {
  filepath_base <- "/media/TOSHIBA EXT/GFO/BushEncroachment/"
}

path_figures <- paste0(filepath_base, "figures/")
path_model <- paste0(filepath_base, "model/")
path_raster <- paste0(filepath_base, "raster/")
path_rdata <- paste0(filepath_base, "rdata/")
path_temp <- paste0(filepath_base, "temp/")
path_vector <- paste0(filepath_base, "vector/")

## required packages
library(raster)
library(lattice)
library(ggplot2)
library(caret)
library(randomForest)
library(rgdal)
library(leaflet)
library(mapview)
library(gpm)
#Load Rasterstack

sen<- get(load(paste0(path_raster, "senstack.RData")))

ogrListLayers(paste0(path_vector, "Training_Areas_20160823.shp"))
plots<-readOGR(dsn = paste0(path_vector, "Training_Areas_20160823.shp","Training_Areas_20160823"))

proj4string(plots)
proj4string(sen)

plots<- spTransform(plots, CRS("+proj=utm +zone=34 +south +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
proj4string(plots)

#Extract raster values using training plots
Plots_exstract<- extract(sen, plots, df=TRUE, na.rm=TRUE)
head(Plots_exstract)
print((plots))


#Merge extracted data to plot data
ID <- as.integer(rownames(plots@data)) + 1 # Add ID
plots$ID <- ID
plots@data <- merge(plots@data, Plots_exstract, by = "ID", all = TRUE)

#Alternative option to merge to datasets
extracteddata$class<-plots$CLASS[extracteddata$ID]
head(extracteddata)

# Save extracted data to disk
save(extracteddata, file = paste0(path_model, "extracteddata.RData"))
extracteddata<-get(load(paste0(path_model, "extracteddata.RData")))


################################Caret Procedure################################
#partition data into training and test data
dim(extracteddata)
set.seed(25)
trainid<-createDataPartition(extracteddata$class,
                             p=0.5,       # partitioning threshold
                             list= FALSE, # 
                             times = 1)   #
#Split data into training set and testing dataset
traindata<-extracteddata[trainid,]
testdata<-extracteddata[-trainid,]
#Specify predictor and response variables
predictors<-traindata[,2:52]
response<-traindata[,53]

#Train Random Forest
RFModel<-train(predictors,response,
               method = "rf",  # Classifier type: rf for Random Forest
               tuneLength = 3, # Control tuning parameter
               trControl = trainControl(method = "cv")) # cv: cross validation

#############################GPM Procedure#################################

