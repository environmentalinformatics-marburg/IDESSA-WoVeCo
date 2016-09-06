# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "F:/GFO/BushEncroachment/"
} else {
  filepath_base <- "/media/TOSHIBA\ EXT//GFO/BushEncroachment/"
}

path_figures <- paste0(filepath_base, "figures/")
path_raster <- paste0(filepath_base, "raster/")
path_rdata <- paste0(filepath_base, "rdata/")
path_temp <- paste0(filepath_base, "temp/")
path_vector <- paste0(filepath_base, "vector/")
path_model <- paste0(filepath_base, "model/")
path_scripts <- paste0(filepath_base, "scripts/molopo_wc/src/")


# Libraries --------------------------------------------------------------------
library(gpm)
library(raster)
library(ggplot2)
library(rgdal)
library(foreach)


# Additional settings ----------------------------------------------------------
rasterOptions(tmpdir = path_temp)



################################################################################
###Prediction Caret
################################################################################
#Load caret rf model
caretrf <- get(load(paste0(path_model,"model_rf_caret.RData")))

#validation using independent test Data
testData<- get(load(paste0(path_model,"testdata_caret.RData")))
testpredict <- predict(caretrf,testData)
kstat(testpredict,testData$class)

#Plot Model characteristics

png(paste0(path_figures,"caret_varImp_all.png"))
plot(varImp(caretrf))
dev.off()
png(paste0(path_figures,"caret_varImp_10.png"))
plot(varImp(caretrf),10)
dev.off()

################################################################################
###Prediction GPM
################################################################################

GPMrf <- get(load(paste0(path_model,"models_rf_2016-08-25.RData")))
###Validation
 var_imp <- compVarImp(GPMrf)
 var_imp_scale <- compVarImp(GPMrf,scale=TRUE)
#Rewrite Sentinel-1 scene names to shorter derivative
 var_imp[[1]]$VARIABLE <- as.character(var_imp[[1]]$VARIABLE)
  var_imp[[1]]$VARIABLE[[which(var_imp[[1]]$VARIABLE=="S1A_IW_GRDH_1SDV_20160110_VH_SIGMA0")]] <-'S1A_summer'
 var_imp[[1]]$VARIABLE[[which(var_imp[[1]]$VARIABLE=="S1A_IW_GRDH_1SDV_20160720_VH_SIGMA0")]] <-'S1A_winter'
 var_imp[[1]]$VARIABLE[[which(var_imp[[1]]$VARIABLE=="S1A_IW_GRDH_1SDV_20160403_VH_SIGMA0")]] <-'S1A_autumn'
  var_imp[[1]]$VARIABLE <- as.factor(var_imp[[1]]$VARIABLE)
 
 #Plot of Variable Importance
 png(paste0(path_figures,"gpm_varImp_all600.png"),
     width = 15, height = 20, # Graph size
     units = "cm",            # unit sice for graph
     res = 600)               # Graph resolution
 plotVarImp(var_imp)
 dev.off()

 #Heatmap 
png(paste0(path_figures,"gpm_varimp_heat.png"))
plotVarImpHeatmap(var_imp_scale, xlab = "Class", ylab = "Band")
dev.off()


# Kappa Response Plot
tstat <- compContTests(GPMrf, mean = TRUE)
tstat
tstat_mean <- merge(tstat[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")
tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]
# ggplot(data = tstat[[2]][[1]], aes(y = Response, x = Kappa)) + geom_point() + geom_smooth()
png(paste0(path_figures,"gpm_varRes_all600.png"),
    width = 10, height = 15, # Graph size
    units = "cm",            # unit sice for graph
    res = 600)               # Graph resolution
ggplot(data = tstat[[2]][[1]], aes(x = Response, y = Kappa))+geom_boxplot()
dev.off()


################################Image reclassification################################################

#read raster
GPMpred<-raster(paste0(path_model,"GPMprediction.tif"))


#levels(GPMpred) <- data.frame("ID"=1:length(GPMrf[[1]][[1]]$model$levels),"Class"=GPMrf[[1]][[1]]$model$levels)
#GPMpred is the predicted raster and GPMrf the gpm based model


#Write a reclassification matrix
m<-c(4.9,8.1,1,
     9.9,11.1,1,
     8.9,9.1,5)
reclasstable<-matrix(m,ncol=3, byrow = TRUE)
# Perform reclassification based on reclassification matrix
GPM_pred_reclass<-reclassify(GPMpred, rcl = reclasstable)
#Assign class names to reclassified image
levels(GPM_pred_reclass)<-data.frame("ID"=c(1,2,3,4,5),
                                     "Class"=c("Bare","Irrigation",
                                               "Low bush encroachment",
                                               "Moderate bush encroachment",
                                               "Severe bush encroachment"))

#Export raster object as tif
writeRaster(GPM_pred_reclass,paste0(path_model,"GPMprediction_reclass.tif"))

#################################Creating Map layout################################################
#Read ancillairy vector data for map overlay
ogrListLayers(paste0(path_vector, "Municipalities.shp"))
localmun<- readOGR(paste0(path_vector,"Municipalities.shp"),"Municipalities")
proj4string(localmun)

#Create vector string defining colours used for each class
colourcodes<-(c("rosybrown1",       #Bare Ground
                "paleturquoise4",   #Irrigation
                "cornsilk",         #Low Bush Encraochment
                "gold3",            #Moderate Bush Encraochment
                "saddlebrown"       #Severe Bush Encraochment
                ))


yat = seq(extent(GPM_pred_reclass)@ymin, 
          extent(GPM_pred_reclass)@ymax, length.out = 5)
xat = seq(extent(GPM_pred_reclass)@xmin, 
          extent(GPM_pred_reclass)@xmax, length.out = 5)

png(paste0(path_figures,"GPM_pred_reclass_Map_600.png"),
                  width = 30, height = 30, # Graph size
                  units = "cm",            # unit sice for graph
                  res = 600)               # Graph resolution
spplot(GPM_pred_reclass, col.regions = colourcodes, maxpixels= 1000000,
       panel = function(...){
         panel.levelplot(...)
         panel.abline(h = yat, v = xat, col = "grey0", lwd = 0.8, lty = 3)
         sp.polygons(localmun, lwd = 2, col="red1")
       },
       scales = list(x = list(at = xat),
                     y = list(at = yat)))
dev.off()




