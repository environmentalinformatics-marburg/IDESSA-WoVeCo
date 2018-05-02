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
library(Rsenal)
library(foreach)


# Additional settings ----------------------------------------------------------
rasterOptions(tmpdir = path_temp)

################################################################################
# Process Sentinel-2 dataset-------------------------------------------------
################################################################################
drs_S2<-dir(path = path_raster, 
            pattern = "S2A_USER_MTD_SAFL2A_PDMC.*data", full.names = TRUE)

sen2<- foreach(i = drs_S2, season = c("summer", "autumn", "winter")) %do% {
  fls <- list.files(i, pattern = ".img$", full.names = TRUE)
  bands <- paste0("B", c(2:8, 11, 12))
  
  # band layers
  bands_id <- sapply(bands, function(i) grep(i, fls))
  bands_id <- unlist(bands_id)
  fls_bands <- fls[bands_id]
  rst_bands <- stack(fls_bands)
  
  # ndvi layers
  ndvi_id <- grep("ndvi", fls)
  fls_ndvi <- fls[ndvi_id]
  rst_ndvi <- stack(fls_ndvi)
  
  # pca layers
  pca_id <- grep("pca", fls)
  fls_pca <- fls[pca_id]
  rst_pca <- stack(fls_pca)
  
  # stack all layers
  rst_all <- stack(rst_bands, rst_ndvi, rst_pca)
  names(rst_all) <- paste(names(rst_all), season, sep = "_")
  
  return(rst_all)
}

sen2 <- stack(sen2)

# Process Sentinel-1 dataset names ---------------------------------------------
fls <- list.files(path = paste0(path_raster, "S1A_IW_GRDH_1SDV"), 
                  pattern = ".img$", full.names = TRUE)
sen1 <- stack(fls)


# Create combined Sentinel 1/2 data stack --------------------------------------
sen <- stack(sen1, sen2)

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

#Create Map
caretpred <- predict(sen,caretrf)
writeRaster(caretpred,paste0(path_model,"caretprediction.tif"))

################################################################################

