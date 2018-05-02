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
###Prediction GPM
################################################################################

GPMrf <- get(load(paste0(path_model,"models_rf_2016-08-25.RData")))
###Validation
 var_imp <- compVarImp(GPMrf)
 var_imp_scale <- compVarImp(GPMrf,scale=TRUE)

png(paste0(path_figures,"gpm_varImp_all.png"))
plotVarImp(var_imp)
dev.off()

png(paste0(path_figures,"gpm_varimp_heat.png"))
plotVarImpHeatmap(var_imp_scale, xlab = "Class", ylab = "Band")
dev.off()
# 
tstat <- compContTests(GPMrf, mean = TRUE)
tstat
# 
#tstat_mean <- merge(tstat[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")
# 
# tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]
# 
# ggplot(data = tstat_mean, aes(x = OCCURENCE, y = Kappa_mean)) + geom_point() + geom_smooth()
# 

##Create Map
GPMpred <- predict(sen,GPMrf[[1]][[1]]$model)
writeRaster(GPMpred,paste0(path_model,"GPMprediction.tif"))
