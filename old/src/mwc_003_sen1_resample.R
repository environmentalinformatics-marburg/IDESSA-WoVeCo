# Prepare combined stack of sentinel 1 and 2 data

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


# Libraries --------------------------------------------------------------------
library(doParallel)
library(raster)


# Parallelization --------------------------------------------------------------
cl <- makeCluster(3)
registerDoParallel(cl)


# Additional settings ----------------------------------------------------------
rasterOptions(tmpdir = path_temp)


# Process Sentinel-1 dataset names -------------------------------------------------
# list available sentinel-1 files
fls <- list.files(path = paste0(path_raster, "S1A_IW_GRDH_1SDV"), 
                  pattern = ".tif$", full.names = TRUE)

sen2 <- raster(paste0(path_raster, 
                      "S2A_USER_MTD_SAFL2A_PDMC_20160114_resampled.data/", 
                      "B2.img"))
# loop over files
sen1 <- foreach(i = fls, .packages = "raster") %dopar% {
  # import current file 
  rst <- raster(i)
  # if required, extend and/or crop raster
  img <- gsub(".tif", ".img", i)
  
  if (file.exists(img)) {
    raster(img)
  } else {
    resample(rst, sen2,  method="ngb", filename = img, 
             format = "HFA", overwrite = TRUE)
  }
}

# sen1 <-stack(sen1)
# sen <- stack(sen1, sen2)
# 
# save(sen, file = paste0(path_rdata, "senstack.RData"))

## deregister parallel backend
stopCluster(cl)
