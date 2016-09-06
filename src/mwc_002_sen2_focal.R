# Check Sentinel-2 data and compute artifical images

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
library(mapview)
library(raster)


# Additional settings ----------------------------------------------------------
rasterOptions(tmpdir = path_temp)


# Import data ------------------------------------------------------------------
# list available files
fls <- list.files(path_raster, pattern = "^ndvi.img$", full.names = TRUE,
                  recursive = TRUE)

# matrix of weights
w <- matrix(c(1, 1, 1,
              1, 1, 1,
              1, 1, 1), nrow = 3, ncol = 3)

# calculate focal mean and standard deviation
lst_fcl <- lapply(fls, function(i) {
  
  # import image
  rst <- raster(i)
  
  # calculate focal mean and standard deviation
  fls_mn <- gsub("ndvi.img", "ndvi_mn.img", i)
  # rst_mn <- if (file.exists(fls_mn)) {
  #   raster(fls_mn) 
  # } else {
  focal(rst, w = w, fun = mean, na.rm = TRUE, pad = TRUE,
        filename = fls_mn, format = "HFA", overwrite = TRUE)
  # }
  
  fls_sd <- gsub("ndvi.img", "ndvi_sd.img", i)
  # rst_sd <- if (file.exists(fls_sd)) {
  # raster(fls_sd) 
  # } else {
  focal(rst, w = w, fun = sd, na.rm = TRUE, pad = TRUE,
        filename = fls_sd, format = "HFA", overwrite = TRUE)
  # }
  
  return(stack(rst_mn, rst_sd))
})


### display data -----
# mapview(lst_fcl[[1]])
