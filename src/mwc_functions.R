# Extract samples of high resolution pixels based on low resolution raster overlay
sample_xres <- function(lowres, highres_files, highres_prj = NULL,
                        path_rdata, path_temp, n = 10000){
  library(raster)
  library(rgdal)
  library(sp)
  
  rasterOptions(tmpdir = path_temp)
  
  
  # Get overall extent of all highres images -----------------------------------
  highres_extent <- SpatialPolygonsDataFrame(
    as(extent(raster(highres_files[[1]])), "SpatialPolygons"),
    data = data.frame(NAME = highres_files[[1]]))
  
  for(f in highres_files[2:length(highres_files)]){
    highres_extent <- bind(highres_extent, SpatialPolygonsDataFrame(
      as(extent(raster(f)), "SpatialPolygons"),
      data = data.frame(NAME = f)))
  }
  
  if(!is.null(highres_prj)){
    projection(highres_extent) <- highres_prj
  }
  highres_extent <- spTransform(highres_extent, CRS(projection(lowres)))

  
  # Intersect low resolution and high resolution images ------------------------
  lowres_extent <- as(extent(lowres), "SpatialPolygons")
  projection(lowres_extent) <- projection(lowres)
  
  mask <- intersect(lowres_extent, highres_extent)
  mask <- disaggregate(mask)
  
  
  # lowres_pixels <- lapply(seq(length(mask)), function(p){
  lowres_pixels <- lapply(seq(2), function(p){
    lowres_pixels <- crop(lowres[[1]], mask[p, ], cellnumbers = TRUE)
    lowres_pixels <-  rasterToPolygons(lowres_pixels, fun=NULL, n=4, na.rm=TRUE, 
                                       digits=12, dissolve=FALSE)
    names(lowres_pixels) <- as.character(mask[p, ]@data[1,1])
    return(lowres_pixels)
  })
  saveRDS(lowres_pixels, file = paste0(path_rdata, "lowres_pixels.rds"))
  # lowres_pixels <- readRDS(paste0(path_rdata, "lowres_pixels.rds"))
  
  
  # Select n low resolution pixel ids by chance --------------------------------
  n <- 100
  n_sample <- round(n / length(lowres_pixels), 0)
  
  samples <- lapply(seq(length(lowres_pixels)), function(s){
    total_pixels <- length(lowres_pixels[[s]])
    set.seed(s)
    lowres_samples <- sample(total_pixels, n_sample)
    lowres_pixels[[s]][lowres_samples,]
  })
  names(samples) <- sapply(lowres_pixels, names)
  
  
  
  # Extract information from high resolution images for each sampled pixel -----
  highres_info <- lapply(seq(samples), function(s){
    act_highres <- raster(highres_files[grep(names(samples[s]), highres_files)])
    act_lowres_pixels <- spTransform(samples[[s]], CRS(projection(act_highres)))
    ext <- extract(act_highres, act_lowres_pixels)
  })
  names(highres_info) <- names(samples)
  saveRDS(highres_info, file = paste0(path_rdata, "highres_info.rds"))
  # highres_info <- readRDS(paste0(path_rdata, "highres_info.rds"))
  

  # Extract information from low resolution images for each sampled pixel ------
  lowres_info <- lapply(seq(samples), function(s){
    act_lowres_pixels <- spTransform(samples[[s]], CRS(projection(act_highres)))
    lowres_crop <- crop(lowres, act_lowres_pixels)
    ext <- extract(lowres_crop, act_lowres_pixels, cellnumbers = TRUE)
  })
  names(lowres_info) <- names(samples)
  saveRDS(lowres_info, file = paste0(path_rdata, "lowres_info.rds"))
  # lowres_info <- readRDS(paste0(path_rdata, "lowres_info.rds"))
  
  
  # Combine high and low resolution information --------------------------------
  info <- list(LowRes_Info = lowres_info, HighRes_Ifo = highres_info,
               LowRes_Pixels = lowres_pixels)
}
  

  


  
  

























#Prepare combined stack of sentinel 1 and 2 data
stack_sen <- function(path_raster, path_temp){
  library(raster)
  library(rgdal)
  
  
  # Additional settings ----------------------------------------------------------
  rasterOptions(tmpdir = path_temp)
  
  
  # Process Sentinel-2 dataset names -------------------------------------------------
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
}

