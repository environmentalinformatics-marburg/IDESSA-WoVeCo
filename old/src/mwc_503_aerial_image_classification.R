# Prepare model dataset (Sentinel and aerial images)

# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "E:/analysis/theunis/"
} else {
  filepath_base <- "/media/memory01/casestudies/hmeyer/SentinelMolopo/GFODat/"
}

path_figures <- paste0(filepath_base, "figures/")
path_raster <- paste0(filepath_base, "raster/")
path_rdata <- paste0(filepath_base, "rdata/")
path_temp <- paste0(filepath_base, "temp/")
path_vector <- paste0(filepath_base, "vector/")
path_model <- paste0(filepath_base, "model/")
path_source <- paste0("/home/hmeyer/hmeyer/molopo_wc/")


# Libraries --------------------------------------------------------------------
library(doParallel)
library(foreach)
library(raster)
library(rgdal)
library(maptools)
library(mapview)
library(rgeos)
library(caret)
library(Rsenal)
library(doParallel)

source(paste0(path_source, "mwc_functions.R"))
source(paste0(path_source, "AIclassification_functions.R"))

# Additional settings ----------------------------------------------------------
rasterOptions(tmpdir = path_temp)

# Load aerial images and RGB model ---------------------------------------------
model <- get(load(paste0(path_model,"googleAImodel.RData")))
highres_files <- list.files(path_rdata,pattern = glob2rx("highres_values_*.rds"))


cl <- makeCluster(detectCores()-10, outfile = "debug.txt")
registerDoParallel(cl)  


highres_cl<- foreach(img = 1:length(highres_files),
                     .packages = c("caret","raster","doParallel","Rsenal","foreach")) %dopar% {
                       aerial_images <- readRDS(paste0(path_rdata,highres_files[img]))
                       classified_df <- data.frame(matrix(ncol=2))
                       names(classified_df) <- c("woody","reliability")
                       
                       for (i in 1:length(aerial_images)){
                         print(paste0("file ",img," raster ",i," in process..."))
                         predicted_AI <- predict_woody(aerial_images[[i]], biome=2, model)
                         ### calculate woody percentage and reliability
                         woody <- sum(values(predicted_AI[[1]])==1)/(sum(values(predicted_AI[[1]])==1)+
                                                                       sum(values(predicted_AI[[1]])==0))
                         reliability <- (sum(values(predicted_AI[[2]])<0.25)+
                                           sum(values(predicted_AI[[2]])>0.75))/ncell(predicted_AI[[2]])
                         
                         classified_df[i,] <- c(woody,reliability)
                         outraster<-stack(aerial_images[[i]],predicted_AI)
                         names(outraster) <- c("R","G","B","pred","prob")
                         writeRaster(outraster,paste0(path_raster,"/classified/classified_",
                                                         substr(highres_files[img],1,
                                                                nchar(highres_files[img])-4),
                                                         "_",i,".tif"))
                         
                       }
                       save(classified_df, file = paste0(path_rdata,"/classified/classified_",
                                                         substr(highres_files[img],1,
                                                                nchar(highres_files[img])-4),".RData"))
                       
                     }
stopCluster(cl)