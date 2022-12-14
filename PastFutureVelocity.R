#### Climate Velocity ####

#This calculates the velocity of climate change for past and future downscaled rasters

cluster <- TRUE

if(cluster == TRUE){
  library(gstat)
  library(pbapply, lib.loc="/nobackup/beegfs/home/ISAD/hw656/R/x86_64-pc-linux-gnu-library/3.5/sp/libs")
  library(pbmcapply, lib.loc="/nobackup/beegfs/home/ISAD/hw656/R/x86_64-pc-linux-gnu-library/3.5/sp/libs")
  library(raster, lib.loc="/nobackup/beegfs/home/ISAD/hw656/R/x86_64-pc-linux-gnu-library/3.5/sp/libs")
  library(rgeos, lib.loc="/nobackup/beegfs/home/ISAD/hw656/R/x86_64-pc-linux-gnu-library/3.5/sp/libs")
  library(data.table, lib.loc="/nobackup/beegfs/home/ISAD/hw656/R/x86_64-pc-linux-gnu-library/3.5/sp/libs")
  library(dismo)
  library(stringr, lib.loc="/nobackup/beegfs/home/ISAD/hw656/R/x86_64-pc-linux-gnu-library/3.5/sp/libs")
  library(CircStats, lib.loc="/nobackup/beegfs/home/ISAD/hw656/R/x86_64-pc-linux-gnu-library/3.5/sp/libs")

  DataFP <- "/nobackup/beegfs/workspace/hw656/Data/"
  ResultsFP <- "/nobackup/beegfs/workspace/hw656/Projects/PaleoClimateVelocity/"
  
  ncores=16
  
} else {
  library(gstat)
  library(pbapply)
  library(pbmcapply)
  library(raster)
  library(rgeos)
  library(data.table)
  library(dismo)
  library(stringr)
  library(CircStats)
  
  DataFP <- "/Users/hannahwauchope/Dropbox/Work/Data/"
  ResultsFP <- "/Users/hannahwauchope/Dropbox/Work/Projects/PaleoClimateVelocity/"
  ncores=6
}

source("~/Documents/GitHub/ClimateVelocity/ClimateVelocity.r")

#### Functions ####
Past <- "DONE"

if(Past!="DONE"){
  TimeSteps <- rev(seq(1, 210, 1))
  #TimeSteps <- c(210, 209, 208, 1)
  pbmclapply(TimeSteps, function(TS){
    print(TS)
    if(TS==1){
      Period <- list.files(paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/BioClim"), pattern=paste0("Downscaled_BioVars_.{2,2}BP.gri"), full.names=TRUE)
    } else {
      Period <- list.files(paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/BioClim"), pattern=paste0("Downscaled_BioVars_", TS-1, ".{2,2}BP.gri"), full.names=TRUE)
      Period <- Period[Period != paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/BioClim/Downscaled_BioVars_", TS-1, "00BP.gri")]
    }
    Period <- c(Period, paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/BioClim/Downscaled_BioVars_", TS, "00BP.gri"))
    Period <- rev(Period)
    YearNames <- c("100BP", "80BP", "60BP", "40BP", "20BP")
    YearVals <- c(0, 20, 40, 60, 80)
    
    PeriodBiovars <- lapply(c(1:19), function(x){
      PeriodStack <- stack(lapply(Period, function(y) raster(y, band=x)))
      PeriodStack
      names(PeriodStack) <- YearNames
      return(PeriodStack)
    })
    
    #With full raster
    ClimVelBioVars <- stack(pblapply(PeriodBiovars, function(x){
      TempVel <- tempTrendPerYear(x, 5, YearVals)
      SpatVel <- spatGrad(x[[1]], th = 0.01)
      ClimVel <- gVoCC(TempVel, SpatVel)
      return(ClimVel)
    }))
    
    writeRaster(ClimVelBioVars, paste0(ResultsFP, "21000to0BP_100yrVelocity_Bioclim/BioclimVelocity_", TS, "00.grd"))
  }, mc.cores=ncores)
}

#### Calculate Velocity (future) ####

ClimMods <- c("ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CanESM5", "CanESM5-CanOE", "CMCC-ESM2", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg", "EC-Earth3-Veg-LR", "GISS-E2-1-G", "GISS-E2-1-H", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL") #"FIO-ESM-2-0", "GFDL-ESM4", "HadGEM3-GC31-LL", 
ssps <- c("126", "245", "370", "585")
timesteps <- c("2021-2040", "2041-2060", "2061-2080")
var <- "bioc" #c("tmin", "tmax", "tavg", "prec", "bioc")
res <- 10 #c(10, 5, 2.5)
YearVals <- c(0, 20, 40)

ModSSPs <- expand.grid(ClimMods, ssps)

pbmclapply(c(1:nrow(ModSSPs)), function(ModSSP){
  Mod <- ModSSPs[ModSSP, 1]
  ssp <- ModSSPs[ModSSP, 2]
  print(paste0(Mod,"   ", ssp))
  
  dir.create(paste0(ResultsFP, "2020to2080_CMPI6_Velocity_Bioclim/", ssp, "/"), showWarnings = FALSE)
  
  Period <- list.files(paste0(DataFP, "CMIP6"), pattern=paste0("wc2.1_", res, "m_", var, "_", Mod, "_ssp", ssp, "_.{9,9}.tif"), full.names=TRUE, recursive=TRUE)
  Period <- sort(Period)
  
  PeriodBiovars <- lapply(c(1:19), function(x){
    PeriodStack <- stack(lapply(Period, function(y) raster(y, band=x)))
    names(PeriodStack) <- c(timesteps)
    return(PeriodStack)
  })
  
  #With full raster
  ClimVelBioVars <- stack(pblapply(PeriodBiovars, function(x){
    TempVel <- tempTrendPerYear(x, 3, YearVals)
    SpatVel <- spatGrad(x[[1]], th = 0.01)
    ClimVel <- gVoCC(TempVel, SpatVel)
    return(ClimVel)
  }))
  
  writeRaster(ClimVelBioVars, paste0(ResultsFP, "2020to2080_CMPI6_Velocity_Bioclim/", ssp, "/BioclimVelocity_", Mod, ".grd"))
}, mc.cores=ncores)

#### My own versions that didn't quite work ####

# #Temporal velocity
# time <- c(1:nlayers(PeriodBiovars[[1]]))
# time <- c(1970, 2020, 2040, 2060)
# fun <- function(x) { lm(x ~ time)$coefficients[2] }
# fun <- function(x) { if (is.na(x[1])){ NA } else {lm(x ~ time)$coefficients[2] }}
# 
# #This is mid way through being a thing, the coefficients part isn't right
# fun <- function(x) { if (is.na(x[1])){ NA } else {
#   rasmod <- lm(x ~ time)
#   return(list(lm(x ~ time)$coefficients[2], lm(x ~ time)$coefficients[2], lm(x ~ time)$coefficients[1]))}}
# 
# TempVelocity <- calc(PeriodBiovars[[1]], fun)
# 
# #Spatial velocity
# w <- matrix(1, 3, 3)  #c(0,1,0,1,0,1,0,1,0), nrow=3)
# x <- focal(PeriodBiovars[[1]][[1]], w, mean, na.rm=TRUE, NAonly=TRUE, pad=TRUE)
# plot(x)

