### I've not annotated this well, and it needs to be converted to terra (from raster)
#BUT this script is to calculate, for past and future velocities of change, where in the distribution of past velocities the future velocity falls (For each gridcell)

#### Velocity Distribution Rasters ####
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

VelNums <- c(210:1)

Vels <- lapply(VelNums, function(TS){
  stack(paste0(ResultsFP, "21000to0BP_100yrVelocity_Bioclim/BioclimVelocity_", TS, "00.grd"))
})

names(Vels) <- VelNums

ClimMods <- c("ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CanESM5", "CanESM5-CanOE", "CMCC-ESM2", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg", "EC-Earth3-Veg-LR", "GISS-E2-1-G", "GISS-E2-1-H", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL") #"FIO-ESM-2-0", "GFDL-ESM4", "HadGEM3-GC31-LL", 
ssps <- c("245", "370", "585") #"126", 

ModSSPs <- expand.grid(ClimMods, ssps)

pblapply(nrow(ModSSPs):1, function(ModSSP){
  Mod <- ModSSPs[ModSSP, 1]
  ssp <- ModSSPs[ModSSP, 2]
  print(paste0(Mod,"   ", ssp))
  if(file.exists(paste0(ResultsFP, "FuturevsPastDistribution/Logbook/", ssp, "_", Mod, "_Begin.csv"))){
    return(NULL)
  }
  write.csv(NULL, paste0(ResultsFP, "FuturevsPastDistribution/Logbook/", ssp, "_", Mod, "_Begin.csv"))
  
  dir.create(paste0(ResultsFP, "FuturevsPastDistribution/", ssp, "/"), showWarnings = FALSE)
  
  VelsFuture <- list(stack(paste0(ResultsFP, "2020to2080_CMPI6_Velocity_Bioclim/", ssp, "/BioclimVelocity_", Mod, ".grd")))
  names(VelsFuture) <- "Future"
  
  VelsAll <- append(VelsFuture, Vels)
  
  fun <- function(x){if (is.na(x[1])){ NA } else {ecdf(abs(x[2:length(x)]))(abs(x[1]))}}
  
  BioStacks <- stack(lapply(seq(1, 37, 2), function(Bio){
    BioVel <- stack(lapply(VelsAll, function(y) y[[Bio]]))
    DistribCalc <- calc(BioVel, fun)
    return(DistribCalc)
  }))
  
  writeRaster(BioStacks, paste0(ResultsFP, "FuturevsPastDistribution/", ssp, "/Distribution_", Mod, ".grd"))
  write.csv(NULL, paste0(ResultsFP, "FuturevsPastDistribution/Logbook/", ssp, "_", Mod, "_End.csv"))
})








