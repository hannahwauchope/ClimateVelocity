#Simple script to convert Palaeoview Pre, Max, Min to Bioclim

#NEEDS TO BE CONVERTED TO TERRA 

### 
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
  
  DataFP <- "/nobackup/beegfs/workspace/hw656/Data/"
  ncores=32
  
} else {
  library(gstat)
  library(pbapply)
  library(pbmcapply)
  library(raster)
  library(rgeos)
  library(data.table)
  library(dismo)
  library(stringr)
  
  DataFP <- "/Users/hannahwauchope/Dropbox/Work/Data/"
  ncores=6
}

ClimFP <- paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/")
TimeSteps <- list.files(path=paste0(ClimFP, "MaxTemp/Jan/"), pattern="*.tif")

dir.create(paste0(ClimFP, "BioClim"), showWarnings=FALSE)

write.csv(NULL, file=paste0(ClimFP, "Logbook/Biovars_Begin.csv"))

pbmclapply(TimeSteps, function(x){
  FileNam <- str_split_fixed(str_split_fixed(x, "[.]",2)[,1], "[_]", 4)[,4]
  if(file.exists(paste0(ClimFP, "/BioClim/", FileNam, ".grd"))){
    return("Done")
  }
  MaxMonths <- stack(lapply(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), function(mon) raster(list.files(paste0(ClimFP, "MaxTemp/", mon), paste0("_", FileNam, ".tif$"), full.names = TRUE))))
  MinMonths <- stack(lapply(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), function(mon) raster(list.files(paste0(ClimFP, "MinTemp/", mon), paste0("_", FileNam, ".tif$"), full.names = TRUE))))
  PreMonths <- stack(lapply(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), function(mon) raster(list.files(paste0(ClimFP, "Precip/", mon), paste0("_", FileNam, ".tif$"), full.names = TRUE))))
  BioRasts <- biovars(PreMonths, MinMonths, MaxMonths)
  writeRaster(BioRasts,paste0(ClimFP, "/BioClim/Downscaled_BioVars_", FileNam, ".grd"), format="raster")
  return("Done")
}, mc.cores=ncores)  

write.csv(NULL, file=paste0(ClimFP, "Logbook/Biovars_End.csv"))

