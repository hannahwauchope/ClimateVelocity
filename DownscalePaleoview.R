##### Script to downscale paleoclimate data (in this case as downloaded from paleoview)

#### Going to use Change Factor method, with kriging. Lovely. (As per this paper: https://www.nature.com/articles/sdata2018254 "PaleoClim, high spatial resolution paleoclimate surfaces for global land areas")
#### It's a very simple method. We've have past rasters at a low resolution (large grid cells), and a present day raster at a high resolution (small grid cells). 
#### We first *upscale* the high res present day raster, to the same res as the past raster, taking the mean of the all the smaller grid cells.
#### We find the difference between the past raster, and this upscaled present day raster, for each grid cell. This is the "change factor". 
#### We then downscale the "change factor" raster, down to the high resolution of the present day raster. To do this, we use "kriging" which linearly interpolates across space
#### Finally, we add the downscaled kriged change factor raster, to the high res present day raster. This gives us a high res past raster.
#### (If someone random reading this is confused I'm afraid you're on your own, but han if you're confused look for "WTF happened w/ downscaling" in your notebook)

#### Present day high res data will be worldclim, 10min res
#### Past data will be data downloaded from paleoview (see this: https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.03031 "PaleoView: a tool for generating continuous climate projections spanning the last 21 000 years at regional and global scales")
#### I've downloaded paleoview data for every 20 years from 0years BP to 21000 years BP (so... a lot of rasters)
#### Going to downscale monthly: precipitation, maxtemp, mintemp and mean temp; plus mean annual temperature

#### Palaeoview gives precipitation data in "an average day in the month", rather than "average monthly total" (which is what is given in worldclim), so we convert this
#NOTE, to streamline this code I've put Mean Annual Temp is in a folder called "January" just because it was easier, but obvs it represents the whole year

#THIS NEEDS TO BE CONVERTED TO TERRA

#### Initialise
cluster <- FALSE 

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

WGSCRS <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
MollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")

#### Set parameters, set up run type (basically a big dataframe of all the rasters I wanna downscale - i.e. monthly for the 4 climate vars, + mean annual temp)
### I've then split each of those combos by three, because this is CHUNKY because of all the many time steps I wanna downscale, so I do them in chunks of three

Months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
MonthDays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
ClimTypes <- c("Precip", "MaxTemp", "MinTemp", "MeanTemp")

Sys.sleep(sample(120, 1))

dir.create(paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/Logbook/"), showWarnings=FALSE)
ClimTypeMonth <- expand.grid(c(1:3), c(1:12), ClimTypes)
ClimTypeMonth$Var3 <- as.character(ClimTypeMonth$Var3)
ClimTypeMonth <- rbind(ClimTypeMonth, data.frame(Var1=c(1:3), Var2=1, Var3="MeanAnnualTemp"))
ClimTypeMonth$Var2 <- as.numeric(as.character(ClimTypeMonth$Var2))

#This is for running on the cluster, it goes through and looks to see which runs have already been started, until it hits one that hasn't begun yet
i <- 1
while(file.exists(paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/Logbook/", ClimTypeMonth[i,3], "_", ClimTypeMonth[i,2], "_", ClimTypeMonth[i,1], "_Begin.csv"))){
 i <- i + 1
}

#It marks that its beginning work on this run
write.csv(NULL, file=paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/Logbook/", ClimTypeMonth[i,3], "_", ClimTypeMonth[i,2], "_", ClimTypeMonth[i,1], "_Begin.csv"))

#And then notes the right chunk, month, and climate type to begin work on
ChunkRun <- ClimTypeMonth[i,1]
ClimMonth <- ClimTypeMonth[i,2]
ClimType <- ClimTypeMonth[i,3]

dir.create(paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/", ClimType, "/"), showWarnings=FALSE)
dir.create(paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/", ClimType, "/", Months[[ClimMonth]], "/"), showWarnings=FALSE)

#### Load Paleoview Data (this has been freshly downloaded from paleoview)
Period <- list.files(paste0(DataFP, "PalaeoView/21000to0BP_20yrInterval_20yrMean/", ClimType, "/", Months[[ClimMonth]]), pattern="*.asc$", full.names=TRUE)
if(cluster==TRUE){
  RasChunk <- if(ChunkRun==1){c(1:350)} else if(ChunkRun==2){c(351:700)} else {c(701:1051)}
  PeriodRas <- pbmclapply(Period[RasChunk], raster, mc.cores=ncores)
} else {
  PeriodRas <- pblapply(Period[1:5], raster)
}

#### Convert palaeoview precip data from "an average day in the month" to "monthly sum"
if(ClimType == "Precip"){
  PeriodRas <- pblapply(PeriodRas, function(x) x * MonthDays[[ClimMonth]])
}

#### Upscale worldclim data to palaeoview res (2.5deg, "agg") and get the right month, and also get the right month at native res (10min)
if(ClimType == "MeanAnnualTemp"){
  WorldClim <- stack(raster(paste0(DataFP, "WorldClim/WC2.1_BioClim_10min_1970-2000/wc2.1_10m_bio_1.tif")))
} else {
  WorldClim <- stack(lapply(list.files(paste0(DataFP, "WorldClim/WC2.1_", ClimType, "_10min_1970-2000"), pattern="*.tif", full.names=TRUE), raster))
}
WorldClimAgg <- resample(WorldClim, PeriodRas[[1]])
WorldClimMonthAgg <- WorldClimAgg[[ClimMonth]]
WorldClimMonth <- WorldClim[[ClimMonth]]
crs(WorldClimMonth) <- WGSCRS
crs(WorldClimMonthAgg) <- WGSCRS

#### Calculate change factor, by subtracting the upscaled worldclim values from the palaeoview values
print("begin CF")
ChangeFactor <- pbmclapply(PeriodRas, function(x){
  crs(x) <- WGSCRS
  CF <- x - WorldClimMonthAgg
  names(CF) <- names(x)
  return(CF)
}, mc.cores=ncores)
print("end CF")

#### Krige the change factor, recombine with high res present day raster, export
## For kriging, you use something called a variogram, all credit to these guys for explaining: https://lazymodellingcrew.com/post/post_05_variowhat_wb/
## Basically, it calculates the interrelationship between two locations that are separated by distance x and direction y

print("BeginKrige")
DownscaledPaleo <- pbmclapply(ChangeFactor, function (i){
  SaveName <- paste0("Downscaled_", ClimType, "_", Months[[ClimMonth]], "_", str_split_fixed(names(i), "_", 3)[,3])
  if(file.exists(paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/", ClimType, "/", Months[[ClimMonth]], "/", SaveName, ".tif"))){
    return(NULL)
  }
  print(paste("Multicore", names(i)))
  CF <- as.data.frame(i, xy=TRUE)
  names(CF) <- c("x", "y", "Climate")
  CF <- CF[complete.cases(CF),]
  CF <- SpatialPointsDataFrame(cbind(CF$x, CF$y), CF, proj4string = WGSCRS)
  
  #NOTE FOR THE FOLLOWING. ACCORDING TO GSTAT, DISTANCES FOR UNPROJECTED DATA (i.e. WGS84) ARE GREAT CIRCLE DISTANCE IN KM
  CFVGM <- variogram(Climate~1, data=CF) #Fit a variogram to change factor in each gridcell (for kriging)
  CFFit <- fit.variogram(CFVGM, model=vgm("Sph")) # fit model
  #plot(CFVGM, CFFit) #have checked plots, they look sensible. 
  
  TenMinUnique <- WorldClimMonth #Make a new raster of same size and resolution as what we want to downscale to (i.e. non-aggregated worldclim)
  values(TenMinUnique) <- 1:length(TenMinUnique) #Change the values so each cell has a unique value
  TenMinPoints <- as.data.frame(rasterToPoints(TenMinUnique)) #Get a point for the centre of each raster gridcell
  TenMinPoints <- SpatialPointsDataFrame(cbind(TenMinPoints$x, TenMinPoints$y), TenMinPoints, proj4string = WGSCRS)
  
  #nmax = 12 as per ECOCLIMATE: A DATABASE OF CLIMATE DATA FROM MULTIPLE MODELS FOR PAST, PRESENT, AND FUTURE FOR MACROECOLOGISTS AND BIOGEOGRAPHERS
  #This value represents the number of cells around which we Krige
  
  CFKrige <- as.data.frame(krige(Climate~1, CF, TenMinPoints, CFFit, nmax=12)) #Obtain kriged values for all 10 min points
  print("finished krige")
  CFKrigeRas <- rasterFromXYZ(CFKrige[,c("coords.x1", "coords.x2", "var1.pred")], crs=WGSCRS) #Convert back into raster
  
  Downscaled <- WorldClimMonth + CFKrigeRas
  writeRaster(Downscaled, paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/", ClimType, "/", Months[[ClimMonth]], "/", SaveName, ".tif"))
  gc()
  return(SaveName)
}, mc.cores=ncores)

write.csv(NULL, file=paste0(DataFP, "PalaeoView/Downscale_21000to0BP_20yrInterval_20yrMean/Logbook/", ClimTypeMonth[i,3], "_", ClimTypeMonth[i,2], "_", ClimTypeMonth[i,1], "_End.csv"))





