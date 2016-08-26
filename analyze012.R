# When a hexagon is clicked on, we want to display the following things:
# 1) A barplot of samples from across those years (xvalues=1986:2016, yvalues=numSamples)
# 2) A lineplot where the allele frequencies through time are displayed, with potential BTS alleles in color and 
#    suspected native alleles shown in grey
# 3) Some measure of genotyping uncertainty (?)


setwd("~/Box Sync/UCLA/Research/Papers/CTS_HybEx/pondsThroughTime/")

##### LIBRARIES #####
library(tigris)
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(dplyr)
library(tidyr)
library(rgbif)
library(viridis)
library(gridExtra)
library(rasterVis)
library(rnaturalearth)


##### FUNCTIONS #####
make_grid <- function(x, cell_diameter, cell_area, clip = FALSE) {
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  ext <- as(extent(x) + cell_diameter, "SpatialPolygons")
  projection(ext) <- projection(x)
  # generate array of hexagon centers
  g <- spsample(ext, type = "hexagonal", cellsize = cell_diameter, 
                offset = c(0.5, 0.5))
  # convert center points to hexagons
  g <- HexPoints2SpatialPolygons(g, dx = cell_diameter)
  # clip to boundary of study area
  if (clip) {
    g <- gIntersection(g, x, byid = TRUE)
  } else {
    g <- g[x, ]
  }
  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  return(g)
}

# We'll use this from the R cookbook
fillNAgaps <- function(x, firstBack=FALSE) {
  ## NA's in a vector or factor are replaced with last non-NA values
  ## If firstBack is TRUE, it will fill in leading NA's with the first
  ## non-NA value. If FALSE, it will not change leading NA's.
  
  # If it's a factor, store the level labels and convert to integer
  lvls <- NULL
  if (is.factor(x)) {
    lvls <- levels(x)
    x    <- as.integer(x)
  }
  goodIdx <- !is.na(x)
  # These are the non-NA values from x only
  # Add a leading NA or take the first good value, depending on firstBack   
  if (firstBack)   goodVals <- c(x[goodIdx][1], x[goodIdx])
  else             goodVals <- c(NA,            x[goodIdx])
  # Fill the indices of the output vector with the indices pulled from
  # these offsets of goodVals. Add 1 to avoid indexing to zero.
  fillIdx <- cumsum(goodIdx)+1
  x <- goodVals[fillIdx]
  # If it was originally a factor, convert it back
  if (!is.null(lvls)) {
    x <- factor(x, levels=seq_along(lvls), labels=lvls)
  }
  x
}


##### Create plot of hexes colored by BTS ancestry #####
# This will be what users will click on when they're interacting with the map

## We'll use this to calculate the percentage BTS within the hexagon and color the hex accordingly
agOverHexNoFillMissing <- function(hexGrid, pointsLayer) {
  rbPal <- colorRamp(c('red','blue'), space = "Lab") # needed to color points based on 0 to 1 values
  
  hexSubset <- hexGrid[unique(sp::over(pointsLayer, hexGrid)),]
  
  agBTS <- over(hexSubset, pointsLayer, returnList = TRUE)
  agBTShex <- plyr::ldply(agBTS, .fun=function(x) x, .id="id") # Assign the hexagon ID to each observation
  agBTShexChar <- mutate(agBTShex,id = as.character(id))
  # Now take all the agBTShex hexagons and aggregate based on id values
  aveAgBTShexes <- aggregate(PercBTS ~ id, data=agBTShexChar, FUN=mean) # Just take the mean of all the values in the hexagon
  aveAgBTShexes <- mutate(aveAgBTShexes, id=as.character(id))
  aveAgBTShexes$col <- rbPal(aveAgBTShexes$PercBTS)
  row.names(aveAgBTShexes) <- aveAgBTShexes$id
  
  aveAgBTShexesSP <- SpatialPolygonsDataFrame(hexSubset, aveAgBTShexes)
  
  #plot(aveAgBTShexesFullSP, col=aveAgBTShexesFullSP$col, add=T, lwd=0.001)
  return(aveAgBTShexesSP)
}

##### Spatial mander data #####
calCounties <- counties(state="CA", cb=TRUE) # Using package "tigris" to download TIGR data
calCountiesUTM <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=km +no_defs") %>% spTransform(calCounties, .)
countiesWithMandersUTM <- calCountiesUTM[calCountiesUTM$COUNTYFP=="001"|calCountiesUTM$COUNTYFP=="009"|calCountiesUTM$COUNTYFP=="013"|calCountiesUTM$COUNTYFP=="019"|calCountiesUTM$COUNTYFP=="033"|calCountiesUTM$COUNTYFP=="039"|calCountiesUTM$COUNTYFP=="041"|calCountiesUTM$COUNTYFP=="043"|calCountiesUTM$COUNTYFP=="047"|calCountiesUTM$COUNTYFP=="053"|calCountiesUTM$COUNTYFP=="067"|calCountiesUTM$COUNTYFP=="069"|calCountiesUTM$COUNTYFP=="077"|calCountiesUTM$COUNTYFP=="079"|calCountiesUTM$COUNTYFP=="083"|calCountiesUTM$COUNTYFP=="085"|calCountiesUTM$COUNTYFP=="087"|calCountiesUTM$COUNTYFP=="095"|calCountiesUTM$COUNTYFP=="097"|calCountiesUTM$COUNTYFP=="099"|calCountiesUTM$COUNTYFP=="101"|calCountiesUTM$COUNTYFP=="107"|calCountiesUTM$COUNTYFP=="113"|calCountiesUTM$COUNTYFP=="029"|calCountiesUTM$COUNTYFP=="081",] # Added Kern and San Mateo County because there are salamanders very close by
# Make a single polygon that consists of all the salamander counties so that county lines don't break up hexagons later
manderCountyUTMPolyLabelPts <- getSpPPolygonsLabptSlots(countiesWithMandersUTM) 
manderCountyUTMIDOneBin <- cut(manderCountyUTMPolyLabelPts[,1], range(manderCountyUTMPolyLabelPts[,1]), include.lowest=TRUE)
manderCountiesUTMDissolved <- unionSpatialPolygons(countiesWithMandersUTM ,manderCountyUTMIDOneBin)


# Read in the genotype data, merge it with the sample metadata, and turn it into a spatial dataframe:
genotypes <- read.csv("HSEM031-040_and_HSET01.SNPs.sr.q30.passOnly.minGQ20maxMiss50p.erinRename.hasMeta.oneHeaderLine.noSampsOver50pmissing.012", header=F, sep="\t", row.names=1)

# Shorten it up just for testing:
genotypes <- genotypes[,1:1000] # Shorten it up just for testing:

individuals <- read.csv("HSEM031-040_and_HSET01.SNPs.sr.q30.passOnly.minGQ20maxMiss50p.erinRename.hasMeta.oneHeaderLine.noSampsOver50pmissing.012.indv", header=F)
loci <- read.csv("HSEM031-040_and_HSET01.SNPs.sr.q30.passOnly.minGQ20maxMiss50p.erinRename.hasMeta.oneHeaderLine.noSampsOver50pmissing.012.pos", header=F)

# Shorten up just for testing:
loci <- loci$V1[1:1000]

colnames(genotypes) <- loci$V1
row.names(genotypes) <- individuals$V1
metaz <- read.csv("allMetaData", sep="\t", header=T)
genotypeWithLatLong <- merge(genotypes, metaz, by.x=0, by.y="Individual")
genotypeWithLatLong$Long <- abs(genotypeWithLatLong$Long) * -1 # Make sure that the longitudes are negative
# Project the samples to UTM:
coordinates(genotypeWithLatLong) <- c("Long", "Lat")
proj4string(genotypeWithLatLong) <- CRS("+proj=longlat +datum=WGS84")
genotypeWithLatLongUTM <- spTransform(genotypeWithLatLong, CRS("+proj=utm +zone=10 +datum=WGS84 +units=km +no_defs"))


fiveStarIDs <- c(109680,109695,38114,38115,38116,38117,38118,38119,38120,38121,38122,38123,38124,38125,38126,38127,38128,38129,38130,38131,38132,38133,38134,38135,38137,38139,38140,26869,109693,109696)
fiveStarRows <- match(fiveStarIDs,rownames(genotypes))
fiveStarGenos <- genotypes[fiveStarRows,]

pureCTSids <- c(109699, 109700, 109701, 109702, 109703, 109704, 109705, 109706, 109707, 109708, 109709, 109710, 109711, 109712, 109713, 19276, 19277, 19278, 19279, 19280, 19281, 19282, 19283, 19284, 19285, 19286 , 106188, 109714, 118294, 124117, 125597, 125600, 125647, 6694, 6700, 8307, 8331, 8335, 8357, 8855, 106497, 119977, 122381, 6693, 106257, 119974, 105862, 106201, 106202, 106276, 118274, 118285, 106496, 106498, 119961, 124114, 8856, 32664, 32720, 32735, 32762, 32815, 33077, 38793, 38807, 38816, "6100M", 32636, 32821, 33085, 32763, 32766, 32778, 32826, 32838, 33084, 33089, 32628, 32635, 32665, 32689, 32717, 32738, 126241, 16761, 108766, 112960, 14379, 14388, 126244, 110585, 110626, 110658, 110683, 110588, 110597, 110619, 110628, 110647, 110651, 110668, 110669 , 106335, 107081, 107084, 107136, 107155, 110375, 110390, 110416, 110431, 110441, 110452, 110459, 110478, 110513, 110543, 11632, 124005, 124012, 124013, 124017, 14347, 14365, 21713, 26364, 26368, 26402, 110335, 110351, 110358, 110397, 110554, 107087, 110518, 110548, 110563, 124001, 124018, 124019, 124020, 124021, 124022, 16722, 10529, 110409, 110430, 110482, 110484, 110503, 110512, 110895, 110330, 110338, 110343, 110350, 110355, 110362, 110376, 110389, 110395, 110402, 110403, 110553, 110557, 110568, 110569, 14353, 16717, 26408, 27744, 110421, 110426, 110433, 110436, 110439, 110454, 110455, 110499, 110504, 110505, 110507, 110522, 110540, 121488, 124002, 21700, 21705, 9800, 9808, 9809, 9817, 9801)
pureCTSrows <- match(pureCTSids, rownames(genotypes))
pureCTSgenos <- genotypes[pureCTSrows,]


# Now lets find all the columns (SNPs) where all the sum of all the genotyped alleles in CTS is 0 and is non-zero
# in the Five Star Fish Farm individuals
pureCTSfix0names <- names(pureCTSgenos[,colSums(pureCTSgenos, na.rm=T)==0])
pureBTSfix1names <- names(fiveStarGenos[,(apply(fiveStarGenos,2,function(x) sum(is.na(match(c(0,1),x))))==2)])

fix0CTSfix1BTS <- intersect(pureCTSfix0names,pureBTSfix1names)
diagnosticGenos <- genotypes[,fix0CTSfix1BTS]
# Now calculate the percentage BTS ancestry in each individual
diagnosticGenos$PercBTS <- (apply(diagnosticGenos, 1, function(x) (sum(x,na.rm=T)/(sum(!is.na(x))*2))))
# Now generate the map
diagnosticGenosWithLatLong <- merge(diagnosticGenos, metaz, by.x=0, by.y="Individual")
diagnosticGenosWithLatLong$Long <- abs(diagnosticGenosWithLatLong$Long) * -1
coordinates(diagnosticGenosWithLatLong) <- c("Long", "Lat")
proj4string(diagnosticGenosWithLatLong) <- CRS("+proj=longlat +datum=WGS84")
diagnosticGenosWithLatLongUTM <- spTransform(diagnosticGenosWithLatLong, CRS("+proj=utm +zone=10 +datum=WGS84 +units=km +no_defs"))


# Make the hexagon grid
manderCountiesUTMDissolvedHex <- make_grid(manderCountiesUTMDissolved, cell_area = 200, clip=TRUE)
# Calculate the average BTSness in the grid cells and color hexagons accordingly
overlappingPoints <- diagnosticGenosWithLatLongUTM[!is.na(sp::over(diagnosticGenosWithLatLongUTM,as(manderCountiesUTMDissolved,"SpatialPolygons"))),]
countiesWithMandersHexcolors <- agOverHexNoFillMissing(manderCountiesUTMDissolvedHex,overlappingPoints)
plot(manderCountiesUTMDissolved, axes=T, main="Hexagons colored by percentage fixed-difference\nancestry for illustrative purposes.")
plot(calCountiesUTM, lwd=1, add=T)
#plot(lakesUTM, col="midnightblue", add=T)
#plot(oceansClippedUTM, col="midnightblue", add=T)
#plot(roadsI_UTM, col="grey", lwd=.7, add=T, lty=2) # Interstate
#plot(roadsU_UTM, col="grey", lwd=0.7, add=T, lty=2) # U.S.
# Plot the hexagons
plot(countiesWithMandersHexcolors, col=rgb(countiesWithMandersHexcolors$col, maxColorValue = 255), add=T, lwd=0.001)
# Legend and scale bar

#save(countiesWithMandersHexcolors, file="countiesWithMandersHexcolors.RData")
#save(manderCountiesUTMDissolved, file="manderCountiesUTMDissolved.RData")
#save(calCountiesUTM, file="calCountiesUTM.RData")






###### Now calculate the allele frequences in each hex through time ######


overlappingPoints <- genotypeWithLatLongUTM[!is.na(sp::over(genotypeWithLatLongUTM,as(manderCountiesUTMDissolved,"SpatialPolygons"))),]

# Use the same hex overlay grid and points: manderCountiesUTMDissolvedHex and overlappingPoints
hexSubset <- manderCountiesUTMDissolvedHex[unique(sp::over(overlappingPoints, manderCountiesUTMDissolvedHex)),]

agBTS <- over(hexSubset, overlappingPoints, returnList = TRUE)
agBTShex <- plyr::ldply(agBTS, .fun=function(x) x, .id="id") # Assign the hexagon ID to each observation
agBTShexChar <- mutate(agBTShex,id = as.character(id))
# Now create a vectors for every allele in every hexagon that stores the allele frequency in that hexagon through time
# agBTShexChar$id is the hexagon identifier

# The resulting dataframe needs to store a vector for each hexagon for each allele that is the length of 
# the sampling period (eg 1986:2015 or 1986:2016)

# Set it up like array[hexagonID][locusID][year]
freqThroughTimeArray <- array(dim=c(length(unique(agBTShexChar$id)), length(3:(ncol(agBTShexChar)-1)), length(min(agBTShexChar$Year):max(agBTShexChar$Year))), dimnames=list(unique(agBTShexChar$id), colnames(agBTShexChar)[3:(ncol(agBTShexChar)-1)], min(agBTShexChar$Year):max(agBTShexChar$Year)))

# Set up like array[hexagonID][year]
hexSampsPerYear <- array(dim=c(length(unique(agBTShexChar$id)), length(min(agBTShexChar$Year):max(agBTShexChar$Year))), dimnames=list(unique(agBTShexChar$id), min(agBTShexChar$Year):max(agBTShexChar$Year)))   
  
  
hexIDs <- unique(agBTShexChar$id)
for (hexID in 1:length(hexIDs)) {
   hexFrame <- agBTShexChar[agBTShexChar$id==hexIDs[hexID],]
   # hexFrame now holds the observations for each sample within hexagon hexIDs[hexID], including $Year
   
   # Now we have to calculate the allele frequency for each year and store it in a vector:
   for (year in 1:length(1986:2015)) {
     hexFrameYear <- hexFrame[hexFrame$Year==(min(agBTShexChar$Year)-1+year),] # The earliest year minus 1 plus the "year" counter
     hexSampsPerYear[hexID,year] <- base::nrow(hexFrameYear)
     
   for (allele in 1:(ncol(hexFrameYear)-4)) {
         alleleCounter = allele+2 # The first two columns are not data
       # The first column is the hexagon id, the second is sample ID, the third is the first
       # genotype column, the second to last is the last genotype column, and the last is the year
       if (length(hexFrameYear[,alleleCounter]) == 0) {
         freqThroughTimeArray[hexID,allele,year] <- NA
       } else {
         numGenotypesCalled <- sum(2 * !is.na(hexFrameYear[,alleleCounter]))
         numAlternateAlleles <- sum(hexFrameYear[,alleleCounter])
         altAlleleFreq <- numAlternateAlleles / numGenotypesCalled
         freqThroughTimeArray[hexID,allele,year] <- altAlleleFreq
       }
       
     }
   }
}

save(freqThroughTimeArray, file="freqThroughTimeArray.RData")
save(genotypes, file="genotypes.RData")
save(hexSampsPerYear, file="hexSampsPerYear.RData")

# We have to deal with the missing values in a sensible way. We want the allele to "appear" on the graph when 
# the hexagon was first sampled, but we do not want any breaks after that until its final sampling event

#for (allele in 1:dim(freqThroughTimeArray)[2]) {
#  for (hex in 1:dim(freqThroughTimeArray)[1]) {
#    # The vector we're working with here is freqThroughTimeArray[hex,allele,], which holds the allele frequency of thenon-ref allele each year from 1986 to 2015
#    firstPresent <- which.min(is.na(freqThroughTimeArray[hex,allele,]))
#  }
#}




#for (hex in 1:1) {
#   Sys.sleep(1)
#   plot(1,ylim=c(0,1), xlim=c(1986,2015))
#   for (allele in 1:1000) {
#     lines(1986:2015, fillNAgaps(freqThroughTimeArray[hex,allele,]), main=allele, type="l", col="black")
#   }
#}




