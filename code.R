library(sp)
library(raster)
setOptions(todisk=FALSE)
setOptions(chunksize = 1e+06, maxmemory = 1e+07)
library(maptools)
library(gstat)

library(lattice)
library(latticeExtra)
library(rasterVis)

library(solaR)

library(parallel)

library(xtable)
library(colorspace)
library(classInt)

## Change the folder name if needed
## This folder must include three subfolders
## with names SIAR, data, and figs
setwd('~/Investigacion/DocsPropios/CMSAF_SIAR/')

################################################################################
## HORIZONTAL PLANE
################################################################################

##########
## SIAR
##########

SIAR <- read.csv('http://solar.r-forge.r-project.org/data/SIAR.csv')

## SIAR webpage has changed and the function readSIAR does not
## work with it. The result of this code is contained in the
## spainMeteo20102011.RData file

## spainMeteo <- apply(SIAR[, c(7, 3, 1)], 1,
##                     function(x){
##                       try(readSIAR(prov=x[2], est=x[3],
##                                start='01/01/2010', end='31/12/2011',
##                                lat=x[1]))
##                     }
##                     )
## idxMeteo <- sapply(spainMeteo, function(x)class(x)=='Meteo')

## spainMeteoOK<-spainMeteo[idxMeteo]
## save(spainMeteo, idxMeteo, spainMeteoOK, 'spainMeteo20102011.RData')

##########
## CMSAF
##########
projSIAR <- CRS('+proj=longlat +ellps=WGS84')

## This code must be run with daily data previously downloaded from
## the CMSAF web service.  Change it to a working folder where the
## CMSAF data is available.
old <- setwd('~/Datos/CMSAF/2010_2011')

listFich <- dir(pattern='SISdm2010')
stackSIS <- stack(listFich)
stackSIS <- stackSIS*24 ##para pasar de W/m2 (irradiancia media) a Wh/m2
projection(stackSIS) <- projSIAR

listFich <- dir(pattern='SISdm2011')
stackSIS <- stack(listFich)
stackSIS <- stackSIS*24 
projection(stackSIS) <- projSIAR

setwd(old)

writeRaster(stackSIS, filename='data/SISd2011', overwrite=TRUE)

G0yCMSAF2010 <- calc(SISd2010, fun=function(x)sum(x, na.rm=1)/1000,
                     filename='data/G0yCMSAF2010')
G0yCMSAF2011 <- calc(SISd2011, fun=function(x)sum(x, na.rm=1)/1000,
                     filename='data/G0yCMSAF2011')
G0yCMSAF <- (G0yCMSAF2010 + G0yCMSAF2011)/2

writeRaster(G0yCMSAF, filename='data/G0yCMSAF')

################################################################################
## EFFECTIVE IRRADIATION
################################################################################


#########
## SIAR
#########
load('data/spainMeteo20102011.RData')

#### Cleaning...

## Number of days
nDays <- function(x)as.numeric(diff(range(indexD(x))))
ndays <- sapply(spainMeteoOK, nDays)

## Yearly sums (kWh)
meanYearlySums <- function(x)mean(aggregate(getG0(x), year, sum, na.rm=1))

spainG0y <- sapply(spainMeteoOK, meanYearlySums)/1000

clean <- (ndays>=600) & (spainG0y>=1000)

spainMeteoOK <- spainMeteoOK[clean]


gefList <- function(x){
  lat <- getLat(x, 'deg')
  cat(x@source, '\n')
  gefFixed <- calcGef(lat=lat, dataRad=x, modeRad='bd', modeTrk='fixed')
  gef2x <- calcGef(lat=lat, dataRad=gefFixed, modeRad='prev', modeTrk='two')
  gefHoriz <- calcGef(lat=lat, dataRad=gefFixed, modeRad='prev', modeTrk='horiz')
  resultFixed <- as.zooD(gefFixed)$Gefd
  result2x <- as.zooD(gef2x)$Gefd
  resultHoriz <- as.zooD(gefHoriz)$Gefd

  result <- cbind(resultFixed, result2x, resultHoriz)

  ## Those days with G0d>Bo0d or with G0d too low are erroneous
  Bo0d <- gefFixed@solD$Bo0d
  G0d <- gefFixed@G0D$G0d
  err <- (G0d>Bo0d) | (G0d < 0.01*Bo0d) 
  result[err] <- NA
  
  aggregate(result, by=year, sum, na.rm=1)/1000
  }

ncores <- detectCores()
gefSIAR <- mclapply(spainMeteoOK, FUN = gefList, mc.cores=ncores)

save(gefSIAR, file='data/gefSIAR.RData')

#########
## CMSAF
#########

gefParallel <- function(data, filename="", nodes=detectCores(), blocks=8,...){
  ## latitude values as a new raster
  y <- init(data, v='y')
  idx <- getZ(data)

  bs <- blockSize(data, minblocks=blocks*nodes)

  gef <- function(g0, idx){
    n <- length(g0)
    lat <- g0[1]
    g0d <- list(file = zoo(data.frame(G0 = g0[2:n]),  idx),
                lat = lat)
    gefFixed <- calcGef(lat=lat, dataRad =  g0d, 
                        modeRad = 'bd',  modeTrk='fixed')
    gef2x <- calcGef(lat = lat, dataRad = gefFixed, 
                     modeRad = 'prev',  modeTrk='two')
    gefHoriz <- calcGef(lat =  lat, dataRad = gefFixed, 
                        modeRad = 'prev',  modeTrk='horiz')
    resultFixed <- as.numeric(as.data.frameY(gefFixed)$Gefd)
    result2x <- as.numeric(as.data.frameY(gef2x)$Gefd)
    resultHoriz <- as.numeric(as.data.frameY(gefHoriz)$Gefd)
    result <- c(resultFixed, result2x, resultHoriz)
    result
  }

  ## List with the indices of blocks for each node
  iCluster <- splitIndices(bs$n, nodes) 
  resCl <- mclapply(iCluster,
                    ## Each node receives an element of iCluster, a set of indices
                    function(icl){
                      resList <- lapply(icl, function(i){
                        ## An element of icl is the number of block to
                        ## be read by each node
                        vals <- getValues(data, bs$row[i], bs$nrows[i])
                        lat <- getValues(y, bs$row[i], bs$nrows[i])
                        vals <- cbind(lat, vals)
                        cat(i, ':', range(lat), '\n')
                        res0 <- apply(vals, MARGIN=1L, FUN=gef, idx)
                        cat(i, ':', range(res0), '\n')
                        res0
                      })
                      ## The result of lapply is a list of matrices
                      ## with 3 rows (one per each tracking system).
                      ## I use cbind() to join columns and t() to transpose to
                      ## get a matrix with 3 columns. Each row is a
                      ## cell of the original raster.
                      t(do.call(cbind, resList))
                    }, mc.cores=nodes)
  ## The result of mclapply is a list with as many elements as nodes
  ## Each element of the list is a matrix with 3 columns (resCl0)
  ## corresponding to a block as defined by bs.
  resCl <- do.call(rbind, resCl)
  
  out <- brick(data, nl=3) ## 3 layers, one for each tracking mode
  layerNames(out)=c('Fixed', 'Two', 'Horiz')
  out <- setValues(out, resCl)
  if (filename!='') out <- writeRaster(out, filename=filename)
  out
}

ncores <- detectCores()

SISd2010 <- brick('data/SISd2010.grd')
SISd2011 <- brick('data/SISd2011.grd')

idx2010 <- fBTd('serie', year=2010)
SISd2010 <- setZ(SISd2010, idx2010)
layerNames(SISd2010) <- as.character(idx2010)

idx2011 <- fBTd('serie', year=2011)
SISd2011 <- setZ(SISd2011, idx2011)
layerNames(SISd2011) <- as.character(idx2011)

gef2010 <- gefParallel(SISd2010, filename='data/gef2010')

gef2011 <- gefParallel(SISd2011, filename='data/gef2011')

gefCMSAF <- (gef2010 + gef2011)/2
layerNames(gefCMSAF)=c('Fixed', 'Two', 'Horiz')

writeRaster(gefCMSAF, filename='data/gefCMSAF')


################################################################################
## KRIGING WITH EXTERNAL DRIFT
################################################################################


##################################################
######  HORIZONTAL IRRADIATION
##################################################

##########
## SIAR
##########

load('data/spainMeteo20102011.RData')

SIAR <- read.csv('http://solar.r-forge.r-project.org/data/SIAR.csv')
SIAR <- SIAR[idxMeteo,]

nDays <- function(x)as.numeric(diff(range(indexD(x))))
ndays <- sapply(spainMeteoOK, nDays)

meanYearlySums <- function(x)mean(aggregate(getG0(x), year, sum, na.rm=1))

spainG0y <- sapply(spainMeteoOK, meanYearlySums)/1000

clean <- (ndays>=600) & (spainG0y>=1000)

spainMeteoOK <- spainMeteoOK[clean]
SIAR <- SIAR[clean,]

g0SIAR <- lapply(spainMeteoOK, getG0)
g0SIAR <- do.call(cbind, g0SIAR)
names(g0SIAR) <- SIAR$Estacion

BTd <- index(g0SIAR)
Bo0d <- lapply(spainMeteoOK, function(rad){
  lat <- getLat(rad, units='deg')
  fSolD(lat=lat, BTd=BTd)$Bo0d
})

Bo0d <- do.call(cbind, Bo0d)
names(Bo0d) <-names(g0SIAR)

coredata(g0SIAR)[g0SIAR<0.01*Bo0d] <- NA
coredata(g0SIAR)[g0SIAR>Bo0d] <- NA

spainG0y <- colMeans(aggregate(g0SIAR, by=year, sum, na.rm=1))/1000

SIAR$G0ySIAR=spainG0y

projSIAR <- CRS('+proj=longlat +ellps=WGS84')
spSIAR <- SpatialPointsDataFrame(SIAR[, c(6, 7)], SIAR[, -c(6, 7)],
                                 proj4str=projSIAR)


##########
## CMSAF
##########

SISd2010 <- brick('data/SISd2010.grd')

idx2010 <- fBTd('serie', year=2010)
SISd2010 <- setZ(SISd2010, idx2010)
names(SISd2010) <- as.character(idx2010)

SISd2011 <- brick('data/SISd2011.grd')

idx2011 <- fBTd('serie', year=2011)
SISd2011 <- setZ(SISd2011, idx2011)
names(SISd2011) <- as.character(idx2011)

G0yCMSAF <- raster('data/G0yCMSAF.grd')
G0yCMSAF2010 <- raster('data/G0yCMSAF2010.grd')
G0yCMSAF2011 <- raster('data/G0yCMSAF2011.grd')

spSIAR$G0yCMSAF <- extract(G0yCMSAF, spSIAR)

spSIAR$difG0y <-spSIAR$G0yCMSAF-spSIAR$G0ySIAR
spSIAR$outTolerance <- abs(spSIAR$difG0y/spSIAR$G0ySIAR)>0.05

##############################
## Comparison SIAR-CMSAF (daily values)
##############################

g0CMSAF2010 <- extract(SISd2010, spSIAR)
g0CMSAF2011 <- extract(SISd2011, spSIAR)
g0CMSAF <- rbind(t(g0CMSAF2010), t(g0CMSAF2011))

g0Dif <- g0CMSAF-g0SIAR

applyStats <- function(p, o, ...){
  d <- p - o
  meanO <- colMeans(o, na.rm=TRUE)
  SDp <- apply(p, 2, sd, na.rm=TRUE)
  SDo <- apply(o, 2, sd, na.rm=TRUE)
  MBD <- colMeans(d, na.rm=TRUE)
  RMSDc <- apply(d, 2, sd, na.rm=TRUE)
  RMSD <- sqrt(MBD^2 + RMSDc^2)
  MAD <- apply(d, 2, FUN=function(x)mean(abs(x), na.rm=TRUE))
  res <- data.frame(sdP=SDp, sdO=SDo, MBD=MBD/meanO, MBDabs=abs(MBD)/meanO,
                    RMSDc=RMSDc/meanO, RMSD=RMSD/meanO, MAD=MAD/meanO)
  res
}
  
g0Stats <- applyStats(g0CMSAF, g0SIAR)
names(g0Stats) <- paste(names(g0Stats), 'G0', sep='.')
spSIAR <- spCbind(spSIAR, g0Stats)



##################################################
#### KED
##################################################

krigeRaster <- function(formula, data, raster, ...){
  latLayer <- init(raster, v='y')
  lonLayer <- init(raster, v='x')

  grd <- as(stack(lonLayer, latLayer, raster), 'SpatialGridDataFrame')
  names(grd) <- c('lon', 'lat', deparse(substitute(raster)))

  proj4string(data) <- proj4string(grd) ## seems to be a bug in proj4string<- with SpatialGrid
  resSP <- krige(formula, data, grd, ...)
  res <- as(resSP, 'RasterStack')
  names(res) <- c('pred', 'var')
  projection(res) <- projection(raster)
  res
}

## Universal kriging
vgmCMSAF <- variogram(G0ySIAR~G0yCMSAF, spSIAR)##, width=10, cutoff=500)
plot(vgmCMSAF)

fitvgmCMSAF <- fit.variogram(vgmCMSAF, vgm(model='Nug')) ##pure nugget effect...no spatial autocorrelation
plot(vgmCMSAF, fitvgmCMSAF)
G0yKrig <- krigeRaster(G0ySIAR~G0yCMSAF, spSIAR, G0yCMSAF, model=fitvgmCMSAF)

##################################################
#### Effective Irradiation
##################################################

##########
## SIAR
##########

load('data/gefSIAR.RData')

gefSIARMean <- do.call(rbind, lapply(gefSIAR, colMeans, na.rm=1))
colnames(gefSIARMean) <- c('FixedSIAR', 'TwoSIAR', 'HorizSIAR')

SIARGef <- cbind(SIAR, gefSIARMean)
spGef <- spCbind(spSIAR, as.data.frame(gefSIARMean))

##########
## CMSAF
##########

gefCMSAF2010 <- brick('data/gef2010')
gefCMSAF2011 <- brick('data/gef2011')
gefCMSAF <- brick('data/gefCMSAF')

gefExtract <- as.data.frame(extract(gefCMSAF, spGef))
names(gefExtract) <- paste(layerNames(gefCMSAF), 'CMSAF', sep='')

spGef <- spCbind(spGef, gefExtract)
spGef$difFixed <- spGef$FixedCMSAF - spGef$FixedSIAR
spGef$difTwo <- spGef$TwoCMSAF - spGef$TwoSIAR
spGef$difHoriz <- spGef$HorizCMSAF - spGef$HorizSIAR


##################################################
#### KED
##################################################

FixedCMSAF <- raster(gefCMSAF, layer='Fixed')

vgmFixed <- variogram(FixedSIAR~G0yCMSAF, data=spGef)##, width=10, cutoff=500)
plot(vgmFixed)

modelFixed <- vgm(model='Nug')
fitvgmFixed <- fit.variogram(vgmFixed, modelFixed)
plot(vgmFixed, fitvgmFixed)

FixedKrig <- krigeRaster(FixedSIAR~G0yCMSAF, spGef, G0yCMSAF, model=fitvgmFixed)


HorizCMSAF <- raster(gefCMSAF, layer='Horiz')

vgmHoriz <- variogram(HorizSIAR~G0yCMSAF, data=spGef)
plot(vgmHoriz)
modelHoriz <- vgm(psill=17000, model='Sph', range=200, nugget=5000)
fitvgmHoriz <- fit.variogram(vgmHoriz, modelHoriz)
plot(vgmHoriz, fitvgmHoriz)

HorizKrig <- krigeRaster(HorizSIAR~G0yCMSAF, spGef, G0yCMSAF, model=fitvgmHoriz)


TwoCMSAF <- raster(gefCMSAF, layer='Two')

vgmTwo <- variogram(TwoSIAR~G0yCMSAF, data=spGef)
plot(vgmTwo)
modelTwo <- vgm(psill=17000, model='Sph', range=200, nugget=5000)
fitvgmTwo <- fit.variogram(vgmTwo, modelTwo)
plot(vgmTwo, fitvgmTwo)

TwoKrig <- krigeRaster(TwoSIAR~G0yCMSAF, spGef, G0yCMSAF, model=fitvgmTwo)


## Extract points from kriging
krigExtract <- data.frame(G0yKrig=extract(raster(G0yKrig, 1), spGef),
                          FixedKrig=extract(raster(FixedKrig, 1), spGef),
                          HorizKrig=extract(raster(HorizKrig, 1), spGef),
                          TwoKrig=extract(raster(TwoKrig, 1), spGef))

spGef <- spCbind(spGef, krigExtract)
spGef$difG0yKrig <- spGef$G0yKrig - spGef$G0ySIAR
spGef$difFixedKrig <- spGef$FixedKrig - spGef$FixedSIAR
spGef$difTwoKrig <- spGef$TwoKrig - spGef$TwoSIAR
spGef$difHorizKrig <- spGef$HorizKrig - spGef$HorizSIAR

save(spSIAR, g0SIAR, g0CMSAF, g0Dif,
     G0yKrig,
     spGef, FixedKrig, HorizKrig, TwoKrig,
     file='data/krig.RData')


## Display variograms together
vgmCMSAF$id <- 'G0'
vgmFixed$id <- 'Fixed'
vgmHoriz$id <- 'Horiz.'
vgmTwo$id <- 'Two'

v <- list(G0=vgmCMSAF, Fixed=vgmFixed, Horiz=vgmHoriz, Two=vgmTwo)
m <- list(G0=fitvgmCMSAF, Fixed=fitvgmFixed, Horiz=fitvgmHoriz, Two=fitvgmTwo)

modelLine <- function(v, m){
  l <- variogramLine(m, max(v$dist), n=nrow(v))
  names(l) <- c('dist.model', 'gamma.model')
  l
  }

vList <- mapply(modelLine, v, m, SIMPLIFY=FALSE)

vm <- cbind(do.call(rbind, v), do.call(rbind, vList))

pV <- xyplot(gamma ~ dist, data=vm, groups=id,
             xlab='Distance (km)', ylab='Semivariance') +
  glayer(panel.text(x[1], y[1], group.value, pos=3))
pM <- xyplot(gamma.model ~ dist.model, data=vm, groups=id, type='l')

trellis.device(pdf, file='figs/variograms.pdf')
pV + pM
dev.off()

################################################################################
## GRAPHICS
################################################################################
## setwd where the github repository is located

##############################
## Bubbles
##############################
panel.bubbles <- function(x, y, subscripts, col, cex, ...){
  panel.points(x, y, col=col[subscripts], cex=cex[subscripts], ...)
  }

bubbles <- function(obj, n=7, style='fisher',
                    pal=brewer.pal(name='Blues', n=9),
                    size=c(0.3, 1.1), pwr.size=1, alpha=0.7,...){
  
  xlim = sp:::bbexpand(bbox(obj)[1, ], 0.04)
  ylim = sp:::bbexpand(bbox(obj)[2, ], 0.04)
  scales <- list(draw = TRUE, cex=0.7)

  xy <- as.data.frame(coordinates(obj))
  z <- stack(obj@data)
  z$ind <- factor(z$ind, levels=names(obj))
  data <- cbind(xy, z)
  
  intervals <- classIntervals(z$values, n=n, style=style)
  nInt <- length(intervals$brks) - 1

  idx <- findCols(intervals)
  
  op <- options(digits=2)
  tab <- classInt:::tableClassIntervals(cols = idx, brks = intervals$brks,
                                        under = "under", over = "over", between = "-", 
                                        cutlabels = TRUE,
                                        intervalClosure = "left",
                                        dataPrecision = NULL)
  options(op)

  rval <- seq(1, 0, length=nInt)
  cex.key <- size[2] - diff(size)*rval^pwr.size 
  cex <- cex.key[idx]
  pal <- colorRampPalette(pal)(nInt)
  col <- findColours(intervals, pal)
  col.key=attr(col, 'palette')

  key <- list(space='right',
              text=list(labels=names(tab), cex=0.85),
              points=list(col=col.key, pch=19, cex=cex.key))

  p <- xyplot(lat~lon|ind, data=data,
         xlab='', ylab='', 
         cex=cex, col=col, key=key, alpha=alpha,
         asp=mapasp(obj),
         scales=sp:::longlat.scales(obj, scales, xlim, ylim),
         panel=panel.bubbles, ...)
  p
}


#########################
## CMSAF daily data
#########################
SISd2010 <- brick('data/SISd2010.grd')

idx2010 <- fBTd('serie', year=2010)
SISd2010 <- setZ(SISd2010, idx2010)
names(SISd2010) <- as.character(idx2010)

SISd2011 <- brick('data/SISd2011.grd')

idx2011 <- fBTd('serie', year=2011)
SISd2011 <- setZ(SISd2011, idx2011)
names(SISd2011) <- as.character(idx2011)


##############################
## hovmoller
##############################
SISd <- stack(c('data/SISd2010.grd',
                'data/SISd2011.grd'))
names(SISd) <- c(names(SISd2010), names(SISd2011))
SISd <- setZ(SISd, c(idx2010, idx2011))

pHov <- hovmoller(SISd, par.settings=BTCTheme, add.contour=FALSE)
trellis.device(pdf, file='figs/hovmoller.pdf')
pHov
dev.off()

##################################################################
## Auxiliary data
##################################################################

## Spanish administrative boundaries
old <- setwd(tempdir())
download.file('http://www.gadm.org/data/shp/ESP_adm.zip', 'ESP_adm.zip')
unzip('ESP_adm.zip')
projSIAR <- CRS(projection(spSIAR))
mapaSHP <- readShapeLines('ESP_adm2.shp', proj4string=projSIAR)
setwd(old)

## Spanish altitude mask
old <- setwd(tempdir())
download.file('http://www.diva-gis.org/data/msk_alt/ESP_msk_alt.zip', 'ESP_msk_alt.zip')
unzip('ESP_msk_alt.zip')
elevMask <- raster('ESP_msk_alt')
setwd(old)


## France, Andorra, Portugal and Morocco boundaries (from Natural Earth Data)
neighbours <- readShapePoly('data/neighbours', proj4string=CRS('+proj=longlat +ellps=WGS84'))


##################################################################
## Kriging results
##################################################################
load('data/krig.RData')

datGef <- as.data.frame(spGef)

G0yCMSAF <- raster('data/G0yCMSAF.grd')
gefCMSAF <- brick('data/gefCMSAF')

elevMask <- aggregate(elevMask, fact=3.6)
elevMask <- resample(elevMask, G0yCMSAF)

G0yCMSAF <- mask(G0yCMSAF, elevMask)
gefCMSAF <- mask(gefCMSAF, elevMask)

FixedCMSAF <- raster(gefCMSAF, layer='Fixed')
HorizCMSAF <- raster(gefCMSAF, layer='Horiz')
TwoCMSAF <- raster(gefCMSAF, layer='Two')

G0yKrig <- mask(G0yKrig, elevMask)
FixedKrig <- mask(FixedKrig, elevMask)
HorizKrig <- mask(HorizKrig, elevMask)
TwoKrig <- mask(TwoKrig, elevMask)

difFixed <- raster(FixedKrig,1)-FixedCMSAF
difTwo <- raster(TwoKrig,1)-TwoCMSAF
difHoriz <- raster(HorizKrig,1)-HorizCMSAF

brickG0y <- stack(G0yCMSAF, raster(G0yKrig, 1))
difG0y <- raster(G0yKrig, 1)- G0yCMSAF

brickFixed <- stack(FixedCMSAF, raster(FixedKrig, 1))
names(brickFixed) <- c('CMSAF', 'SIAR+CMSAF')

difFixed <- raster(FixedKrig, 1)- FixedCMSAF

brickHoriz <- stack(HorizCMSAF, raster(HorizKrig, 1))
names(brickHoriz) <- c('CMSAF', 'SIAR+CMSAF')

difHoriz <- raster(HorizKrig, 1)- HorizCMSAF

brickTwo <- stack(TwoCMSAF, raster(TwoKrig, 1))
names(brickTwo) <- c('CMSAF', 'SIAR+CMSAF')

difTwo <- raster(TwoKrig, 1)- TwoCMSAF

## intervals with equal number of stations
nIntervals <- 8
latEqualCount <- co.intervals(datGef$lat,
                              nIntervals,
                              overlap=0.1)
latEqualCount[1,1] <- ymin(G0yCMSAF)
latEqualCount[nIntervals,2] <- ymax(G0yCMSAF)

nStations <- apply(latEqualCount, 1,
                   function(x)sum(datGef$lat >= x[1] & datGef$lat <= x[2]))

latCutLabels <- apply(latEqualCount, 1, function(x){
  x <- signif(x, 3)
  paste('(', x[1], ',', x[2],']',
        collapse='', sep='')})

latEqualCount <- cbind(latEqualCount, seq_len(nIntervals))
latEqualCountRaster <- reclassify(init(G0yCMSAF, v='y'), latEqualCount)

brickDif <- stack(difG0y/G0yCMSAF * 100,
                  difFixed/FixedCMSAF * 100,
                  difHoriz/HorizCMSAF * 100,
                  difTwo/TwoCMSAF * 100,
                  latEqualCountRaster)
names(brickDif) <- c('G0', 'Fixed', 'Horiz', 'Two', 'latitude')


brickKrig <- stack(raster(G0yKrig, 1),
                 raster(FixedKrig, 1),
                 raster(HorizKrig, 1),
                 raster(TwoKrig, 1))
names(brickKrig) <- c('G0', 'Fixed', 'Horiz', 'Two')


########################
## Statistical summaries
########################
stats <- function(p, o, ...){
  d <- p - o
  meanO <- mean(o, na.rm=TRUE)
  SDp <- sd(p, na.rm=TRUE)
  SDo <- sd(o, na.rm=TRUE)
  MBD <- mean(d, na.rm=TRUE)
  RMSDc <- sd(d, na.rm=TRUE)
  RMSD <- sqrt(MBD^2 + RMSDc^2)
  MAD <- mean(abs(d), na.rm=TRUE)
  res <- data.frame(sdP=SDp, sdO=SDo, MBD=MBD/meanO*100, ##MBDabs=abs(MBD)/meanO*100,
                    RMSDc=RMSDc/meanO*100, RMSD=RMSD/meanO*100, MAD=MAD/meanO*100)
  res
}

datGef <- as.data.frame(spGef)

statPoints <- with(datGef,
                   rbind(stats(G0yCMSAF, G0ySIAR),
                         stats(FixedCMSAF, FixedSIAR),
                         stats(HorizCMSAF, HorizSIAR),
                         stats(TwoCMSAF, TwoSIAR)))
rownames(statPoints) <- c('G0', 'Fixed', 'N-S Horiz', 'Two-axis')
print(xtable(statPoints, digits=2), booktabs=TRUE)

intG0y <- seq(0, .35, by=0.05)
cG0 <- with(datGef, cut(abs(difG0y/G0ySIAR), intG0y))
table(cG0)/length(cG0)*100

statPointsKrig <- with(datGef, rbind(stats(G0yKrig, G0ySIAR),
                                     stats(FixedKrig, FixedSIAR),
                                     stats(HorizKrig, HorizSIAR),
                                     stats(TwoKrig, TwoSIAR)))
rownames(statPointsKrig) <- c('G0', 'Fixed', 'N-S Horiz', 'Two-axis')
print(xtable(statPointsKrig, digits=2), booktabs=TRUE)


statKrigCMSAF <- cellStats(stack(difG0y,
                                 difFixed,
                                 difHoriz,
                                 difTwo),
                           stat=function(d, ...){
                               MBD <- mean(d, na.rm=TRUE)
                               RMSDc <- sd(d, na.rm=TRUE)
                               RMSD <- sqrt(MBD^2 + RMSDc^2)
                               MAD <- mean(abs(d), na.rm=TRUE)
                             data.frame(MBD, RMSDc, RMSD, MAD)
                               })
statKrigCMSAF <- do.call(rbind, statKrigCMSAF)
meanCMSAF <- cellStats(stack(G0yCMSAF, FixedCMSAF, HorizCMSAF, TwoCMSAF), mean, na.rm=1)
sdCMSAF <- cellStats(stack(G0yCMSAF, FixedCMSAF, HorizCMSAF, TwoCMSAF), sd, na.rm=1)
sdKrig <- cellStats(stack(raster(G0yKrig,1),
                            raster(FixedKrig, 1),
                            raster(HorizKrig, 1),
                            raster(TwoKrig, 1)),
                      sd, na.rm=1)
statKrigCMSAF <- sweep(statKrigCMSAF, 1, meanCMSAF, FUN='/')*100
statKrigCMSAF <- cbind(sdKrig, sdCMSAF, statKrigCMSAF)
rownames(statKrigCMSAF) <- c('G0', 'Fixed', 'N-S Horiz', 'Two-axis')

print(xtable(statKrigCMSAF, digits=2), booktabs=TRUE)


################################################################################
## GRAPHICS
################################################################################

##############################
## Comparison SIAR-CMSAF
##############################
meanDif <- zoo(rowMeans(g0Dif/g0SIAR, na.rm=1), index(g0Dif))
medianDif <- zoo(apply(g0Dif/g0SIAR, 1, median, na.rm=1), index(g0Dif))
quantDif <- zoo(t(apply(g0Dif/g0SIAR, 1, quantile, probs=c(0.05, 0.95), na.rm=1)), index(g0Dif))

trellis.device(pdf, file='figs/difG0d.pdf')
p <- xyplot(g0Dif/g0SIAR, superpose=TRUE, 
            lwd=0.2, alpha=0.2, col='black',
            ylim=c(-1, 5),
            ylab=expression(G[d]^{CMSAF}*(0) / G[d]^{SIAR}*(0) -1),##~ ~(Wh/m^2)),
            auto.key=FALSE)


p +
  xyplot(medianDif, lwd=2, col='darkred') +
  xyplot(quantDif, superpose=TRUE, lwd=0.8, col='blue') +
  layer(panel.abline(h=0, col='gray', alpha=0.4))
dev.off()


g0MonthMean <- aggregate(g0SIAR, by=month, mean)
g0RMSDm <- aggregate(g0Dif, by=month, FUN=function(x)sqrt(mean(x^2, na.rm=1)))/g0MonthMean
g0MBDm <- aggregate(g0Dif, by=month, FUN=function(x)mean(x, na.rm=1))/g0MonthMean

spRMSDm <- SpatialPointsDataFrame(coords=spSIAR,
                                  data=as.data.frame(t(coredata(g0RMSDm))))
spMBDm <- SpatialPointsDataFrame(coords=spSIAR,
                                 data=as.data.frame(t(coredata(g0MBDm))))
names(spMBDm) <- names(spRMSDm) <- make.names(index(g0RMSDm))##month.abb

spRMSDm <- SpatialPointsDataFrame(coords=spSIAR,
                                  data=as.data.frame(t(coredata(g0RMSDm))))
spMBDm <- SpatialPointsDataFrame(coords=spSIAR,
                                 data=abs(as.data.frame(t(coredata(g0MBDm)))))
names(spMBDm) <- names(spRMSDm) <- month.abb##make.names(index(g0RMSDm))

trellis.device(pdf, file='figs/bubbleRMSDm.pdf')
bubbles(spRMSDm, n=10, size=c(0.3, 1.1), pwr.size=0.5,
        layout=c(4, 3),
        style='fisher', alpha=1,
        pal=rev(BTC(n=9))) +
  layer_({
    sp.lines(mapaSHP, lwd=0.3)
    sp.polygons(neighbours, col='black', fill='lightgray')
    })
dev.off()

system('pdfcrop figs/bubbleRMSDm.pdf figs/bubbleRMSDm_crop.pdf') ##without margins

trellis.device(pdf, file='figs/bubbleMBDm.pdf')
bubbles(spMBDm, n=10, size=c(0.3, 1.1), pwr.size=0.5,
        layout=c(4, 3),
        style='fisher', alpha=1,
        pal=rev(BTC(n=9))) + ##brewer.pal(name='Blues', n=9)) + 
  layer_({
    sp.lines(mapaSHP, lwd=0.3)
    sp.polygons(neighbours, col='black', fill='lightgray')
    })

dev.off()
system('pdfcrop figs/bubbleMBDm.pdf figs/bubbleMBDm_crop.pdf') ##without margins


trellis.device(pdf, file='figs/bubbleMBDG0.pdf')
bubbles(spSIAR['MBDabs.G0'], n=10, size=c(0.3, 1.1), pwr.size=0.5,
        strip=strip.custom(strip.levels=FALSE),
        style='fisher', alpha=1,
        pal=rev(BTC(n=9))) + ##brewer.pal(name='Blues', n=9)) + 
  layer_({
    sp.lines(mapaSHP, lwd=0.3)
    sp.polygons(neighbours, col='black', fill='lightgray')
    })

dev.off()
system('pdfcrop figs/bubbleMBDG0.pdf figs/bubbleMBDG0_crop.pdf') ##without margins

trellis.device(pdf, file='figs/bubbleMADG0.pdf')
bubbles(spSIAR['MAD.G0'], n=10, size=c(0.3, 1.1), pwr.size=0.5,
        strip=strip.custom(strip.levels=FALSE),
        style='fisher', alpha=1,
        pal=rev(BTC(n=9))) + ##brewer.pal(name='Blues', n=9)) + 
  layer_({
    sp.lines(mapaSHP, lwd=0.3)
    sp.polygons(neighbours, col='black', fill='lightgray')
    })

dev.off()
system('pdfcrop figs/bubbleMADG0.pdf figs/bubbleMADG0_crop.pdf') ##without margins

trellis.device(pdf, file='figs/bubbleRMSDG0.pdf')
bubbles(spSIAR['RMSD.G0'], n=10, size=c(0.3, 1.1), pwr.size=0.5,
        strip=strip.custom(strip.levels=FALSE),
        style='fisher', alpha=1,
        pal=rev(BTC(n=9))) + ##brewer.pal(name='Blues', n=9)) + 
  layer_({
    sp.lines(mapaSHP, lwd=0.3)
    sp.polygons(neighbours, col='black', fill='lightgray')
    })

dev.off()
system('pdfcrop figs/bubbleRMSDG0.pdf figs/bubbleRMSDG0_crop.pdf') ##without margins

trellis.device(pdf, file='figs/bubbleRMSDcG0.pdf')
bubbles(spSIAR['RMSDc.G0'], n=10, size=c(0.3, 1.1), pwr.size=0.5,
        strip=strip.custom(strip.levels=FALSE),
        style='fisher', alpha=1,
        pal=rev(BTC(n=9))) + ##brewer.pal(name='Blues', n=9)) + 
  layer_({
    sp.lines(mapaSHP, lwd=0.3)
    sp.polygons(neighbours, col='black', fill='lightgray')
    })

dev.off()
system('pdfcrop figs/bubbleRMSDcG0.pdf figs/bubbleRMSDcG0_crop.pdf') ##without margins

trellis.device(pdf, file='figs/dotplotDif.pdf')
dotplot(difG0y/G0ySIAR + difFixed/FixedSIAR + difHoriz/HorizSIAR + difTwo/TwoSIAR ~
        reorder(Estacion, difG0y/G0ySIAR),
        data=datGef,
        horizontal=FALSE,
        ylab='Relative differences between SIAR and CMSAF',
        xlab='',
        auto.key=list(corner=c(0.9, 0.1),
          columns=2,
          cex=0.9,
          text=c(expression({G[y]^{CMSAF}}(0)/{G[y]^{SIAR}}(0) - 1),
            expression({G[fixed]^{CMSAF}}/{G[fixed]^{SIAR}} - 1),
            expression({G[horiz]^{CMSAF}}/{G[horiz]^{SIAR}} - 1),
            expression({G[two]^{CMSAF}}/{G[two]^{SIAR}} - 1))),
        alpha=0.6,
        ##scales=list(x=list(rot=90, cex=0.3)),
        scales=list(x=list(draw=FALSE)),
        panel=function(...){
          panel.abline(h=c(-0.05, 0, 0.05), col.line='darkgray')
          panel.dotplot(...)
        })
dev.off()

trellis.device(pdf, file='figs/dotplotDifKrig.pdf')
dotplot(difG0yKrig/G0ySIAR + difFixedKrig/FixedSIAR + difHorizKrig/HorizSIAR + difTwoKrig/TwoSIAR ~
        reorder(Estacion, difG0yKrig/G0ySIAR),
        data=datGef,
        horizontal=FALSE,
        ylab='Relative differences between SIAR and KED',
        xlab='',
        auto.key=list(corner=c(0.9, 0.1),
          columns=2,
          cex=0.9,
          text=c(expression({G[y]^{KED}}(0)/{G[y]^{SIAR}}(0) - 1),
            expression({G[fixed]^{KED}}/{G[fixed]^{SIAR}} - 1),
            expression({G[horiz]^{KED}}/{G[horiz]^{SIAR}} - 1),
            expression({G[two]^{KED}}/{G[two]^{SIAR}} - 1))),
        alpha=0.6,
        ##scales=list(x=list(rot=90, cex=0.3)),
        scales=list(x=list(draw=FALSE)),
        panel=function(...){
          panel.abline(h=c(-0.05, 0, 0.05), col.line='darkgray')
          panel.dotplot(...)
        })
dev.off()


spDif <- spGef[c('difG0y', 'difFixed', 'difHoriz', 'difTwo')]
spDif@data <- spDif@data/spGef[,c('G0ySIAR', 'FixedSIAR', 'HorizSIAR', 'TwoSIAR')]@data
names(spDif) <- c('G0', 'Fixed', 'N-S Horiz', 'Two axis')

trellis.device(pdf, file='figs/bubbleDif.pdf')
bubbles(spDif, n=8, size=c(0.3, 1.1), pwr.size=0.5,
        style='fisher', alpha=1,
        pal=rev(BTC(n=9))) + ##brewer.pal(name='Blues', n=9)) + 
  layer_({
    sp.lines(mapaSHP, lwd=0.3)
    sp.polygons(neighbours, col='black', fill='lightgray')
    })

dev.off()
system('pdfcrop figs/bubbleDif.pdf figs/bubbleDif_crop.pdf') ##without margins

spDifKrig <- spGef[c('difG0yKrig', 'difFixedKrig', 'difHorizKrig', 'difTwoKrig')]
spDifKrig@data <- spDifKrig@data/datGef[,c('G0ySIAR', 'FixedSIAR', 'HorizSIAR', 'TwoSIAR')]
names(spDifKrig) <- c('G0', 'Fixed', 'N-S Horiz', 'Two axis')

trellis.device(pdf, file='figs/bubbleDifKrig.pdf')
bubbles(spDifKrig, n=8, size=c(0.3, 1.1), pwr.size=0.5,
        style='fisher', alpha=1,
        pal=rev(BTC(n=9))) + ##brewer.pal(name='Blues', n=9)) + 
  layer_({
    sp.lines(mapaSHP, lwd=0.3)
    sp.polygons(neighbours, col='black', fill='lightgray')
    })
dev.off()
system('pdfcrop figs/bubbleDifKrig.pdf figs/bubbleDifKrig_crop.pdf') ##without margins


##########
## MAPA SHP
##########

##Plot the stations in a map
p <- spplot(spSIAR['Comunidad'],
       col.regions=brewer.pal(n=12, 'Paired'),  
       key.space='right', scales=list(draw=TRUE))

terrainTheme <- rasterTheme(region=terrain.colors(15))
terrainTheme$panel.background$col <- 'lightblue' ##Sea

trellis.device(pdf, file='figs/mapaSIAR.pdf')
levelplot(elevMask, par.settings=terrainTheme,
          panel=panel.levelplot.raster, interpolate=TRUE,
          colorkey=list(raster=TRUE, interpolate=TRUE),
          margin=FALSE) +
  layer(sp.polygons(neighbours, col='black', fill='lightgray')) +
  layer(sp.points(spSIAR[spSIAR$outTolerance,], pch=16, col='red', cex=0.4, alpha=0.8)) +
  layer(sp.points(spSIAR[!spSIAR$outTolerance,], pch=16, col='darkblue', cex=0.4, alpha=0.8)) + 
  layer(sp.lines(mapaSHP, lwd=0.3))
dev.off()


system('pdfcrop figs/mapaSIAR.pdf figs/mapaSIAR_crop.pdf') ##without margins

##########
## G0
##########

radTheme <- modifyList(rasterTheme(),
                       list(panel.background=list(col=adjustcolor('lightskyblue1', alpha=0.3))))
trellis.device(pdf, file='figs/G0yCMSAF.pdf')
levelplot(G0yCMSAF, panel=panel.levelplot.raster, par.settings=radTheme) +
  layer(sp.polygons(neighbours, col='black', fill='lightgray')) +
  layer(sp.points(spSIAR, pch=19, cex=0.3, col='black')) +
  layer(sp.lines(mapaSHP, lwd=0.5))
dev.off()

trellis.device(pdf, file='figs/difG0y.pdf')
xyplot(difG0y/G0ySIAR~lat, data=datGef,
       ylim=c(-0.25, 0.45),
       groups=outTolerance, auto.key=FALSE,
       xlab='Latitude',
       ylab=expression(frac({G[y]^{CMSAF}}(0),{G[y]^{SIAR}}(0)) - 1),
       panel=function(...){
         panel.abline(h=c(-0.05, 0.05), col='darkgray')
         panel.xyplot(...)
       })
dev.off()

##################################################
#### ANALISIS ESTADÍSTICO RADIACION HORIZONTAL
##################################################

trellis.device(pdf, file='figs/difG0yKrigSIAR.pdf')
xyplot(difG0yKrig/G0ySIAR~lat, data=datGef,
       ylim=c(-0.25, 0.35),
       groups=outTolerance,
       auto.key=FALSE,
       xlab='Latitude',
       ylab=expression(frac({G[y]^{KED}}(0),{G[y]^{SIAR}}(0)) - 1),
       panel=function(...){
         panel.abline(h=c(-0.05, 0.05), col='darkgray')
         panel.xyplot(...)
       })
dev.off()


trellis.device(pdf, file='figs/G0yKrig.pdf')
levelplot(G0yKrig, layer='pred', panel=panel.levelplot.raster, par.settings=radTheme) +
  layer(sp.polygons(neighbours, col='black', fill='lightgray')) +
  layer(sp.points(spGef, pch=19, cex=0.3, col='black')) +
  layer(sp.lines(mapaSHP))
dev.off()

RdBuTheme <- modifyList(RdBuTheme(),
                        list(panel.background=list(col=adjustcolor('lightskyblue1', alpha=0.3))))

## Same breaks for all plots
minDif <- min(minValue(brickDif)[1:4])
maxDif <- max(maxValue(brickDif)[1:4])
maxAbsDif <- max(abs(c(minDif, maxDif)))
rangeDif <- c(-maxAbsDif, maxAbsDif)
atDif <- pretty(rangeDif, 15)
## rangeG0y <- c(minValue(difG0y/G0yCMSAF), maxValue(difG0y/G0yCMSAF))
## maxG0y <- max(abs(rangeG0y))
## rangeG0y <- c(-maxG0y, maxG0y)

trellis.device(pdf, file='figs/difG0yKrigCMSAF.pdf')
levelplot(brickDif, layers='G0', par.settings=RdBuTheme,
          at=atDif, panel=panel.levelplot.raster) +
  layer(sp.polygons(neighbours, col='black', fill='lightgray')) +
  layer(sp.points(spGef, pch=19, cex=0.3, col='black')) +
  layer(sp.lines(mapaSHP))
dev.off()

trellis.device(pdf, file='figs/difG0yKrigCMSAFBWplot.pdf')
bwplot(G0~latitude, data=brickDif,
       xlab='Latitude', ylab='',
       scales=list(x=list(rot=30,
                     labels=latCutLabels)))
dev.off()

##################################################
#### ANALISIS ESTADÍSTICO RADIACION EFECTIVA
##################################################


#### Fixed
trellis.device(pdf, file='figs/difFixedKrigSIAR.pdf')
xyplot(difFixedKrig/FixedSIAR~lat, data=datGef,
       ylim=c(-0.25, 0.35),
       groups=outTolerance,
       auto.key=FALSE,
       xlab='Latitude',
       ylab=expression(frac({G[efy]^{KED}},{G[efy]^{SIAR}}) - 1),
       panel=function(...){
         panel.abline(h=c(-0.05, 0.05), col='darkgray')
         panel.xyplot(...)
       })
dev.off()


trellis.device(pdf, file='figs/difFixed.pdf')
xyplot(difFixed/FixedSIAR~lat, data=datGef,
       ylim=c(-0.25, 0.45),
       groups=outTolerance, auto.key=FALSE,
       xlab='Latitude',
       ylab=expression(frac({G[efy]^{CMSAF}},{G[efy]^{SIAR}}) - 1),
       panel=function(...){
         panel.abline(h=c(-0.05, 0.05), col='darkgray')
         panel.xyplot(...)
         })
dev.off()

levelplot(brickFixed) +
  layer(sp.polygons(neighbours, col='black', fill='lightgray')) +
  layer(sp.points(spGef, pch=19, cex=0.3, col='black')) +
  layer(sp.lines(mapaSHP))

## rangeFixed <- c(minValue(difFixed/FixedCMSAF), maxValue(difFixed/FixedCMSAF))
## maxFixed <- max(abs(rangeFixed))
## rangeFixed <- c(-maxFixed, maxFixed)

trellis.device(pdf, file='figs/difFixedKrigCMSAF.pdf')
levelplot(brickDif, layers='Fixed', par.settings=RdBuTheme,
          at=atDif, panel=panel.levelplot.raster) +
  layer(sp.polygons(neighbours, col='black', fill='lightgray')) +
  layer(sp.points(spGef, pch=19, cex=0.3, col='black')) +
  layer(sp.lines(mapaSHP))
dev.off()

trellis.device(pdf, file='figs/difFixedKrigCMSAFBWplot.pdf')
bwplot(Fixed~latitude, data=brickDif,
       xlab='Latitude', ylab='',
       scales=list(x=list(rot=30,
                     labels=latCutLabels)))
dev.off()


#### Horiz
trellis.device(pdf, file='figs/difHorizKrigSIAR.pdf')
xyplot(difHorizKrig/HorizSIAR~lat, data=datGef,
       ylim=c(-0.25, 0.35),
       groups=outTolerance,
       auto.key=FALSE,
       xlab='Latitude',
       ylab=expression(frac({G[efy]^{KED}},{G[efy]^{SIAR}}) - 1),
       panel=function(...){
         panel.abline(h=c(-0.05, 0.05), col='darkgray')
         panel.xyplot(...)
       })
dev.off()

trellis.device(pdf, file='figs/difHoriz.pdf')
xyplot(difHoriz/HorizSIAR~lat, data=datGef,
       ylim=c(-0.25, 0.45),
       groups=outTolerance, auto.key=FALSE,
       xlab='Latitude',
       ylab=expression(frac({G[efy]^{CMSAF}},{G[efy]^{SIAR}}) - 1),
              panel=function(...){
         panel.abline(h=c(-0.05, 0.05), col='darkgray')
         panel.xyplot(...)})
dev.off()

levelplot(brickHoriz) +
  layer(sp.polygons(neighbours, col='black', fill='lightgray')) +
  layer(sp.points(spGef, pch=19, cex=0.3, col='black')) +
  layer(sp.lines(mapaSHP))

## rangeHoriz <- c(minValue(difHoriz/HorizCMSAF), maxValue(difHoriz/HorizCMSAF))
## maxHoriz <- max(abs(rangeHoriz))
## rangeHoriz <- c(-maxHoriz, maxHoriz)

trellis.device(pdf, file='figs/difHorizKrigCMSAF.pdf')
levelplot(brickDif, layers='Horiz', par.settings=RdBuTheme,
          at=atDif, panel=panel.levelplot.raster) +
  layer(sp.polygons(neighbours, col='black', fill='lightgray')) +
  layer(sp.points(spGef, pch=19, cex=0.3, col='black')) +
  layer(sp.lines(mapaSHP))
dev.off()

trellis.device(pdf, file='figs/difHorizKrigCMSAFBWplot.pdf')
bwplot(Horiz~latitude, data=brickDif,
       xlab='Latitude', ylab='',
       scales=list(x=list(rot=30,
                     labels=latCutLabels)))
dev.off()

#### Two
trellis.device(pdf, file='figs/difTwoKrigSIAR.pdf')
xyplot(difTwoKrig/TwoSIAR~lat, data=datGef,
       ylim=c(-0.25, 0.35),
       groups=outTolerance,
       auto.key=FALSE,
       xlab='Latitude',
       ylab=expression(frac({G[efy]^{KED}},{G[efy]^{SIAR}}) - 1),
       panel=function(...){
         panel.abline(h=c(-0.05, 0.05), col='darkgray')
         panel.xyplot(...)
       })
dev.off()

trellis.device(pdf, file='figs/difTwo.pdf')
xyplot(difTwo/TwoSIAR~lat, data=datGef,
       ylim=c(-0.25, 0.45),
       groups=outTolerance, auto.key=FALSE,
       xlab='Latitude',
       ylab=expression(frac({G[efy]^{CMSAF}},{G[efy]^{SIAR}}(0)) - 1),
              panel=function(...){
         panel.abline(h=c(-0.05, 0.05), col='darkgray')
         panel.xyplot(...)})
dev.off()

levelplot(brickTwo) +
  layer(sp.polygons(neighbours, col='black', fill='lightgray')) +
  layer(sp.points(spGef, pch=19, cex=0.3, col='black')) +
  layer(sp.lines(mapaSHP))

## rangeTwo <- c(minValue(difTwo/TwoCMSAF), maxValue(difTwo/TwoCMSAF))
## maxTwo <- max(abs(rangeTwo))
## rangeTwo <- c(-maxTwo, maxTwo)

trellis.device(pdf, file='figs/difTwoKrigCMSAF.pdf')
levelplot(brickDif, layers='Two', par.settings=RdBuTheme,
          at=atDif, panel=panel.levelplot.raster) +
  layer(sp.polygons(neighbours, col='black', fill='lightgray')) +
  layer(sp.points(spGef, pch=19, cex=0.3, col='black')) +
  layer(sp.lines(mapaSHP))
dev.off()

trellis.device(pdf, file='figs/difTwoKrigCMSAFBWplot.pdf')
bwplot(Two~latitude, data=brickDif,
       xlab='Latitude', ylab='',
       scales=list(x=list(rot=30,
                     labels=latCutLabels)))
dev.off()

trellis.device(pdf, file='figs/difKrigCMSAFBWplot.pdf', width=10)
bwplot(G0 + Fixed + Horiz + Two ~ latitude, data=brickDif,
       outer=TRUE, layout=c(4, 1),
       xlab='Latitude', ylab='',
       strip=strip.custom(factor.levels =
         c('Horizontal Irradiation', 'Fixed Plane',
           'NS Hor. Axis Tracker', 'Two Axis'),
         par.strip.text=list(cex=0.8)),
       scales=list(x=list(rot=45, cex=0.6,
                     labels=latCutLabels)))
dev.off()
