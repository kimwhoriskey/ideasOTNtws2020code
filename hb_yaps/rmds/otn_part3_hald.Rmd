---
title: YAPS - Yet Another Positioning Solver
author:
  affiliation: Technical University of Denmark (DTU)
  email: hba@aqua.dtu.dk
  name: Henrik Baktoft
date: Feb 19, 2020
output:
  html_notebook:
    toc: yes
    toc_float: no
  pdf_document:
    toc: yes
  html_document:
    df_print: paged
    toc: yes
subtitle: |
  | OTN Telemetry Workshop Series, Dalhousie University, Halifax, Canada
  | Part 3 - Danish pike in a large lake
---
<div style="float: right;">
  [![][yaps_logo]](https://github.com/baktoft/yaps)   ![][otn_logo]  
</div>


# Setup and needed libraries
Set timezone to UTC and load needed libraries
```{r, echo=TRUE, message=FALSE}
rm(list=ls())
Sys.setenv(TZ='UTC')
library(data.table)
library(leaflet)
library(sp)
library(lubridate)
library(ggplot2)
library(caTools)
library(viridis)
library(yaps)
source('../yaps_OTN_helperFuncs.R')
```

# Description
These data are included in the `yaps` package (run `?hald` for info). Data were collected as part of a larger study including brown trout (smolt and adult), pike and eel in the Danish lake Hald. Lake Hald is large and deep (by Danish standards; 3.42 km^2^;  max depth 34 m), so depth is relevant in order to get best possible track estimation. A total of 70 hydrophones (Thelma TBR700) were deployed to cover the entire lake. The data set contains a single week of detections including a single tow track and tracks of a few pike (*Esox lucius*).  

This example focus on:

* Tag 5138: tow track, fixed BI (10 s), **gps data available**
* Tag 1335: pike, random BI (10-30 s), **known sequence**
* Tag 1315: pike, random BI (60-120 s), **known sequence**

An appropiate sync_model is included in the package and will be used for this demonstration to save time. Code to obtain this sync_model can be found in the appendix below.  

## A look at the data
Notice the object `burst_seqs` - these are vectors of intervals between consequtive pings from the two transmitters included in this example data. Not all manufacturers are willing to provide these data to their users even though it is just a stream of random numbers - luckily Thelma Biotel are.  
`data.table`s `hydros` and `detections` are as usual.  
The object `sync_model` contains the already optimized synchronization model applicable for these data.
```{r}
names(hald)
hydros <- hald$hydros
detections <- hald$detections
sync_model <- hald$sync_model
str(hald$burst_seqs)
```


## The study site
Study site geometry is rather complex. (Code to make the map is in the .Rmd file).  
Note, how a sequential synchronization process in this case could lead to error propagation.
Also note, there are lots of impossible sync tag - hydro combinations.
```{r echo=FALSE, eval=TRUE}
hydros_map <- hald$hydros
coordinates(hydros_map) <- ~x+y
proj4string(hydros_map) <- CRS("+init=epsg:32632")
hydros_latlon <- spTransform(hydros_map, CRS("+init=epsg:4326"))
m <- leaflet(data=hydros_latlon, options = leafletOptions(minZoom = 0, maxZoom = 18), width="100%", height=700)
m <- addTiles(m, group="OSM")
m <- addCircles(m, radius=5, label=as.character(hydros_latlon$idx), labelOptions = labelOptions(noHide = T, textOnly = TRUE))
m <- addCircles(m, data=hydros_latlon[!is.na(hydros_latlon$sync_tag), ], col="red", radius=2)
m <- addMeasure(m, primaryLengthUnit="meters")
m <- addProviderTiles(m, providers$Esri.WorldImagery, group="Esri.WorldImagery")
m <- addLayersControl(m, baseGroups = c("OSM (default)", "Esri.WorldImagery"),    options = layersControlOptions(collapsed = FALSE)  )
m
```


## A look at the sync model
There might be room for improvement, but it will suffice for now. Additionally, the hope for this study was not to achieve sub-meter accuracy and precision given the large extent of the study area.
```{r}
plotSyncModelResids(sync_model, by='overall')
plotSyncModelResids(sync_model, by='quantiles')
# plotSyncModelResids(sync_model, by='sync_tag')
# plotSyncModelResids(sync_model, by='hydro')
```


As part of fitting the synchronization model, we estimated positions of two hydros as their gps-obtained positions were unceratin. We will use these estimated positions for the track estimation.
```{r}
which(sync_model$inp_synced$dat_tmb_sync$fixed_hydros_vec !=1)

hydros_yaps <- data.table::data.table(sync_model$pl$TRUE_H)
colnames(hydros_yaps) <- c('hx','hy','hz')

cbind(hydros[c(32, 70), c('x','y')], hydros_yaps[c(32,70), c('hx','hy')])
```


## Applying the sync model to all detections
This is easily done, since we already have the sync_model
```{r}
detections_synced <- applySync(toa=detections, hydros=hydros, sync_model=sync_model)
```


# Estimating tracks using `yaps`
## Tag 5138 - tow track - **fixed** burst interval
### Compiling the data
Remember this is a fixed burst interval transmitter (BI = 10 s).  

Let's have a look at the data. Note, this plot only makes sense for transmitters where the burst interval is set at whole seconds (random or not).
```{r}
dat_tow <- detections_synced[tag==5138]
gps <- fread("../data/hald_gps_test_track.csv")

with(dat_tow, table(hydro_idx))
ggplot(data=dat_tow) + geom_point(aes(x=eposync, y=eposync-floor(eposync), col=factor(hydro_idx)))

```


Get TOA matrix for YAPS and have a look
```{r}
toa_tow <- getToaYaps(synced_dat=dat_tow, hydros=hydros, rbi_min=10, rbi_max=10)
plotToa(toa_tow)
plotNobs(toa_tow)
```



### Get `inp` and run YAPS
```{r}
inp_tow <- getInp(toa=toa_tow, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType="sbi", sdInits=1, ss_data_what="est", ss_data=0, biTable=NULL)
yaps_tow <- runYaps(inp_tow, maxIter=500, getPlsd=TRUE, getRep=TRUE, silent=TRUE)
```

### Results
```{r}
plotYaps(inp_tow, yaps_tow)
lines(Y~X, data=gps, lty=2)

# zoom in 
plotYaps(inp_tow, yaps_tow, xlim=c(520500, 521500), ylim=c(6247400, 6247800))
lines(Y~X, data=gps, lty=2)

plotYaps(inp_tow, yaps_tow, type="coord_X")
plotYaps(inp_tow, yaps_tow, type="coord_Y")

# standard error of position estimates are hardly visible
plot(yaps_tow$plsd$X, type="l")
plot(yaps_tow$plsd$Y, type="l")

# Fixed burst interval - how fixed? Note y-axis is in meter
bi <- diff(yaps_tow$pl$top)
plot((bi - mean(bi)) * 1460, ylab="BI variation (m)", type="h")

# estimated speed of sound - sanity check
yaps_tow$pl$ss
getSS(temp=c(15, 16, 17))

# look at the residual
matplot(t((yaps_tow$rep$mu_toa - inp_tow$datTmb$toa))*1460, ylab="TOA residuals (m)")

```


### Tag 5138 tow track - downsampling max_hydro <- 3
Looks good, but this data set is highly saturated! To showcase the benefit of fixed burst interval, we downsample and re-run `yaps`
```{r}
plotNobs(toa_tow)

max_hydro <- 3
toa_tow3<- downSampleToa(toa_tow, max_hydro)

inp_tow3 <- getInp(toa=toa_tow3, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType="sbi", sdInits=1, ss_data_what="est", ss_data=0, biTable=NULL)
yaps_tow3 <- runYaps(inp_tow3, maxIter=500, getPlsd=TRUE, getRep=TRUE, silent=TRUE)

plotYaps(inp_tow3, yaps_tow3, xlim=c(520500, 521500), ylim=c(6247400, 6247800))
lines(Y~X, data=gps, lty=2)

plotYaps(inp_tow3, yaps_tow3, type="coord_X")
plotYaps(inp_tow3, yaps_tow3, type="coord_Y")

```

### Tag 5138 tow track - downsampling max_hydro <- 1
For fixed burst interval transmitters, we can go really low on detections.
```{r}
plotNobs(toa_tow)

max_hydro <- 1
toa_tow1<- downSampleToa(toa_tow, max_hydro)

inp_tow1 <- getInp(toa=toa_tow1, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType="sbi", sdInits=1, ss_data_what="est", ss_data=0, biTable=NULL)
yaps_tow1 <- runYaps(inp_tow1, maxIter=500, getPlsd=TRUE, getRep=TRUE, silent=TRUE)

plotYaps(inp_tow1, yaps_tow1, xlim=c(520500, 521500), ylim=c(6247400, 6247800))
lines(Y~X, data=gps, lty=2)

plotYaps(inp_tow1, yaps_tow1, type="coord_X")
plotYaps(inp_tow1, yaps_tow1, type="coord_Y")

```

### Comparing standard errors of position estimates
```{r}

par(mfrow=c(2,2))
plotYaps(inp_tow, yaps_tow, type="coord_X", main="X all data")
plotYaps(inp_tow, yaps_tow, type="coord_Y", main="Y all data")
plotYaps(inp_tow1, yaps_tow1, type="coord_X", main="X max_hydro=1")
plotYaps(inp_tow1, yaps_tow1, type="coord_Y", main="Y max_hydro=1")
par(mfrow=c(1,1))

par(mfrow=c(3,1))
ylim <- c(0, 5)
plot(yaps_tow$plsd$X,  type="h", ylim=ylim, main="Full")
plot(yaps_tow3$plsd$X, type="h", ylim=ylim, main="max_hydro=3")
plot(yaps_tow1$plsd$X, type="h", ylim=ylim, main="max_hydro=1")
par(mfrow=c(1,1))

```

### Note on the downsampling
One reason we can estimate a decent track even with max_hydro = 1 is, that multiple hydrophones are detecting the signals. It will, of course, not work if there only was a single hydrophone in the water. 


## Tag 1335 - pike - random BI sequence - assuming **secret** sequence
Remember this transmitter is a random burst interval (BI 10 - 30 s), but we know the burst sequence. However, first we assume, we don't and treat it like a random BI trasnmitter with secret burst sequence.
### Compiling the data
```{r}
dat1335 <- detections_synced[tag==1335]
ping_type <- 'rbi'
rbi_min <- 10
rbi_max <- 30

# get number of detections per hydro per hour
dat1335[, hour := floor_date(ts, unit="hour")]
summ <- dat1335[, .N, by=c('hydro_idx', 'hour')]
ggplot(data=summ) + geom_tile(aes(x=hour, y=hydro_idx, fill=N)) + scale_fill_viridis()

# Get TOA matrix for YAPS
toa1335_rbi <- getToaYaps(synced_dat=dat1335, hydros=hydros_yaps, rbi_min=rbi_min, rbi_max=rbi_max)

# Have a look at the raw data...
plotToa(toa1335_rbi)
plotNobs(toa1335_rbi)

```
More than 30.000 pings will take pretty long time to process, so we use a small a subset of the data.
```{r}
ping_1 <- 501
n_ping <- 1000
toa1335_rbi_yaps <- toa1335_rbi[ping_1:(ping_1 + n_ping - 1), ]
plotToa(toa1335_rbi_yaps)
plotNobs(toa1335_rbi_yaps)

```



### Get inp and run YAPS
```{r}
ping_type <- 'rbi'
inp1335_rbi <- getInp(toa=toa1335_rbi_yaps, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType=ping_type, rbi_min=rbi_min, rbi_max=rbi_max, sdInits=1, ss_data_what="est", ss_data=0, biTable=NULL)
inp1335_rbi$inits[1] <- rnorm(1, 1, 1)
yaps1335_rbi <- runYaps(inp1335_rbi, maxIter=500, getPlsd=TRUE, getRep=TRUE, silent=TRUE)

```

### Results
```{r}
# basic plotting
plotYaps(inp1335_rbi, yaps1335_rbi)

# zoom in 
plotYaps(inp1335_rbi, yaps1335_rbi, xlim=c(520500, 521500), ylim=c(6247400, 6247800))

plotYaps(inp1335_rbi, yaps1335_rbi, type="coord_X")
plotYaps(inp1335_rbi, yaps1335_rbi, type="coord_Y")


# uncertainties... standard error of position estimates
plot(yaps1335_rbi$plsd$X, type="h")
plot(yaps1335_rbi$plsd$Y, type="h")

# estimated speed of sound - sanity check
yaps1335_rbi$pl$ss
getSS(temp=c(15, 17))

# look at the residual
matplot(t((yaps1335_rbi$rep$mu_toa - inp1335_rbi$datTmb$toa))*1460, ylab="TOA residuals (m)")

```


## Tag 1335 - pike - random BI - utilizing **known** sequence
To utilize the known sequence of burst intervals, we need to align the observed data with the sequence. For this, we can use the function `alignBurstSeq()` from `yaps`. We can also take advantage of knowing the burst sequence when compiling the TOA-marix and we need extra info to feed into `getInp()`.

### Compiling the data
```{r}
dat1335 <- detections_synced[tag==1335]
ping_type <- 'pbi'
rbi_min <- 10
rbi_max <- 30

seq1335 <- hald$burst_seqs$seq1335 # just a long sequence of random numbers...
length(seq1335)

dat1335_kbi <- alignBurstSeq(synced_dat=dat1335, burst_seq=seq1335, seq_lng_min=25, rbi_min=rbi_min, rbi_max=rbi_max, plot_diag=TRUE)

toa1335_kbi_list <- buildToaKnownSeq(seq=seq1335, aligned_dat=dat1335_kbi, hydros=hydros_yaps)

toa1335_kbi <- toa1335_kbi_list$toa
seq1335_kbi <- toa1335_kbi_list$seq

```



More than 30 K pings will still take a long time to finish - we use the same subset as above.
```{r}
ping_1 <- 501
n_ping <- 1000
toa1335_kbi_yaps <- toa1335_kbi[ping_1:(ping_1 + n_ping - 1), ]
biTable1335 <- seq1335_kbi[ping_1:(ping_1 + n_ping - 1) ]
plotToa(toa1335_kbi_yaps)
plotNobs(toa1335_kbi_yaps)

```

### Get inp and run YAPS
Note, that ping_type now is `'pbi'`and we include the BItable in `getInp()`. Everything else is the same as above.
```{r}
inp1335_kbi <- getInp(toa=toa1335_kbi_yaps, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType=ping_type, rbi_min=rbi_min, rbi_max=rbi_max, sdInits=1, ss_data_what="est", ss_data=0, biTable=biTable1335)

yaps1335_kbi <- runYaps(inp1335_kbi, maxIter=500, getPlsd=TRUE, getRep=TRUE, silent=TRUE)

```
### Results
Quick look at the results
```{r}
# basic plotting
plotYaps(inp1335_kbi, yaps1335_kbi)

# zoom in 
plotYaps(inp1335_kbi, yaps1335_kbi, xlim=c(520500, 521500), ylim=c(6247400, 6247800))

plotYaps(inp1335_kbi, yaps1335_kbi, type="coord_X")
plotYaps(inp1335_kbi, yaps1335_kbi, type="coord_Y")


# uncertainties... standard error of position estimates
plot(yaps1335_kbi$plsd$X, type="h")
plot(yaps1335_kbi$plsd$Y, type="h")

# estimated speed of sound - sanity check
yaps1335_kbi$pl$ss
getSS(temp=c(15, 17))

# look at the residual
matplot(t((yaps1335_kbi$rep$mu_toa - inp1335_kbi$datTmb$toa))*1460, ylab="TOA residuals (m)")

```


### Tag 1335 pike - downsampling max_hydro <- 3
Knowing the burst interval sequence, gives opportunity to estimate tracks even when data density is low.
```{r}
ping_type <- 'pbi'

plotNobs(toa1335_kbi_yaps)
max_hydro <- 3
toa1335_kbi_down3 <- downSampleToa(toa1335_kbi_yaps, max_hydro)

inp1335_kbi_down3 <- getInp(toa=toa1335_kbi_down3, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType=ping_type, rbi_min=rbi_min, rbi_max=rbi_max, sdInits=1, ss_data_what="est", ss_data=0, biTable=biTable1335)

yaps1335_kbi_down3 <- runYaps(inp1335_kbi_down3, maxIter=500, getPlsd=TRUE, getRep=TRUE, silent=TRUE)

plotYaps(inp1335_kbi_down3, yaps1335_kbi_down3)
# zoom in 
plotYaps(inp1335_kbi_down3, yaps1335_kbi_down3, xlim=c(520500, 521500), ylim=c(6247400, 6247800))

plotYaps(inp1335_kbi_down3, yaps1335_kbi_down3, type="coord_X")
plotYaps(inp1335_kbi_down3, yaps1335_kbi_down3, type="coord_Y")
```
### Tag 1335 pike - downsampling max_hydro <- 1
Trying at very low data density
```{r}
ping_type <- 'pbi'

plotNobs(toa1335_kbi_yaps)
max_hydro <- 1
toa1335_kbi_down1 <- downSampleToa(toa1335_kbi_yaps, max_hydro)

inp1335_kbi_down1 <- getInp(toa=toa1335_kbi_down1, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType=ping_type, rbi_min=rbi_min, rbi_max=rbi_max, sdInits=1, ss_data_what="est", ss_data=0, biTable=biTable1335)

yaps1335_kbi_down1 <- runYaps(inp1335_kbi_down1, maxIter=500, getPlsd=TRUE, getRep=TRUE, silent=TRUE)

# zoom in 
plotYaps(inp1335_kbi_down1, yaps1335_kbi_down1, xlim=c(520500, 521500), ylim=c(6247400, 6247800))

plotYaps(inp1335_kbi_down1, yaps1335_kbi_down1, type="coord_X")
plotYaps(inp1335_kbi_down1, yaps1335_kbi_down1, type="coord_Y")
```

### Comparing secret vs known BI sequences
```{r}
par(mfrow=c(2,3))
ylim=c(0,5)
plot(yaps1335_kbi$plsd$X, type="h", ylim=ylim)
plot(yaps1335_kbi_down1$plsd$X, type="h", ylim=ylim)
plot(yaps1335_rbi$plsd$X, type="h", ylim=ylim)
plot(yaps1335_kbi$plsd$Y, type="h", ylim=ylim)
plot(yaps1335_kbi_down1$plsd$Y, type="h", ylim=ylim)
plot(yaps1335_rbi$plsd$Y, type="h", ylim=ylim)
par(mfrow=c(1,1))

```

### Tag 1335 - the full track
The full track can be estimated using this code:
```{r eval=FALSE}
ping_type <- 'pbi'
rbi_min <- 10
rbi_max <- 30
seq1335 <- hald$burst_seqs$seq1335 # just a long sequence of random numbers...
dat1335_kbi <- alignBurstSeq(synced_dat=dat1335, burst_seq=seq1335, seq_lng_min=10, rbi_min=rbi_min, rbi_max=rbi_max, plot_diag=TRUE)
toa1335_kbi_list <- buildToaKnownSeq(seq=seq1335, aligned_dat=dat1335_kbi, hydros=hydros_yaps)
toa1335_kbi <- toa1335_kbi_list$toa
seq1335_kbi <- toa1335_kbi_list$seq
biTable1335 <- seq1335_kbi

ping_1 <- 1
n_ping <- nrow(toa1335_kbi)
toa1335_kbi_yaps <- toa1335_kbi[ping_1:(ping_1 + n_ping - 1), ]
biTable1335 <- seq1335_kbi[ping_1:(ping_1 + n_ping - 1) ]
plotToa(toa1335_kbi_yaps)
plotNobs(toa1335_kbi_yaps)

inp1335_kbi <- getInp(toa=toa1335_kbi_yaps, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType=ping_type, rbi_min=rbi_min, rbi_max=rbi_max, sdInits=1, ss_data_what="est", ss_data=0, biTable=biTable1335)
inp1335_kbi$inits[1] <- 1
yaps1335_kbi <- runYaps(inp1335_kbi, maxIter=1000, getPlsd=TRUE, getRep=TRUE, silent=FALSE)
```


Lazy people can just grab it here:
```{r}
load('../data/yaps1335_kbi_fullTrack.rObj')
plotYaps(inp1335_kbi, yaps1335_kbi)

plotYaps(inp1335_kbi, yaps1335_kbi, type="coord_X")
plotYaps(inp1335_kbi, yaps1335_kbi, type="coord_Y")

```


## Tag 1315 - pike - random BI - assuming **secret** sequence
This transmitter is a random BI 60 - 120 s.
```{r}
dat1315 <- detections_synced[tag==1315]
rbi_min <- 60
rbi_max <- 120
dat1315[, hour := floor_date(ts, unit="hour")]
summ <- dat1315[, .N, by=c('hydro_idx', 'hour')]
ggplot(data=summ) + geom_tile(aes(x=hour, y=hydro_idx, fill=N)) + scale_fill_viridis()

```

```{r}
# First, assume burst sequence is unknown
# Get TOA matrix for YAPS
ping_type <- 'rbi'
toa1315_rbi <- getToaYaps(synced_dat=dat1315, hydros=hydros_yaps, rbi_min=rbi_min, rbi_max=rbi_max)

# Have a look at the raw data...
plotToa(toa1315_rbi)
plotNobs(toa1315_rbi)
nobs <- getNobs(toa1315_rbi)
sum(nobs>=3) / length(nobs)
```

```{r}
# < 10.000 pings is doable, but we will sub sample to save time
ping_1 <- 1
n_ping <- 1000
toa1315_rbi_yaps <- toa1315_rbi[ping_1:(ping_1 + n_ping - 1), ]
plotToa(toa1315_rbi_yaps)
plotNobs(toa1315_rbi_yaps)
nobs <- getNobs(toa1315_rbi_yaps)
sum(nobs>=3) / length(nobs)

# looks a bit more challenging than 1335


```
```{r}
# get inp to feed into YAPS and run YAPS
inp1315_rbi <- getInp(toa=toa1315_rbi_yaps, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType=ping_type, rbi_min=rbi_min, rbi_max=rbi_max, sdInits=1, ss_data_what="est", ss_data=0, biTable=NULL)
inp1315_rbi$inits[1] <- 1
yaps1315_rbi <- runYaps(inp1315_rbi, maxIter=1000, getPlsd=TRUE, getRep=TRUE, silent=TRUE)
```


```{r}
# basic plotting
plotYaps(inp1315_rbi, yaps1315_rbi, type="coord_X")
plotYaps(inp1315_rbi, yaps1315_rbi, type="coord_Y")
plotYaps(inp1315_rbi, yaps1315_rbi)

# uncertainties... standard error of position estimates
plot(yaps1315_rbi$plsd$X, type="l")
plot(yaps1315_rbi$plsd$Y, type="l")
plotNobs(toa1315_rbi_yaps)

# estimated speed of sound - sanity check
yaps1315_rbi$pl$ss

# look at the residual
matplot(t((yaps1315_rbi$rep$mu_toa - inp1315_rbi$datTmb$toa))*1450)


```



## Tag 1315 - pike - random BI - utilizing **known** sequence

```{r}
ping_type <- 'pbi'
seq1315 <- hald$burst_seqs$seq1315 # just a long sequence of random numbers...

# we need to figure out where our data is on this sequence
dat1315_kbi <- alignBurstSeq(synced_dat=dat1315, burst_seq=seq1315, seq_lng_min=10, rbi_min=rbi_min, rbi_max=rbi_max, plot_diag=TRUE)

```



```{r}
toa1315_kbi_list <- buildToaKnownSeq(seq=seq1315, aligned_dat=dat1315_kbi, hydros=hydros_yaps)
str(toa1315_kbi_list)
toa1315_kbi <- toa1315_kbi_list$toa
seq1315_kbi <- toa1315_kbi_list$seq
biTable1315 <- seq1315_kbi

# subset data...
ping_1 <- 1
n_ping <- 1000
toa1315_kbi_yaps <- toa1315_kbi[ping_1:(ping_1 + n_ping - 1), ]
biTable1315 <- seq1315_kbi[ping_1:(ping_1 + n_ping - 1) ]
plotToa(toa1315_kbi_yaps)
plotNobs(toa1315_kbi_yaps)

sum(getNobs(toa1315_kbi_yaps) >= 3) / nrow(toa1315_kbi_yaps)

```

```{r}
# get inp to feed into YAPS and run YAPS
inp1315_kbi <- getInp(toa=toa1315_kbi_yaps, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType=ping_type, rbi_min=rbi_min, rbi_max=rbi_max, sdInits=1, ss_data_what="est", ss_data=0, biTable=biTable1315)
inp1315_kbi$inits[1] <- 1
yaps1315_kbi <- runYaps(inp1315_kbi, maxIter=1000, getPlsd=TRUE, getRep=TRUE, silent=TRUE)

```


```{r}
# basic plotting
plotYaps(inp1315_kbi, yaps1315_kbi, type="coord_X")
plotYaps(inp1315_kbi, yaps1315_kbi, type="coord_Y")
plotYaps(inp1315_kbi, yaps1315_kbi)

# uncertainties... standard error of position estimates - compare to ping_type = 'rbi'
par(mfrow=c(2,2))
ylim <- c(0,50)
plot(yaps1315_kbi$plsd$X, type="h", ylim=ylim)
plot(yaps1315_rbi$plsd$X, type="h", ylim=ylim)
plot(yaps1315_kbi$plsd$Y, type="h", ylim=ylim)
plot(yaps1315_rbi$plsd$Y, type="h", ylim=ylim)
par(mfrow=c(1,1))

```



### The full track
The full track can be estimated using this code - will probably take some 5 - 10 minutes.
```{r eval=FALSE}
ping_type <- 'pbi'
rbi_min <- 60
rbi_max <- 120
seq1315 <- hald$burst_seqs$seq1315 # just a long sequence of random numbers...
dat1315_kbi <- alignBurstSeq(synced_dat=dat1315, burst_seq=seq1315, seq_lng_min=10, rbi_min=rbi_min, rbi_max=rbi_max, plot_diag=TRUE)
toa1315_kbi_list <- buildToaKnownSeq(seq=seq1315, aligned_dat=dat1315_kbi, hydros=hydros_yaps)
toa1315_kbi <- toa1315_kbi_list$toa
seq1315_kbi <- toa1315_kbi_list$seq
biTable1315 <- seq1315_kbi

ping_1 <- 1
n_ping <- nrow(toa1315_kbi)
toa1315_kbi_yaps <- toa1315_kbi[ping_1:(ping_1 + n_ping - 1), ]
biTable1315 <- seq1315_kbi[ping_1:(ping_1 + n_ping - 1) ]
plotToa(toa1315_kbi_yaps)
plotNobs(toa1315_kbi_yaps)

inp1315_kbi <- getInp(toa=toa1315_kbi_yaps, hydros=hydros_yaps, E_dist="t", n_ss=5, pingType=ping_type, rbi_min=rbi_min, rbi_max=rbi_max, sdInits=1, ss_data_what="est", ss_data=0, biTable=biTable1315)
inp1315_kbi$inits[1] <- 1
yaps1315_kbi <- runYaps(inp1315_kbi, maxIter=1000, getPlsd=TRUE, getRep=TRUE, silent=TRUE)
```


Lazy people can just grab it here:
```{r}
load('../data/yaps1315_kbi_fullTrack.rObj')
plotYaps(inp1315_kbi, yaps1315_kbi)

plotYaps(inp1315_kbi, yaps1315_kbi, type="coord_X")
plotYaps(inp1315_kbi, yaps1315_kbi, type="coord_Y")

```





# Appendix
## Obtaining the sync_model for the `hald` data
If you try to sync these data (which is encouraged), be aware of impossible sync tag - hydro combos.
Also note that positions of hydro idx 32 and 70 are uncertain, so they sould be estimated, i.e. removed from the fixed_hydros_idx vector. These settings should give a good sync model

```{r eval=FALSE}
	max_epo_diff <- 30
	min_hydros <- 5
	time_keeper_idx <- 1
	fixed_hydros_idx <- 1:70
	fixed_hydros_idx <- fixed_hydros_idx[!fixed_hydros_idx %in% c(32, 70)]
	n_offset_day <- 2
	n_ss_day <- 1
	inp_sync <- getInpSync(sync_dat, max_epo_diff, min_hydros, time_keeper_idx, fixed_hydros_idx, n_offset_day, n_ss_day, keep_rate=0.2)
	sync_model <- getSyncModel(inp_sync, silent=TRUE)
```









[yaps_logo]: yaps_logo_hex_100px.png
[otn_logo]: otn_logo.png

