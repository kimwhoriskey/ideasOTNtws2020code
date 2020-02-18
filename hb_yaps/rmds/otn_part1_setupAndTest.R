# YAPS - Yet Another Positioning Solver
# OTN Telemetry Workshop Series, Dalhousie University, Halifax, Canada
# Part 1 - installation, test and a quick example
# Henrik Baktoft, DTU - hba@aqua.dtu.dk

## Installation
# If not already installed, please install `TMB` (Template Model Builder) and devtools. Then make sure you have the latest version of `yaps` installed. These steps are only needed once.

# install.packages('devtools')
# install.packages("TMB", type = "source")
# devtools::install_github('baktoft/yaps', ref='dev_ows')

# Load needed libaries - might throw messages about incorrect version of Matrix and/or registered S3 methods being overwritten - these should be safe to ignore fo now.
library(data.table)
library(dplyr)
library(sp)
library(leaflet)

# Test that `yaps` is working. This should show a short known track (in black) and estimated by `yaps` in red.
library(yaps)
testYaps(silent=TRUE)


## Getting started - the `ssu1` example data set

# This is a tiny data set collected as part of a feasibility study using YAPS on Vemco PPM style data to track fish in shallow parts of Florida Bay, USA. Data collected by J.S. Rehage, J.R. Rodemann, R.S. Corujo and N. Viadero. Included in `yaps` with permission from [J.S. Rehage](https://case.fiu.edu/about/directory/profiles/rehage-jennifer-schopf.html), FIU Florida International University.
# Have a look at the data - details can be found in `?ssu1`
names(ssu1)
head(ssu1$hydros)

# Pretty self explanatory. Coordinates are in UTM - YAPS will (most probably) not work well with lat/lon data. Column `sync_tag` indicate serial number of special transmitters co-located with the hydrophones; data from these are used in the synchronization process. Column `idx` is an index running from `1:nrow(hydros)`.

plot(y~x, data = ssu1$hydros, asp=1)
points(y~x, data=ssu1$hydros[!is.na(sync_tag)], col="red", pch=20)

head(ssu1$detections)
# Almost self explanatory. Each row is a detection of a transmitter (`tag`) on a hydrophone identified by `serial`. Column `ts` is the (non-synced) timestamp of the detection in timezone UTC. Column `epo` is `ts` converted to UNIX epoch using `as.numeric(ts)` and `frac` is fractions of second for the detection, i.e. the exact time of detection is given by `epofrac = epo + frac`.

head(ssu1$gps)

# Let's see where we are in the world (mainly included to have an excuse to use the [`leaflet`](https://cran.r-project.org/web/packages/leaflet/index.html) package - really awesome to make quick slippy-maps - code is in the .rmd-file)
hydros <- ssu1$hydros
coordinates(hydros) <- ~x+y
proj4string(hydros) <- CRS("+init=epsg:32617")
hydros_latlon <- spTransform(hydros, CRS("+init=epsg:4326"))
m <- leaflet(data=hydros_latlon, options = leafletOptions(minZoom = 0, maxZoom = 18), width="100%", height=700)
m <- addTiles(m, group="OSM")
m <- addCircles(m, radius=5, label=as.character(hydros_latlon$idx), labelOptions = labelOptions(noHide = T, textOnly = TRUE))
m <- addMeasure(m, primaryLengthUnit="meters")
m <- addProviderTiles(m, providers$Esri.WorldImagery, group="Esri.WorldImagery")
m <- addLayersControl(m, baseGroups = c("OSM (default)", "Esri.WorldImagery"),    options = layersControlOptions(collapsed = FALSE)  )
m


## Synchronizing the array
# The code below is identical to that presented in our pre-print [Opening the black box of high resolution fish tracking using yaps](https://www.researchgate.net/publication/338010182_Opening_the_black_box_of_high_resolution_fish_tracking_using_yaps), which also include detailed description of the parameters in `getInpSync()`.

# set sync parameters 
max_epo_diff <- 120
min_hydros <- 2
time_keeper_idx <- 5
fixed_hydros_idx <- c(2:3, 6, 8, 11, 13:17)
n_offset_day <- 2
n_ss_day <- 2

# get input data ready for getSyncModel()
inp_sync <- getInpSync(sync_dat=ssu1, max_epo_diff, min_hydros, time_keeper_idx, 
    fixed_hydros_idx, n_offset_day, n_ss_day, keep_rate=0.25)

# fit the sync model
sync_model <- getSyncModel(inp_sync, silent=TRUE)

# Plot model residuals and model check plots to ensure the synchronization process was successful...
plotSyncModelResids(sync_model, by='overall')    
plotSyncModelResids(sync_model, by='quantiles')
plotSyncModelResids(sync_model, by='sync_tag')      
plotSyncModelResids(sync_model, by='hydro')         

plotSyncModelCheck(sync_model, by="sync_bin_sync")  
plotSyncModelCheck(sync_model, by="sync_bin_hydro") 
plotSyncModelCheck(sync_model, by="sync_tag")       
plotSyncModelCheck(sync_model, by="hydro")          

# Apply the synchronization model to all data
detections_synced <- applySync(toa=ssu1$detections, hydros=ssu1$hydros, sync_model)

## Running `yaps` to estimate the track
# Prepare to estimate track using `yaps` on newly synchronized `ssu1` data
hydros_yaps <- data.table::data.table(sync_model$pl$TRUE_H)
colnames(hydros_yaps) <- c('hx','hy','hz')

# Specify focal tag and tag specific min and max burst intervals
focal_tag <- 15266
rbi_min <- 20
rbi_max <- 40

# Extract relevant data from the synced data
synced_dat_ssu1 <- detections_synced[tag == focal_tag]

# Compile TOA-matrix to use for yaps
toa_ssu1 <- getToaYaps(synced_dat_ssu1, hydros_yaps, rbi_min, rbi_max)

# Compile all input data needed for yaps
inp_ssu1 <- getInp(hydros_yaps, toa_ssu1, E_dist="Mixture", n_ss=2, pingType="rbi", 
    sdInits=1, rbi_min=rbi_min, rbi_max=rbi_max, ss_data_what="est", ss_data=0)

# Run yaps to obtain estimated track
yaps_out_ssu1 <- runYaps(inp_ssu1, silent=TRUE)

## Basic plotting of estimated track
# plot the estimated track
plotYaps(inp=inp_ssu1, yaps_out=yaps_out_ssu1, type="map")

# Add gps track for direct comparison
lines(utm_y~utm_x, data=ssu1$gps, lty=2)

par(mfrow=c(2,1))
plotYaps(inp=inp_ssu1, yaps_out=yaps_out_ssu1, type="coord_X")
lines(utm_x~ts, data=ssu1$gps, lty=2)
plotYaps(inp=inp_ssu1, yaps_out=yaps_out_ssu1, type="coord_Y")
lines(utm_y~ts, data=ssu1$gps, lty=2)

