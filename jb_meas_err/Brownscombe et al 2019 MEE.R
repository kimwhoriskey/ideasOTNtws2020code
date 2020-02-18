# Brownscombe et al. (2020) A practical method to account for variation in detection range in acoustic telemetry 
# arrays to accurately quantify the spatial ecology of aquatic animals. Methods in Ecology and Evolution 11: 82-94

# This R script covers the methods used in the above paper to derive a detection range correction factor to correct 
# animal detection data for spatiotemporal variations in the detection ranges of a passive acoustic receiver array 

# This method involves collection of data from:
# 1. Range testing to derive the Maximum Range (MR) and Midpoint (50% detection efficiency)
# 2. Detection efficiency of a reference tag placed at the Midpoint throughout the study 
# 3. Animal detections

# In order to run the code below you will need to set your working directory to where you placed the data. You may
# also need to install the necessary packages before calling them. 

library(ggplot2)
library(ggmap)
library(gganimate)
library(reshape2)
library(dplyr)
library(randomForest)



# 1. Range testing data ####

detrange <- read.csv("detection_range_data.csv") #set your working directory to where these datafiles are on your computer
head(detrange)

# prior to this you need to take raw detection data and count the number of detections per hour to create this dataset. 
# Detection efficiency (DE) is the number of detections divided by the expected number of total detections in the sampling period
# based on the tag delay rate + the burst interval of the tag (with Vemco eqiupment). DE is expressed as a %, DEp as a proportion

#plot detection efficiency vs distance from the receiver by habitat type:
detrange$habitat <- factor(detrange$habitat, levels=c("bank","channel","basin"))
detrange$station <- factor(detrange$station, levels=c("BTT37", "BTT27", "BTT17","BTT49","BTT18","BTT54","BTT7","BTT13","BTT46"))
library(ggplot2)
ggplot(detrange, aes(x=distm, y=DE, col=habitat))+geom_point()+facet_wrap(~station)+
  coord_cartesian(ylim=c(0,100))+xlab("Distance (m)")+ylab("Detection efficiency (%)")


# Model this relationship for each site with a 3rd order polynomial function with the y-intercept set at 100 (% DE)
detrange$intercept <- 100 #set y intercept assuming 100% DE at 0 meters from the receiver
df_list <- split(detrange, f = paste(detrange$station) ) #break data into lists by station to model each seperately 

# apply models:
models <- lapply(df_list, function(dat) {
  lm(DE~-1+I(distm^3)+offset(intercept), data=dat)
})

# format result
DRres <- 
  do.call(
    rbind,
    lapply(models, function(x) {
      coefs <- summary(x)
      c( r2 = coefs$r.squared,
         adjr2 = coefs$adj.r.squared,
         AIC = AIC(x)
      )
    })
  )

DRres <- as.data.frame(DRres)
DRres$site <- as.factor(rownames(DRres))
DRres #assess model fit here. You might be better off with a different model for some or all sites


# predictions:
dist <- data.frame(distm=seq(from=0, to=1000, by=0.1))
dist$intercept <- 100

DRpreds <- 
  do.call(
    cbind,
    lapply(models, function(x) {
      round(predict.lm(x, dist),0)
    })
  )


DRpreds <- cbind(dist,DRpreds)
head(DRpreds)

# make long format for plotting:
library(reshape2)
DRpredslong <- melt(DRpreds, measure.vars=c(3:11))

library(dplyr)
DRpredslong <- DRpredslong %>% rename(station=variable, DE=value)
head(DRpredslong)

ggplot(detrange, aes(x=distm, y=DE, col=habitat))+
  stat_smooth(data=DRpredslong, aes(x=distm, y=DE),method="lm", se=TRUE, fill=NA,
             formula=y ~ 1+poly(x, 3, raw=TRUE),colour="black")+
  geom_point()+
  coord_cartesian(ylim=c(0,100), xlim=c(0,500))+
  facet_wrap(~station)+xlab("Distance from receiver (m)")+ylab("Detection efficiency (%)")


# pull out MR and Midpoint, where Midpoint=50% DE and MR=5% DE
head(DRpredslong)
DRpreds2 <- DRpredslong %>% filter(DE==50 | DE==5) %>% group_by(station, DE) %>% summarise(distm=round(mean(distm),0)) %>% as.data.frame()
DRpreds2$variable <- ifelse(DRpreds2$DE==5, "MR", "Midpoint")
head(DRpreds2)

Midpoint <- DRpreds2 %>% filter(variable=="Midpoint")
MR <- DRpreds2 %>% filter(variable=="MR")
Midpoint
MR


# add MR and Midpoint to plot:

ggplot(detrange, aes(x=distm, y=DE))+
  stat_smooth(data=DRpredslong, aes(x=distm, y=DE),method="lm", se=TRUE, fill=NA,
              formula=y ~ 1+poly(x, 3, raw=TRUE),colour="black")+
  geom_point(size=2)+
  geom_point(data=Midpoint, aes(x=distm, y=DE),col="green3", size=3, pch=17)+
  geom_point(data=MR, aes(x=distm, y=DE),col="red", size=3, pch=15)+
  coord_cartesian(ylim=c(0,100), xlim=c(0,500))+
  facet_wrap(~station)+xlab("Distance from receiver (m)")+ylab("Detection efficiency (%)")



# Midpoint can be used to guide the placement of a reference tag at each sentinel receiver site
# Maximum Range (intercept) can be combined with detection efficiency variance below to derive the detection range correction factor







# 2. Reference tag data #### 

refdata <- read.csv("reference_tag_data.csv")
str(refdata)
refdata$rugosity <- as.factor(refdata$rugosity) #make sure all variables are correct data type


# this is an hourly detection dataset. Similar to the range testing data above, you will need to summarise the number of 
# detections per study hour from your raw detections database. Here is some code to do that:

# detection data is presence only, so it won't have hours with zero detections in it. So start by generating an hourly dataset 
# for your study period:

refdata2 <- refdata
refdata2$datetime <- as.POSIXct(refdata2$hour, format="%Y-%m-%d %H:%M:%S") #change 'hour' to your datetime column. pay close attention to the format
refdata2$hour <- strftime(refdata2$datetime, format="%Y-%m-%d %H")
refdata2$hour <- as.POSIXct(refdata2$datetime,  format="%Y-%m-%d %H:%M:%S")
head(refdata2)

datetime <- data.frame(hour=seq(from=min(refdata2$hour), to=max(refdata2$hour), by="hour"))
str(datetime)

# unique stations:
stations <- data.frame(station=unique(refdata$station))
head(stations)

# merge the two, creating every hour and station combo:
datetimestation <- merge(datetime, stations)

# summarise the number of detections by study hour (you can also do this by other scales such as day)
refdata3 <- refdata2 %>% group_by(hour, station) %>% summarise(dets=length(station)) %>% as.data.frame()
head(refdata3)

# merge the hourly station data set with the detections:
refdata4 <- merge(datetimestation, refdata3, all.x=TRUE, by=c("hour","station"))

# and make the NAs zeros:
refdata4$dets[is.na(refdata4$dets)]<-0
head(refdata4)

# In this example the data are already summarised, so all detection counts return as one. In a raw dataset it should
# produce a range of detections. At this point you can assign your relevant environmental variables to the dataset (refdata4)
# as we have already with refdata:
head(refdata)

# calculate detection efficiency for each study hour
# our tag delay was 270 to 330 plus a 3 second burst period. So we should expect
round(3600/(333),0)
round(3600/(273),0)
# 11-13 detections per hour maximum
hist(refdata$dets)
# so we'll consider 11 or more detections an hour to be 100%

refdata$DE <- refdata$dets/11*100 
refdata$DE <-ifelse(refdata$DE>100,100, refdata$DE)
hist(refdata$DE) #frequency distribution of DE values
ggplot(refdata, aes(x=DE))+geom_histogram()+facet_wrap(~station) #by station

# calculate DEv - the difference between detection efficiency (DE) for each hour and the overall mean DE at each station:
refmeans <- refdata %>% group_by(station) %>% summarise(DEmean = mean(DE))
refdata$DEmean <- refmeans$DEmean[match(refdata$station, refmeans$station)]
refdata$DEv <- refdata$DE-refdata$DEmean


# If the mean DE at each receiver isn't 50 (which it will likely never be exactly), there's potential variance
# is unequal between receivers. This can be corrected using this equation, scaling the variance to +/- 50%
head(refdata)
refdata$DEvc <- ifelse(refdata$DE-refdata$DEmean<0, ((refdata$DE-refdata$DEmean)/refdata$DEmean)*50,
                                             ((refdata$DE-refdata$DEmean)/(100-refdata$DEmean))*50)

head(refdata)
hist(refdata$DEvc)




# Integrate MR with the data set to calculate DRc:
head(MR)
refdata$MR <- MR$distm[match(refdata$station, MR$station)]
str(refdata)


#calculate DRc: detrange + detrange*(DEvc/100)
refdata$DRc <- refdata$MR+refdata$MR*(refdata$DEvc/100)
hist(refdata$DRc)
min(refdata$DRc)


#examine the relationship between MR, DEvc and DRc:
ggplot(refdata, aes(x=MR, y=DRc, col=DEvc))+geom_point()+xlab("MR (m)")+
  ylab("DRc")+labs(colour="DEvc")+theme(legend.position=c(0.15,0.8))+
  theme(axis.line.x=element_line(),
        axis.line.y=element_line())



# DRc is a correction factor that integrates the estimated maximum range (MR) and time varying receiver 
# detection range (DEv). We have this data for each sentinel receiver, but need to predict it across the 
# entire receiver array. 



# Model DRc  ####
# based on environmental conditions. Consider including all factors you are interested in 
# exploring the effects of on animal positions (eg temperature, diel period) 

#Random forests model:

str(refdata)
refdata$rugosity <- as.factor(refdata$rugosity) #make sure variable types are consistent with training data


library(randomForest)
set.seed(1210)
z=formula(DRc~habitat+benthos+rugosity+depth+Diel+tidestate+tidemeters)
Refforest <- randomForest(formula=z, data = refdata, replace=FALSE, na.action=na.omit,
                          importance=TRUE, do.trace=1000, ntree=1000)
print(Refforest)
Refforest$importance
# model is performing alright, explaining 76% of variation in non training data. Most important predictors 
# based on %IncMSE are benthos type, depth, habitat, and rugosity 

# This model can be used to predict DRc on other stations, as done below












#animal detections ####

# DRc can be integrated into analyses of animal detection data in various ways as outlined
# in the associated manuscript. Below is a basic application of DRc to correct the number of animal
# detections per hour. 

permdets <- read.csv("permit_detections.csv")

str(permdets)
permdets$rugosity <- as.factor(permdets$rugosity) #variable types have to be consistent with training data

#predict DRc
permdets$DRc <- predict(Refforest, permdets)


#plot DRc spatially:

#summarise by site:
DRsite <- permdets %>% group_by(ID) %>% summarise(mDRc=mean(DRc), sd=sd(DRc), UCI=mDRc+sd*1.96,LCI=mDRc-sd*1.96,lat=mean(lat),lon=mean(lon))
head(DRsite)

#if you register an API with Google you can get a key to access their map content 
library(ggmap)
#register_google(key = "INSERT KEY HERE")
#FLmap <- get_googlemap(center = c(lon=-81.9, lat=24.6),
#                       zoom = 10, size = c(640, 640), scale = 4,
#                       maptype = c("satellite"))

FLmap <- readRDS('FLmap.rds')

ggmap(FLmap, extent='normal')+
  scale_x_continuous(limits=c(-82.05, -81.75))+
  scale_y_continuous(limits=c(24.45, 24.68))+
  ylab("Latitude") +
  xlab("Longitude")+
  geom_point(data=DRsite, aes(x=lon, y=lat, size=mDRc),fill="yellow", alpha=0.5, pch=21)



#animate spatiotemporally:
str(permdets)
permdets$hourpos <- as.POSIXct(permdets$hour, tz="US/Eastern", format="%Y-%m-%d %H")
permdets$lnDet <- log(permdets$Det+1)


#look at first couple weeks:
permdetssub <- permdets %>% filter(hourpos<"2018-06-08 00:00:00")


library(gganimate)

p <- ggmap(FLmap) + 
  theme(text = element_text(size = 17))+
  geom_point(data = permdetssub, mapping = aes(x = lon, y = lat, size=DRc, col=lnDet, group=ID),alpha=0.7)+
  scale_color_continuous(low = "yellow", high = "red")+
  scale_size_continuous(name="DRc", range=c(0.01,10))+
  labs(title = "{frame_time}")+
  scale_x_continuous(limits=c(-82.05, -81.75))+
  scale_y_continuous(limits=c(24.45, 24.68))+
  transition_time(hourpos)

length(unique(permdetssub$hourpos))
pAnim <- animate(p, duration=20, nframes=167, width= 600, height=600) 
#watch it:
pAnim

#save it:
anim_save("pAnim.gif")









#correct number of detections
permdets$Detc <- permdets$Det/(permdets$DRc/mean(permdets$DRc))
head(permdets)


#plot differences in presence only data (the equation above doesn't change zeros)
permdetsP <- permdets %>% filter(Det>0)

#make long format for plotting
permdetsmelt <- melt(permdetsP, measure.vars=c("Det","Detc"))
head(permdetsmelt)
permdetsmelt$logvalue <- log(permdetsmelt$value)

library(ggplot2)
permdetsmelt$variable <- factor(permdetsmelt$variable, levels=c("Det","Detc"))
permdetsmelt$habitat <- factor(permdetsmelt$habitat, levels=c("bank","channel","basin"))
ggplot(permdetsmelt, aes(x=Diel, y=logvalue,col=variable))+stat_summary(fun.y=mean, geom="point",position=position_dodge(0.2))+
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",width=0.1,position=position_dodge(0.2))+ylab("log(detections)")+
  facet_wrap(~habitat)+xlab("Diel period")+coord_cartesian(ylim=c(0,3))

ggplot(permdetsmelt, aes(x=tidemeters, y=logvalue,col=variable))+geom_smooth()+
  ylab("log(detections)")+
  facet_wrap(~habitat)+xlab("Tide height (m)")+coord_cartesian(ylim=c(0,3))







