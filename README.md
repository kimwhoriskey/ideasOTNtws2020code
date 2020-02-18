<img src="ideasOTN.png" width="40%">
 
## Telemetry Workshop Series 2020
## Coding Workshops

The following code will be used in the Telemetry Workshop Series held at Dalhousie Univerisity from Jan 17-20, 2020. This series of lectures and workshops focuses on how to conduct a study with acoustic telemetry, how to effectively support and analyze the data, and appropriate methods for communicating your results and collaborating with agencies outside of academia. The coding lectures include the following topics: 

- Introduction to visualizing and analyzing telemetry data in R with dplyr and ggplot2 (author: Rob Lennox)
- Introduction to the glatos package (author: Ryan Gosse)
- Practically accounting for measurement error in acoustic telemetry (author: Jake Brownscombe)
- Predicting spatially continuous locations from discrete detection data using YAPS (author: Henrik Baktoft)
- Introduction to HMMs using moveHMM and markmodmover (author: Kim Whoriskey)

## Getting Started

Click the download button, or clone the repo using: 

```
git clone https://github.com/kimwhoriskey/ideasOTNtws2020code
``` 

## Prereqs
You will need R to work through these coding lectures. You may also find RStudio helpful (or some other environment). We will be using several packages - you can run the following code to install them. 

For Rob's lecture: 
```
install.packages("tidyverse")
install.packages("marmap")
install.packages("lubridate")
install.packages("gganimate")
install.packages("gifski")
install.packages("rgeos")
install.packages("argosfilter")
install.packages("nlme")
install.packages("igraph")
install.packages("ggraph")
install.packages("mgcv")
```

For Jake's lecture: 
```
install.packages("ggplot2")
install.packages("ggmap")
install.packages("gganimate")
install.packages("reshape2")
install.packages("dplyr")
install.packages("randomForest")
```

For Kim's lecture: 
```
install.packages("ggplot2")
install.packages("dplyr") 
install.packages("TMB")
install.packages("rgdal")
devtools::install_github("TheoMichelot/moveHMM", build_vignettes=TRUE)
devtools::install_github("lawlerem/markmodmover", build_vignettes=TRUE)
install.packages('circular')
install.packages('gridExtra')
```

Have fun!
