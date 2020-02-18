install.packages("devtools")

# Tidyverse (data cleaning and arrangement)
install.packages('tidyverse')

# Mapping spatial data
install.packages('raster')
install.packages('mapdata')
install.packages('maptools')
install.packages('maps')
install.packages('ggplot2')
install.packages('ggmap')

# VTrack - Tools for Telemetry Analysis
devtools::install_github("rossdwyer/VTrack")

# GLATOS - acoustic telemetry package that does filtering, vis, array simulation, etc.
install.packages('remotes')
library(remotes)
install_url("https://gitlab.oceantrack.org/GreatLakes/glatos/repository/master/archive.zip",
            build_opts = c("--no-resave-data", "--no-manual"))  