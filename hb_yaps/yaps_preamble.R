# YAPS - Yet Another Positioning Solver
# ideasOTN workshop, Dalhousie University, Halifax, Canada - 2020-02-19

# Please make sure that these packages are installed.
install.packges(c('data.table','devtools', 'dplyr','sp','leaflet','lubridate','ggplot2','caTools','viridis')

# Please take a look at the TMB pages for info on installing TMB: https://github.com/kaskr/adcomp/wiki/Download
# Most often, this line works
install.packages("TMB", type = "source")

# Then install yaps from github, load it and check it is working:
# For this workshop, make sure to grab the 'dev_ows' branch (ref='dev_ows')!
devtools::install_github('baktoft/yaps', ref='dev_ows')
library(yaps)
testYaps()
# If the last line returned a plot of a simple track with overlapping black and red lines, everything should be working.

# You are also encouraged to take a look at the yaps readme: https://github.com/baktoft/yaps
