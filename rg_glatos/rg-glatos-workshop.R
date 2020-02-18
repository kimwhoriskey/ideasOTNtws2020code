library(glatos)

dets_file <- system.file("extdata", "blue_shark_detections.csv",
                         package = "glatos")

dets <- read_otn_detections(dets_file)
dets <- min_lag(dets)
dets <- false_detections(dets, tf = 3600)

library(tidyverse)

dets_f <- dets %>% filter(passed_filter != FALSE)

det_summary <- summarize_detections(dets, location_col = 'station')

events <- detection_events(dets_f, location_col = 'station')

event_summary <- events %>% 
  group_by(animal_id, location) %>% 
  summarise(
    number_of_events=n_distinct(first_detection),
    total_detections=sum(num_detections),
    total_time=sum(res_time_sec)
  )

rei <- residence_index(events)
rei <- residence_index(events, calculation_method = 'time_interval', time_interval_size = '1 hour')

abacus_plot(dets, location_col = 'station')

library(raster)

Canada <- getData('GADM', country="CAN", level=1)
NS <- Canada[Canada$NAME_1=="Nova Scotia",]

detection_bubble_plot(dets, location_col = 'station', map = NS,
                      background_xlim = c(-66, -62),
                      background_ylim = c(42, 46))

detection_bubble_plot(dets, location_col = 'station', map = NS,
                      background_xlim = c(-63.5, -63),
                      background_ylim = c(44, 44.5))


# This code will not work with the current version of glatos, the instructor will show it off

# dets_path <- system.file("extdata", "blue_shark_detections.csv",
#                          package = "glatos")
# deploy_path <- system.file("extdata", "hfx_deploy_simplified.xlsx",
#                            package = "glatos")
# tag_path <- system.file("extdata", "otn_nsbs_tag_metadata.xls",
#                         package = "glatos")

# dets <- read_otn_detections(dets_path)
# tags <- prepare_tag_sheet(tag_path, 5, 2)
# deploy <- prepare_deploy_sheet(deploy_path)

# ATTdata <- convert_otn_to_att(dets, tags, deploymentSheet = deploy)

# COAdata <- 
#   COA(ATTdata, 
#       timestep = 60) ## timestep used to estimate centers of activity (in minutes)

# View(COAdata)

