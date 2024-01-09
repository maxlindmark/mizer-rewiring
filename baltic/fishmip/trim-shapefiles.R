library(dplyr)
library(sf)
library(ggplot2)

home <- here::here()

baltic <- c('25', '26', '27', '28', '29', '32')

subdivs <- st_transform(
  st_read(paste0(home, "/baltic/fishmip/ICES_areas/ICES_Areas_20160601_cut_dense_3857.shp")), 4326) %>%
  filter(SubDivisio %in% baltic)

st_write(subdivs, paste0(home, '/baltic/fishmip/trim_shapefile.shp'))

# Test plot?

new_shapefile<- read_sf(paste0(home, '/baltic/fishmip/trim_shapefile.shp'))

ggplot(new_shapefile) +
  geom_sf(aes(fill = SubDivisio)) +
  theme_minimal()