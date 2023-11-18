#
# Watershed delineation & pour point analysis
# Jake Mora - Mycology community project
#

# following: https://vt-hydroinformatics.github.io/rgeowatersheds.html
# WHITEBOX https://github.com/opengeos/whiteboxR

# Genrate D8 pointer grid
# fill in pits and breach depressions
# make pour points
# GPS points to flow network


# Our process will look like this:
#   Read in DEM
#   Fill single cell sinks then breach breach larger sinks
#   Create D8 flow accumulation and D8 pointer grids
#   Read in pour points
#   Create stream raster
#   Snap pour points to stream raster
#   *Run watershed function

# library(tidyverse)
# library(raster)
# library(sf)
# library(whitebox)
# library(tmap)
# library(stars)
# library(rayshader)
# library(rgl)



# DEM
library(terra);library(whitebox);library(sf);library(dplyr);library(ggplot2)
pth <- "C:/Users/10504912/OneDrive - Utah Valley University/Projects/students/JacobMora/data"
dem <- rast(file.path(pth,"dem/N40W112_SRTM/N40W112.hgt"))
plot(dem)

# watershed
shed = read_sf(file.path(pth,"shp/Utah_Watersheds_Area/Watersheds_Area.shp"))
head(shed)
plot(shed[,2])
shed %>% filter(HU_8_NAME=="Provo") %>% select(HU_10_NAME) %>% plot()
pshed = shed %>% filter(HU_8_NAME=="Provo") %>% st_union() %>% st_sf()
plot(pshed)
plot(st_transform(provo,st_crs(pshed)),add=T,col='blue')
opshed = shed %>% filter(HU_10_NAME=="Outlet Provo River") %>% st_union() %>% st_sf()
# st_write(opshed,file.path(pth,"shp/results/provo_huc8_watershed.shp"))

# checks in other code about size

pdem = dem %>% crop(opshed) %>% mask(opshed)
plot(pdem)
plot(st_transform(provo,st_crs(pshed)),add=T,col='blue')

## stamen map
library(ggmap)
pv_border <- c(bottom = 40.2, top = 40.56,left = -111.75, right = -111.43)
# bbox = c(left = -111.75, bottom = 40.2, right = -111.43, top = 40.56)
map_pv <- get_stamenmap(pv_border, zoom = 12, maptype = "terrain")
ggmap(map_pv) 

# 
elev_df <- pdem %>%
  as.data.frame(xy = TRUE) %>%
  rename(elevation = "N40W112")
head(elev_df)

ggmap(map_pv) +
  geom_raster(data = elev_df, aes(x = x, y = y, fill = elevation), alpha=0.5) +
  coord_cartesian() 
  # geom_sf(provo, mapping = aes(), color = 'blue', fill = NA)

writeRaster(pdem, file.path(pth,"results/provo_huc8.tif"))



###############################################################################################




#### Watershed modeling
## 1. FeaturePreservingSmoothing
wbt_feature_preserving_smoothing(
  dem = file.path(pth,"results/provo_huc8.tif"),
  output = file.path(pth,"results/smoothed_1clip.tif"),
  filter = 9
)
## 2. BreachDepressions
wbt_breach_depressions(dem = file.path(pth,"results/smoothed_1clip.tif"), 
                       output = file.path(pth,"results/breached_1clip.tif"))
## 3. DInfFlowAccumulation
wbt_d_inf_flow_accumulation(input = file.path(pth,"results/provo_huc8.tif"), output = file.path(pth,"results/flow_accum_1clip.tif"))

terra::plot(terra::rast(file.path(pth,"results/flow_accum_1clip.tif")))


# done





### tutorial - https://vt-hydroinformatics.github.io/rgeowatersheds.html 
## Prepare DEM for hydrolgy analysis
wbt_breach_depressions_least_cost(
  dem = file.path(pth,"results/provo_huc8.tif"),
  output = file.path(pth,"dem/N40W112_SRTM/huc8_provo_breached10.tif"),
  dist = 5, 
  fill = TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = file.path(pth,"dem/N40W112_SRTM/huc8_provo_breached10.tif"),
  output = file.path(pth,"dem/N40W112_SRTM/huc8_provo_breached10_filled.tif")
)

## create flow accumulation and pointers
wbt_d8_flow_accumulation(input = file.path(pth,"dem/N40W112_SRTM/huc8_provo_breached10_filled.tif"),
                         output = file.path(pth,"dem/N40W112_SRTM/huc8_provo_D8FA.tif"))

wbt_d8_pointer(dem = file.path(pth,"dem/N40W112_SRTM/huc8_provo_breached10_filled.tif"),
               output = file.path(pth,"dem/N40W112_SRTM/huc8_provo_D8pointer.tif"))

# ## setting pour points, must be snapped
# wbt_extract_streams(flow_accum = file.path(pth,"dem/N40W112_SRTM/huc8_provo_D8FA.tif"),
#                     output = file.path(pth,"dem/N40W112_SRTM/raster_streams.tif"),
#                     threshold = 6000)
# 
# # start here
# # points are:
# file.path(pth,"shp/results/provo_huc8_watershed.shp")

# points from Jacob

dfp = read.csv("data/sample_sites/aquatic-plant-proj-metadata.csv")
head(dfp)
dfp = dfp %>% mutate(lat = as.numeric(unlist(lapply(strsplit(lat,"_"), function(x) x[[1]][1])) ) ) %>% 
  mutate(lon = (as.numeric(unlist(lapply(strsplit(lon,"_"), function(x) x[[1]][1])) ) *-1) ) 

ppoints <- data.frame(
  Lon= dfp$lon, Lat = dfp$lat
)
head(ppoints)
plot(ppoints)
#

# ppointsSP <- SpatialPoints(ppoints, proj4string = CRS("+proj=longlat +datum=WGS84"))
# 
# ppointsSP <-st_as_sf(ppoints, coords=c("Lon","Lat"), proj4string = st_crs("+proj=longlat +datum=WGS84"))
# plot(ppointsSP)
# shapefile(ppointsSP, filename = "shp/results/mora-pourpoints.shp", overwrite = TRUE)
# raster::shapefile(spts, filename = "shp/results/mora-pourpoints.shp", overwrite = TRUE)

# use sp
library(sp)
spts <- SpatialPoints(ppoints, proj4string = CRS("EPSG:4326"))
is(spts, "SpatialPoints")
plot(st_transform(provo,crs(spts)) )
plot(spts,add=T, col='red')
spts

# st_write(st_as_sf(spts), "data/shp/pourpts/mora-pourpoints.shp", overwrite = TRUE)

# coords dont match
# spts2 <- spts
# crs(spts2) = crs(rast("data/dem/N40W112_SRTM/raster_streams.tif"))

# st_write(st_as_sf(spts2), "data/shp/pourpts/mora-pourpoints-2.shp", overwrite = TRUE)




# # # # ## #

# snap pour points

wbt_extract_streams(flow_accum = file.path(pth,"dem/N40W112_SRTM/huc8_provo_D8FA.tif"),
                    output = file.path(pth,"dem/N40W112_SRTM/raster_streams2.tif"),
                    threshold = 1000)

wbt_jenson_snap_pour_points(pour_pts = "data/shp/pourpts/mora-pourpoints.shp",
                            streams = "data/dem/N40W112_SRTM/raster_streams.tif",
                            output = "data/shp/pourpts/snappedpp.shp",
                            snap_dist = 0.0005) #careful with this! Know the units of your data

#### check data

pp <- shapefile("data/shp/pourpts/snappedpp.shp")
plot(pp)

library(raster)
d8 <- raster("data/dem/N40W112_SRTM/huc8_provo_D8pointer.tif")
plot(d8)
plot(log(d8))

stream <- raster("data/dem/N40W112_SRTM/raster_streams2.tif")
plot(stream,col='red')

df8 = raster("data/dem/N40W112_SRTM/huc8_provo_D8FA.tif")
plot(df8)
plot(df8,ylim=c(40.26,40.33),xlim=c(-111.68,-111.62))

# is working but there are gaps in the stream??


# revisit previous steps



# # # # ## # DELINEATE WATERSHED

wbt_watershed(d8_pntr = "data/dem/N40W112_SRTM/huc8_provo_D8pointer.tif",
              pour_pts = "data/shp/pourpts/snappedpp.shp",
              output = "datadem/N40W112_SRTM/provo_mora_watersheds.tif")

library(raster)
ws <- raster::raster("datadem/N40W112_SRTM/provo_mora_watersheds.tif")


tm_shape(hillshade)+
  tm_raster(style = "cont",palette = "-Greys", legend.show = FALSE)+
  tm_shape(ws)+
  tm_raster(legend.show = TRUE, alpha = 0.5, style = "cat")+
  tm_shape(pp)+
  tm_dots(col = "red")

#NEST
