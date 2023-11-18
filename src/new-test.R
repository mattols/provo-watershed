
#
#
#
# blank slate
#
library(terra);library(whitebox);library(sf);library(dplyr);library(ggplot2);library(raster)



### tutorial - https://vt-hydroinformatics.github.io/rgeowatersheds.html 
## Prepare DEM for hydrolgy analysis
wbt_breach_depressions_least_cost(
  dem = file.path(pth,"results/provo_huc8.tif"),
  output = file.path(pth,"dem/test2/huc8_provo_breached_d10.tif"),
  dist = 10, 
  fill = TRUE)

wbt_fill_depressions_wang_and_liu(
  dem = file.path(pth,"dem/test2/huc8_provo_breached_d10.tif"),
  output = file.path(pth,"dem/test2/huc8_provo_breached10_filled10.tif")
)

## create flow accumulation and pointers
wbt_d8_flow_accumulation(input = file.path(pth,"dem/test2/huc8_provo_breached10_filled10.tif"),
                         output = file.path(pth,"dem/test2/huc8_provo_D8FA.tif"))

wbt_d8_pointer(dem = file.path(pth,"dem/test2/huc8_provo_breached10_filled10.tif"),
               output = file.path(pth,"dem/test2/huc8_provo_D8pointer.tif"))



# # # # ## #

# snap pour points

wbt_extract_streams(flow_accum = file.path(pth,"dem/test2/huc8_provo_D8FA.tif"),
                    output = file.path(pth,"dem/test2/raster_streams2.tif"),
                    threshold = 1000)

wbt_jenson_snap_pour_points(pour_pts = "data/shp/pourpts/mora-pourpoints.shp",
                            streams = "data/dem/test2/raster_streams.tif",
                            output = "data/shp/pourpts/snappedpp.shp",
                            snap_dist = 0.0005) #careful with this! Know the units of your data

#### check data

pp <- shapefile("data/shp/pourpts/snappedpp.shp")
plot(pp)

library(raster)
d8 <- raster("data/dem/test2/huc8_provo_D8pointer.tif")
plot(d8)
plot(log(d8))

stream <- raster("data/dem/test2/raster_streams2.tif")
plot(stream,col='red')
plot(stream,col='red',ylim=c(40.26,40.33),xlim=c(-111.68,-111.62))

df8 = raster("data/dem/test2/huc8_provo_D8FA.tif")
plot(df8)
plot(df8,ylim=c(40.26,40.33),xlim=c(-111.68,-111.62))
plot(pp,add=T)
plot(st_transform(provo,crs(pp)),add=T)


# WTF - things are not lined up
slp = terrain(raster("data/results/provo_huc8.tif"))
asp = terrain(raster("data/results/provo_huc8.tif"),opt="aspect")
hsh = hillShade(slp, asp, 35, 290)
plot(hsh, col=grey.colors(30))
plot(hsh, col=grey.colors(30),ylim=c(40.26,40.33),xlim=c(-111.68,-111.62))

# not able to complete last step - no error and no file created





# # # # ## # DELINEATE WATERSHED

wbt_watershed(d8_pntr = "data/dem/test2/huc8_provo_D8pointer.tif",
              pour_pts = "data/shp/pourpts/snappedpp.shp",
              output = "datadem/test2/provo_mora_watersheds.tif")

ws <- raster("datadem/test2/provo_mora_watersheds.tif")
plot(ws)
