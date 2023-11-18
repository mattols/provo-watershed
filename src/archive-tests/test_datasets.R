#
# Jake Mora
# incorporate "urbanization" in mycology communities study
#

# ideas
# buffer to impervious structures
# distance weighting
# compare water quality data

# Questions
# seasonal timing?
# location

# Utah GIS data
# https://gis.utah.gov/data/water/

# load river shapefile
library(terra);library(sf);library(dplyr);library(ggplot2)
path0 <- "C:\\Users\\10504912\\OneDrive - Utah Valley University\\Projects\\students\\JacobMora\\testdata"
list.files(path0)
riv <- read_sf(file.path(path0,"Utah_Streams_NHD/StreamsNHDHighRes.shp")) 

# parse for provo
head(riv)
provo = riv %>% filter(riv$GNIS_Name=="Provo River") %>% st_union %>% st_sf()
plot(provo)

# randomly select points
sf_linestring <- st_cast(provo, "LINESTRING")
# sampled_points <- st_line_sample(sf_linestring, n = 10)
sampled_points <- st_line_sample(sf_linestring, n = 1, type = 'regular')
sampled_points <- st_line_sample(sf_linestring, density = 0.0003, type = 'regular')

# plot multilinestring and point object
ggplot() +
  geom_sf(data=provo, color="blue") +
  geom_sf(data = sampled_points, color = "purple") +
  theme(panel.grid=element_line(color="transparent"))

# NLCD
# https://www.mrlc.gov/data/nlcd-2019-land-cover-conus
nlcd <- terra::rast(file.path(path0,"nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img"))
nlcd
pv_prj = st_transform(provo,st_crs(nlcd))
nlcd2 <- crop(nlcd, pv_prj) %>% project(.,crs(provo))
plot(nlcd2)
points(sampled_points, col='purple')

# 
df_tere_topo <- nlcd2 %>%
  as.data.frame(xy = TRUE) %>%
  rename(category = "NLCD Land Cover Class")
head(df_tere_topo)

ggplot() +
  geom_raster(data = df_tere_topo, aes(x = x, y = y, fill = category))+
  geom_sf(provo, mapping = aes(), color = 'black', fill = NA)

d <- draw()
# SpatExtent : -12418767.0598624, -12414883.0547005, 4911507.36594426, 4915468.28209956 (xmin, xmax, ymin, ymax)
d <- ext(-12418767.0598624, -12414883.0547005, 4911507.36594426, 4915468.28209956)

nlcd2 %>% crop(d) %>%  as.data.frame(xy = TRUE) %>%
  rename(category = "NLCD Land Cover Class") %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = category))+
  geom_sf(provo, mapping = aes(), color = 'black', fill = NA)


# save
st_write(provo,file.path(pth,"shp/results/provo_riv_agg.shp"))
st_write(sampled_points,file.path(pth,"shp/results/pv_sample_points.shp"))


# NLCD
# - impervous descriptor
# https://www.mrlc.gov/data?f%5B0%5D=category%3ALand%20Cover&f%5B1%5D=category%3AUrban%20Imperviousness&f%5B2%5D=year%3A2021
# https://www.mrlc.gov/nlcd-2021-science-research-products


# just a bit more
# inspect that watershed is correct
### SH from shed

sh <- st_crop(shed, st_transform(provo, st_crs(shed)) )

plot(nlcd2)
plot(st_transform(sh[,3], crs(nlcd2) ), col=NA, add=T)
coords = st_coordinates(st_cast(st_transform(sh[,3], crs(nlcd2) ),"POINT"))
text(coords, sh$HU_10_NAME)

nlcd2 %>% project(crs(pshed)) %>% crop(pshed) %>%
  plot()
plot(pshed, add=T, col=NA)
plot(opshed, add=T, border='blue')

# a bit small for upper provo - could add buffer to river clip for nlcd2