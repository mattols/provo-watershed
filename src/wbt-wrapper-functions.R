#
# Leveraging WTB for watershed delineation
# Jake Mora & Geoff Zahn
# goal: impervious surface % at each measurement point
#

library(terra);library(whitebox);library(sf);library(dplyr);library(ggplot2)
library(stars)

# wrapper to run WTB workflow
delineate_watershed <- function(dempath, outpath, breach_dist=5, fill_opt=T, 
                                stream_thresh=6000, snap_distance=0.5, flname_iter=NULL){
  # 
  # Calculate watershed based wtb toolbox 
  # tutorial - https://vt-hydroinformatics.github.io/rgeowatersheds.html 
  #
  
  # define
  shedoutpath <- dirname(dempath)
  #
  if (is.null(flname_iter)){
    outpath <- file.path(outpath,gsub(":","-",Sys.time()))
  }else{  outpath <- file.path(outpath,paste0("p_",flname_iter))}
  dir.create(outpath)
  
  # 1. breach and fill holes in DEM
  wbt_breach_depressions_least_cost(
    dem = dempath,
    output = file.path(outpath,"breached.tif"),
    dist = breach_dist, 
    fill = fill_opt)
  wbt_fill_depressions_wang_and_liu(
    dem = file.path(outpath,"breached.tif"),
    output = file.path(outpath,"breached_filled.tif")
  )
  
  # 2. create flow accumulation and pointer
  wbt_d8_flow_accumulation(input = file.path(outpath,"breached_filled.tif"),
                           output = file.path(outpath,"D8FA.tif"))
  wbt_d8_pointer(dem = file.path(outpath,"breached_filled.tif"),
                 output = file.path(outpath,"D8pointer.tif"))
  
  # 2.5 iterative single points (longitudinal) or multiple
  if(!is.null(flname_iter)){
    # alternative single point
    pts=st_read(file.path(shedoutpath,"wb_points.shp"))
    pourpnts = file.path(outpath, paste0("point",flname_iter,".shp") )
    st_write(pts[flname_iter,],pourpnts)
  }else{pourpnts=file.path(shedoutpath,"wb_points.shp")} # all points
  
  # 3. snap pour points to stream
  wbt_extract_streams(flow_accum = file.path(outpath,"D8FA.tif"),
                      output = file.path(outpath,"raster_streams.tif"),
                      threshold = stream_thresh)
  wbt_jenson_snap_pour_points(pour_pts =  pourpnts,
                              streams = file.path(outpath,"raster_streams.tif"),
                              output = file.path(outpath,"snappedpp.shp"),
                              snap_dist = snap_distance) #careful with this! Know the units of your data
  # 4. delineate watersheds
  wbt_watershed(d8_pntr = file.path(outpath,"D8pointer.tif"),
                pour_pts = file.path(outpath,"snappedpp.shp"),
                output = file.path(outpath,"watershed.tif"))
  
  print(cat("\n Completed delineation with...
      breach dist:", breach_dist,
      "
      fill opt:", fill_opt,
      "
      stream threshold:", stream_thresh,
      "
      snap distance:",snap_distance, "\n") )
  
  print(outpath)
}

# prep function
watershed_clip <- function(dempath, pointpath, shedpolypath, rivpath){
  #
  # clipping and projecting watershed area
  #
  
  # define
  shedoutpath <- dirname(dempath)
  # read
  shed = read_sf(shedpolypath)
  riv = read_sf(rivpath)
  dem = rast(dempath)
  # subset
  # full provo
  # pshed = shed %>% filter(HU_8_NAME=="Provo") %>% st_union() %>% st_sf()
  # lower provo
  pshed = shed %>% filter(HU_10_NAME=="Outlet Provo River") %>% st_union() %>% st_sf()
  
  # check projection of points
  if(crs(st_read(pointpath)) != crs(dem) ){
    st_read(pointpath) %>% st_transform(crs(rast(dempath))) %>%
      st_write(., file.path(shedoutpath,"wb_points.shp"))
  }else{ st_write(st_read(pointpath), file.path(shedoutpath,"wb_points.shp")) }
  
  # transform
  pshed = st_transform(pshed, crs(dem) )
  
  # write
  pdem = dem %>% crop(pshed)
  pmdem = pdem %>% mask(pshed)
  writeRaster(pdem, file.path(shedoutpath, "pshed.tif"))
  writeRaster(pmdem, file.path(shedoutpath, "pmask_shed.tif"))
  # poly write
  st_write(pshed,file.path(shedoutpath, "pshed.shp"))
  st_write(st_transform(riv, crs(shed)), file.path(shedoutpath, "priver.shp"))
  
  # print paths
  print(file.path(shedoutpath, "pshed.shp"))
  print(file.path(shedoutpath, "pshed.tif"))
  print(file.path(shedoutpath, "pmask_shed.tif"))
}

# basic plot to observe raster results
assess_results <- function(outpath, shedoutpath, clipit=TRUE){
  #
  # observe results of input raster layers
  #
  
  # read in watershed and river
  shed = read_sf(file.path(shedoutpath, "pshed.shp"))
  priver = read_sf(file.path(shedoutpath, "priver.shp"))
  
  # read files
  stk_list <- lapply( list(file.path(outpath,"watershed.tif"),
                   file.path(outpath,"D8pointer.tif"),
                   file.path(outpath,"D8FA.tif"),
                   file.path(outpath,"raster_streams.tif") ), function(x) rast(x) )
  stk <- rast( stk_list )
  names(stk) <- c("wb", "d8p","d8fa","strm")
  
  # watershed layer check
  plot(stk[[1]], main=paste("WB Rast - NaN values:",( sum(is.na(values(stk[[1]])))/length(values(stk[[1]])) )*100, "% of image" ))
  plot(shed, add=T, border='red',col=NA)
  plot(priver,add=T, col='deepskyblue')
  print(paste("NaN values:",( sum(is.na(values(stk[[1]])))/length(values(stk[[1]])) )*100, "% of image" ) )
  
  
  # clip
  # d <-  ext(xmin, xmax, ymin, ymax)
  d <- ext(443245.27310154, 449205.901673278, 4459263.94446199, 4468502.91874819)
  if(clipit){stk = crop(stk,d)}
  
  # plot all
  plot(stk)
  # #
}

# assess point snapping
assess_point_snap <- function(outpath, shedoutpath){
  #
  #  point snapping
  #
  
  # read
  priver = read_sf(file.path(shedoutpath, "priver.shp"))
  o_pnts = read_sf(file.path(shedoutpath, "wb_points.shp"))
  snp_pnts = read_sf(file.path(outpath, "snappedpp.shp"))
  d8fa = rast( file.path(outpath,"D8FA.tif") )
  # crop
  d8fa = crop(d8fa, snp_pnts)
  # plot
  plot(d8fa, main="Point alignment with streams")
  plot(o_pnts, add=T, col='blue', pch=3,cex=1)
  plot(snp_pnts, add=T,col='red', pch=4,cex=1)
  plot(priver, add=T, col='deepskyblue', cex=1)
  
}

watershed_hillshade <- function(outpath, shedoutpath, Zen = 40, Asp = 270, 
                                clipit=TRUE, multipts=NULL){
  #
  # graphics to display polygonized watersheds over hillshade
  #
  
  dem0 <- rast( file.path(shedoutpath,"pshed.tif") )
  # alt <- disagg(dem0, 10, method="bilinear") # smaller resolution?
  slope <- terrain(dem0, "slope", unit="radians")
  aspect <- terrain(dem0, "aspect", unit="radians")
  hill <- shade(slope, aspect, Zen, Asp)
  plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
  # plot(alt, col=rainbow(25, alpha=0.35), add=TRUE)
  
  # watershed 
  wb <- rast( file.path(outpath,"watershed.tif") )
  wsshape <- stars::st_as_stars(wb) %>% st_as_sf(merge = T)
  plot(wsshape,add=T,col=NA)
  
  # add single or multi points
  if(is.null(multipts)){
    plot(st_read(file.path(shedoutpath,"wb_points.shp") ),add=T, pch=4,col='red')
  }else{
    pnt_i = file.path(outpath, paste0("point",multipts,".shp") )
    plot(st_read(pnt_i ),add=T, pch=4,col='red')
  }
  
  # unable to plot text
  # wsshp = wsshape %>% st_centroid() %>% 
  #   # this is the crs from d, which has no EPSG code:
  #   st_transform(., crs(wsshape)) %>%
  #   # since you want the centroids in a second geometry col:
  #   st_geometry()
  # text(st_coordinates(st_centroid(wsshape))[1:2,], labels=(1:length(wsshape)))
  # 
  
  # clip
  # d <-  ext(xmin, xmax, ymin, ymax)
  d <- ext(443245.27310154, 449205.901673278, 4459263.94446199, 4468502.91874819)
  if(clipit){
    # single or multipoints
    if(is.null(multipts)){
      plot(crop(hill,d), col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
      plot(st_crop(wsshape,d),add=T,col=NA)
      plot(st_read(file.path(shedoutpath,"wb_points.shp") ),add=T, pch=4,col='red')
    }else{
      pnt_i = file.path(outpath, paste0("point",multipts,".shp") )
      pnt_buff = st_read(pnt_i) %>% st_buffer(4000) # 4km buffer 
      
      plot(crop(hill,pnt_buff), col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
      plot(st_crop(wsshape,pnt_buff),add=T,col=NA)
      plot(st_read(pnt_i),add=T, pch=4,col='red')
    }
  }
}

delineate_ws_iter <- function(dempath,outpath){
  #
  # runs watershed delineation for each point individually
  #
  strt = Sys.time()
  # define
  shedoutpath <- dirname(dempath)
  #
  pts_dim = dim(st_read(file.path(shedoutpath,"wb_points.shp") ))
  for (i in 1:pts_dim[1]){
    # each point
    print(paste("...beginning",i))
    delineate_watershed(dempath, outpath, breach_dist=5, fill_opt=T, stream_thresh=6000,
                        snap_distance=50, flname_iter = i)
    # plot results
    watershed_hillshade(file.path(outpath,paste0("p_",i)),
                        shedoutpath, Zen = 40, Asp = 270, 
                        clipit=TRUE, multipts = i)
    # must add i to any outpath plots
  }
  print("Finished")
  print(Sys.time() - strt)
}

# MAIN paths
dempath0 = "data/dem/n41w112_10m/n41w112_10m.tif"
outpath = "data/results"
pointpath = "data/shp/pourpts/mora-pourpoints.shp"

# watershed path
shedpolypath = "data/shp/Utah_Watersheds_Area/Watersheds_Area.shp"
rivpath = "data/shp/results/provo_riv_agg.shp"
shedoutpath = file.path(outpath, "dem_shed")
# # create watershed
# dir.create(shedoutpath)
# watershed_clip(dempath0, pointpath, shedpolypath, rivpath)
dempath = "data/results/dem_shed/pshed.tif"
dempath_mask = "data/results/dem_shed/pmask_shed.tif"

# check crs of DEM for inputs
crs(rast(dempath))
res(rast(dempath)) # in meters

# Test 1 
rpaths <- delineate_watershed(dempath, outpath, breach_dist=5, fill_opt=T, stream_thresh=1000,
                                snap_distance=50)
# Completed delineation with...
  # breach dist: 5 
  # fill opt: TRUE 
  # stream threshold: 1000 
  # snap distance: 50
# "data/results/2023-11-17 15-16-32"
assess_results(rpaths, shedoutpath, clipit=TRUE)
assess_point_snap(rpaths, shedoutpath)
#

# Test 2
rpaths <- delineate_watershed(dempath, outpath, breach_dist=5, fill_opt=T, stream_thresh=6000,
                              snap_distance=50)
 # rpaths <- "data/results/2023-11-17 15-21-05"
assess_results(rpaths, shedoutpath, clipit=TRUE)
assess_point_snap(rpaths, pointpath, shedoutpath)
watershed_hillshade(rpaths, shedoutpath, Zen = 40, Asp = 270, clipit=TRUE)
#
wb = rast( file.path(rpaths,"watershed.tif") )
plot(wb)
#

# # # # # # # # # #
# possibly need to run one point at a time. recreate individual watershed points
# new function?
dempath = "data/results/dem_shed/pshed.tif"
outpath = "data/results/ipoints"
delineate_ws_iter(dempath,outpath)
# observe
id = 10
watershed_hillshade(file.path(outpath,paste0("p_",id)), shedoutpath, Zen = 40, Asp = 270, 
                                clipit=TRUE, multipts=id)


# stream (6000 or 1000) does not include all upstream but only the nearby watershed
# perhaps need to add watersheds together


# eventually create a tool that combines all individual watersheds into a single layer



# extra
# pshedpath = "data/results/dem_shed/pshed.shp"

