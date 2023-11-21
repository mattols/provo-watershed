# 
# watershed-stats
#
library(tidyr)

# once delineation is complete and checked
# 1. create watershed polys - merge all into single layer
# 2. clip to NLCD/impervious - generate summary table


# generate watershed poly from individual point folders
watershed_poly_pts <-  function(outpath, shedoutpath, multipts=1){
  #
  # read in and save polygons
  
  # define outpath
  pointpath <-  file.path(outpath, paste0("p_",multipts) )
  
  # read data & transform
  wb <- rast( list.files(pointpath, full.names = T, pattern = "watershed.tif"))
  wsshape <- stars::st_as_stars(wb) %>% st_as_sf(merge = T)
  
  # lower provo watershed  (static)
  pvws <- st_read(list.files(shedoutpath, full.names = T, pattern = "pshed.shp"))
  
  #
  dem0 <- rast( file.path(shedoutpath,"pshed.tif") )
  plot(dem0, main=paste("Point", multipts))
  plot(pvws, add=T, col=NA, border='red')
  plot(wsshape,add=T,col=NA, border='blue')
  
  # clip
  interpv = st_intersection(pvws, wsshape)
  
  # save
  pt_id = rep(1:10,each=3)[multipts]
  savepath_w = file.path(dirname(outpath), "watersheds")
  st_write(wsshape,file.path(savepath_w, paste0("p",pt_id,"_full_shed.shp")) )
  st_write(interpv,file.path(savepath_w, paste0("p",pt_id,"_wshed.shp")) )
  
  # clip nlcd
  savepath_nlcd = file.path(dirname(outpath), "nlcd_clip")
  nlcd21 = rast("data/nlcd/nlcd/nlcd_2021_landcover_proj.tif") %>% crop(interpv) %>% mask(interpv)
  nlcdi = rast("data/nlcd/nlcd/nlcd_2021_impervious_proj.tif") %>% crop(interpv) %>% mask(interpv)
  nlcdi2 = rast("data/nlcd/nlcd/nlcd_2021_impervious_desc_proj.tif") %>% crop(interpv) %>% mask(interpv)
  # write
  writeRaster(nlcd21, file.path(savepath_nlcd, paste0("p",pt_id,"_nlcd.tif")) )
  writeRaster(nlcdi, file.path(savepath_nlcd, paste0("p",pt_id,"_imper.tif")) )
  writeRaster(nlcdi2, file.path(savepath_nlcd, paste0("p",pt_id,"_imdesc.tif")) )
  
}

# build final table with all ws
extract_lc_table <- function(df.table, outpath="data/results"){
  # 
  # extract variables for each site
  # ws_area, nlcd2021, imper2021, imperdesc
  
  df.table = df.table %>% mutate(lat = as.numeric(unlist(lapply(strsplit(lat,"_"), function(x) x[[1]][1])) ) ) %>% 
    mutate(lon = (as.numeric(unlist(lapply(strsplit(lon,"_"), function(x) x[[1]][1])) ) *-1) ) 
  
  # select first sample location & clean
  df.table1 = df.table[seq(1,30,3),]
  df.table1$site_ID = gsub("[A-Z]","",df.table1$site_ID)
  df.table1$collect_date <- gsub("\\.","",df.table1$collect_date)
  df.table1$process_date <- gsub("\\.","",df.table1$process_date)
  
  # iterate through points
  for (i in 1:nrow(df.table1)){
    print(paste("...extracting data from", i, "of", nrow(df.table1)))
    # extract info from layers
    interpv <- st_read(file.path(outpath,"watersheds", paste0("p",i,"_wshed.shp")) )
    v.area <- interpv %>% st_area() %>% as.numeric() / 1000^2 # in km^2
    if(i==1){df.table1$area_km2 = NA}
    df.table1$area_km2[i] = v.area
    
    # extract land cover
    nl21 <- rast( file.path(outpath,"nlcd_clip", paste0("p",i,"_nlcd.tif")) )
    # nl21i <- rast( file.path(outpath,"nlcd_clip", paste0("p",i,"_imper.tif")) )
    nl21i <- rast( file.path(outpath,"nlcd_clip", paste0("p",i,"_imdesc.tif")) )
    
    # summarize NLCD values
    nltbl = nl21 %>% as.data.frame() %>% na.omit() %>%
      group_by(`NLCD Land Cover Class`) %>% count() %>% mutate(n = (n*res(nl21)[1]^2)/1000^2 )
      # pivot_wider(names_from = `NLCD Land Cover Class`, values_from = n)
    # summarize Impervious values
    imtbl = nl21i %>% as.data.frame() %>% na.omit() %>%
      group_by(Class_Names) %>% count() %>% mutate(n = (n*res(nl21i)[1]^2)/1000^2 )
      # pivot_wider(names_from = Class_Names, values_from = n)
    
    # save table
    if (i==1){
      # df.nl = cbind(df.table1, nltbl)
      # df.i = cbind(df.table1, imtbl)
      df.nl = nltbl
      df.i = imtbl
      
    }else{
      # df.nl[i, 7:ncol(df.nl)] <- nltbl
      # df.i[i, 7:ncol(df.i)] <- imtbl
      # account for any new classes
      df.nl <- full_join(df.nl,nltbl,by="NLCD Land Cover Class")
      df.i <- full_join(df.i,imtbl,by= "Class_Names")
    }
  }
  # # rename and join tables
  colnames(df.nl)[2:11] <- paste0("p",1:10)
  colnames(df.nl)[2:11] <- paste0("p",1:10)
  
  # pivot
  # dfnl2 = df.nl %>% pivot_wider(names_from = `NLCD Land Cover Class`, values_from = n)
  dfnl2 = df.nl %>% t() %>% as.data.frame()
  colnames(dfnl2) = dfnl2[1,];dfnl2 <- dfnl2[2:nrow(dfnl2),];dfnl2[is.na(dfnl2)] <- 0.0
  # dfi2 = df.i %>% pivot_wider(names_from = Class_Names, values_from = n)
  dfi2 = df.i %>% t() %>% as.data.frame()
  colnames(dfi2) = dfi2[1,];dfi2 <- dfi2[2:nrow(dfi2),];dfi2[is.na(dfi2)] <- 0.0
  
  # bind
  dfinal.nlcd = cbind(df.table1, dfnl2)
  dfinal.imperv = cbind(df.table1, dfi2)

  # save tables
  write.csv(dfinal.nlcd,file.path(outpath, "tables", "NLCD_watersheds.csv"), row.names=FALSE)
  write.csv(dfinal.imperv,file.path(outpath, "tables", "Impervious_watersheds.csv"), row.names=FALSE)
  
  print("Complete")
  
}

# # # # # # # # # # 
# save files (watersheds and masked landcover tifs)
# only need 1 from each site! (3 samples per site)
outpath = "data/results/ipoints"
shedoutpath = "data/results/dem_shed"
for (i in seq(1,30,3)){ 
  print(paste("Running", i, "of", length(seq(1,30,3)) ))
  watershed_poly_pts(outpath, shedoutpath, i)
}
#

# test table variables
outpath = "data/results"
i = 1
nl21 <- rast( file.path(outpath,"nlcd_clip", paste0("p",i,"_nlcd.tif")) )
sum(!is.na(values(nl21)))
ncell(nl21)
summary(nl21)
tbl1 = nl21 %>% as.data.frame() %>% na.omit() %>% 
  mutate(nlcd = `NLCD Land Cover Class`) %>% select(nlcd) %>% group_by(nlcd) %>% count()
sum(tbl1$n)==sum(!is.na(values(nl21)))
tbl1 %>% t()

#
# # # # # # # # #
# extract tables
outpath = "data/results"
dfp = read.csv("data/sample_sites/aquatic-plant-proj-metadata.csv")
extract_lc_table(df.table = dfp, outpath="data/results")


# data frame test
dfp = read.csv("data/sample_sites/aquatic-plant-proj-metadata.csv")
head(dfp)
dfp %>% select(site_ID) %>% as.numeric()
gsub("[A-Z]","",dfp$site_ID)



