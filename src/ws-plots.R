#
# watershed plots
#

# a few quick figures
library(reshape2);library(ggplot2)

# Figure 1
figure01_rast <- function(outpath, shedoutpath, pti = 1){
  #
  # watershed and land cover
  # create hillshade
  dem0 <- rast( file.path(shedoutpath,"pshed.tif") )
  slope <- terrain(dem0, "slope", unit="radians")
  aspect <- terrain(dem0, "aspect", unit="radians")
  hill <- shade(slope, aspect, 40, 270)
  # plot(alt, col=rainbow(25, alpha=0.35), add=TRUE)
  
  # add watershed
  interpv <- st_read(file.path(outpath,"watersheds", paste0("p",pti,"_wshed.shp")) )
  
  # add NLCD
  nl21 <- rast( file.path(outpath,"nlcd_clip", paste0("p",pti,"_nlcd.tif")) )
  nl21i <- rast( file.path(outpath,"nlcd_clip", paste0("p",pti,"_imdesc.tif")) )

  # final plot
  png(file.path(outpath,"figures",paste0("Watershed_Site",pti,"x.png")), width=11,height=5, units='in', res=300)
  par(mfrow=c(1,3))
  par(mar = c(2,0,4,0))
  par(oma=c(1,0,1,0))
  plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
  plot(interpv,add=T,col=NA, border='firebrick', lwd=1.5)
  par(mar = c(0,0,0,3))
  # par(oma=c(0,0,0,3))
  plot(nl21, axes=F, cex.lab=0.75)
  plot(nl21i, axes=F, cex.lab=0.75)
  # plot(interpv,add=T,col=NA, border='blue')
  dev.off()
  # reset
  par(mfrow=c(1,1))
  par(mar=c(5,3,3,2))
  par(oma=c(0,0,0,0))
  
  # try with layouts
  # https://bookdown.org/ndphillips/YaRrr/arranging-plots-with-parmfrow-and-layout.html
  # layout(matrix(c(1,2),ncol=2), heights = 1, widths=1)  
  # layout.show()
  
  #
}

# figure 2
figure02_nlcd <-  function(outpath){
  # use tables
  #
  df1 = read.csv(file.path(outpath, "tables", "NLCD_watersheds.csv") )
  df2 = read.csv(file.path(outpath, "tables", "Impervious_watersheds.csv") )
  
  # # # Land-cover # # #
  # sqkm coverage
  dfm1 <- df1 %>% select(-c(lat,lon,collect_date,process_date, area_km2)) %>% melt(.,id.vars="site_ID")
  dfm1 %>% ggplot(aes(x=site_ID,y=value,colour=variable)) + geom_point() +
    geom_line() + theme_classic()
  df1 %>% select(site_ID, Developed..Open.Space, Developed..Low.Intensity,Developed..Medium.Intensity, Developed..High.Intensity,
                  Hay.Pasture, Cultivated.Crops)  %>% melt(.,id.vars="site_ID") %>%
    ggplot(aes(x=site_ID,y=value,colour=variable)) + geom_point() +
    geom_line() + theme_classic() + ylab("Land cover type (km^2)") +
    scale_x_continuous(breaks=seq(1,10)) + xlab("Site ID")
  
  # percent coverage
  
  
  # # # IMPERVIOUS # # #
  dfm2 <- df2 %>% select(-c(lat,lon,collect_date,process_date, area_km2)) %>%
    mutate(site_ID = site_ID) %>% melt(.,id.vars="site_ID") %>% 
    filter(variable!="Unclassified")
  dfm2 %>% ggplot(aes(x=site_ID,y=value,colour=variable)) + geom_point() +
    geom_line() + theme_classic() + ylab("Impervious surface type (km^2)") +
    scale_x_continuous(breaks=seq(1,10)) + xlab("Site ID")
  
  # percent coverage
  # divide by area
  df2 %>% select(-c(lat,lon,collect_date,process_date, area_km2)) %>%
    group_by(site_ID) %>% rowSums()
  data.frame(site=1:10,impervious=allim) %>% 
    mutate(imperc = (impervious / df2$area_km2) * 100 ) %>%
    ggplot(aes(x=site,y=imperc)) + geom_line() + theme_classic()
  # ??? becomes larger but more impervious area
  
  
  # cumulative
  allim = df2 %>% select(-c(lat,lon,collect_date,process_date, area_km2, site_ID, Unclassified)) %>%
    rowSums()
  data.frame(site=1:10,impervious=allim) %>% 
    ggplot(aes(x=site,y=impervious)) + geom_line() + theme_classic()
  # ~325 km2 max ?
  
  # might be best to keep as m^2
  
}


# totals

allim = df2 %>% select(-c(lat,lon,collect_date,process_date, area_km2, site_ID, Unclassified)) %>%
  rowSums()
alllcdev = df1 %>% select(site_ID, Developed..Open.Space, Developed..Low.Intensity,Developed..Medium.Intensity, Developed..High.Intensity,
                          ) %>% rowSums()

# new plot
library(gridExtra)
p1 = df1 %>% select(site_ID, Developed..Open.Space, Developed..Low.Intensity,Developed..Medium.Intensity, Developed..High.Intensity) %>% 
  mutate(Total.Developed = alllcdev) %>%  melt(.,id.vars="site_ID") %>% mutate(Classification = variable) %>%
  ggplot(aes(x=site_ID,y=value,colour=Classification, linetype=Classification)) + geom_point() +
  geom_line() + theme_classic() + ylab("Land cover type (km^2)") +
  scale_x_continuous(breaks=seq(1,10)) + xlab("") + 
  scale_color_hue(l=50,c=80) + scale_linetype_manual(values=c(rep("solid",4),"dashed"))
  # scale_color_manual(values = )
p2 = df2 %>% select(-c(lat,lon,collect_date,process_date, area_km2)) %>%
  mutate(site_ID = site_ID) %>% mutate(total.impervious = allim) %>% melt(.,id.vars="site_ID") %>% 
  mutate(Classification = variable) %>% filter(variable!="Unclassified") %>% 
  ggplot(aes(x=site_ID,y=value,colour=Classification, linetype=Classification)) + geom_point() +
  geom_line() + theme_classic() + ylab("Impervious surface type (km^2)") +
  scale_x_continuous(breaks=seq(1,10)) + xlab("Site ID") + scale_color_hue(l=50,c=80) + 
  scale_linetype_manual(values=c(rep("solid",5),"dashed"))
p1
p2
# plot both
png(file.path(outpath,"figures",paste0("LandCover_Impervious_f",1,".png")), width=10,height=5, units='in', res=300)
grid.arrange(p1,p2,nrow=2)
dev.off()
#
# corrected for lat
p1 = df1 %>% select(lat, Developed..Open.Space, Developed..Low.Intensity,Developed..Medium.Intensity, Developed..High.Intensity) %>% 
  mutate(Total.Developed = alllcdev) %>%  melt(.,id.vars="lat") %>% mutate(Classification = variable) %>%
  ggplot(aes(x=lat,y=value,colour=Classification, linetype=Classification)) + geom_point() +
  geom_line() + theme_classic() + ylab("Land cover type (km^2)") +
  xlab("") + 
  scale_color_hue(l=50,c=80) + scale_linetype_manual(values=c(rep("solid",4),"dashed"))
# scale_color_manual(values = )
# scale_color_manual(values = )
p2 = df2 %>% select(-c(site_ID, lon,collect_date,process_date, area_km2)) %>%
  mutate(total.impervious = allim) %>% melt(.,id.vars="lat") %>% 
  mutate(Classification = variable) %>% filter(variable!="Unclassified") %>% 
  ggplot(aes(x=lat,y=value,colour=Classification, linetype=Classification)) + geom_point() +
  geom_line() + theme_classic() + ylab("Impervious surface type (km^2)") +
  xlab("Latitude") + scale_color_hue(l=50,c=80) + 
  scale_linetype_manual(values=c(rep("solid",5),"dashed"))
#
png(file.path(outpath,"figures",paste0("LandCover_Impervious_f",2,".png")), width=10,height=5, units='in', res=300)
grid.arrange(p1,p2,nrow=2)
dev.off()
#


# general map
pv_ws <- st_read(file.path(outpath,"dem_shed", "pshed.shp") )
pts0 <- st_read(file.path(outpath,"dem_shed", "wb_points.shp") )
pts1 <- pts0[seq(1,30,3),]
ptsdf <- pts1 %>% as.data.frame(); ptsdf <- ptsdf$geometry %>% st_coordinates() %>% as.data.frame();ptsdf$Site = (1:10)
# plot
png(file.path(outpath,"figures",paste0("Sample_Points_Map.png")), width=7,height=8.9, units='in', res=300)
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(pv_ws,add=T,col=NA, border='blue', lwd=2,lty=3)
text(ptsdf$X[c(1:4,8)], ptsdf$Y[c(1:4,8)], ptsdf$Site[c(1:4,8)], col='black', adj=c(1.45,0.25))
plot(pts0,add=T, col='firebrick',pch=4,cex=0.75)
text(ptsdf$X[c(5:7,10)], ptsdf$Y[c(5:7,10)], ptsdf$Site[c(5:7,10)], col='black', adj=c(0,1))
text(ptsdf$X[9], ptsdf$Y[9], ptsdf$Site[9], col='black', adj=c(-1.2,1.55))
text(441000,4460000, "Lower Provo River \nWatershed", cex=0.75, col='blue', srt=55)
text(457800,4476000, "Deer Creek Reservoir", cex=0.75, col='black', srt=55)
legend("topright",as.expression(bquote(bold("USGS 10-meter DEM"))),bty="n",cex=0.9, adj=c(0,2))
dev.off()
#



#
outpath = "data/results"
shedoutpath = "data/results/dem_shed"
figure01_rast(outpath, shedoutpath, pti = 8)
