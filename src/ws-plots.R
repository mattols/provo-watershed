#
# watershed plots
#



# Figure 1
figure01_hillsh <- function(outpath, shedoutpath){
  #
  #
  
  dem0 <- rast( file.path(shedoutpath,"pshed.tif") )
  # alt <- disagg(dem0, 10, method="bilinear") # smaller resolution?
  slope <- terrain(dem0, "slope", unit="radians")
  aspect <- terrain(dem0, "aspect", unit="radians")
  hill <- shade(slope, aspect, Zen, Asp)
  plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
  # plot(alt, col=rainbow(25, alpha=0.35), add=TRUE)
  
}