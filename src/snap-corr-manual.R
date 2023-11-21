#
# snap point correction (manual)
# one point only (P1)
#

# temporary fix for one point
# can create new function later to adjust pour point if necessary
# semi-automatic workflow

# # # # # #
# Correction PT1
# Manually snap point to main stream rather than smaller inlet
snp1 <- st_read(list.files(list.files(outpath,full.names = T)[1],full.names = T)[12] )
rs1 <-  rast(list.files(list.files(outpath,full.names = T)[1],full.names = T)[3] )
rs1 %>% crop(st_buffer(snp1, 100)) %>% plot()
plot(snp1,add=T,col='red',pch=4)
#
# extract rowcol numbers - find highest nearby flow (only look north south of pnt)
ext_info <- extract(rs1, snp1, cells=T)
rowcolnum <- rowColFromCell(rs1, ext_info$cell )
rs1[rowcolnum]==ext_info[,2] # check value is same
new_row <- which.max(rs1[,rowcolnum[,2]][,1])
new_row - rowcolnum[,1] # check distance (up two cells)
newsnp1 <- cellFromRowCol(rs1,new_row,rowcolnum[,2]) %>% xyFromCell(rs1,.) %>% as.data.frame() %>%
  st_as_sf(.,coords=c("x","y"), crs = crs(snp1))
# visual test
rs1 %>% crop(st_buffer(snp1, 100)) %>% plot()
plot(snp1,add=T,col='blue',pch=1)
plot(newsnp1,add=T,col='red',pch=4)
#

# RE-RUN
# rename filepath
file.rename(list.files(outpath,full.names = T)[1], file.path(outpath,paste0(basename(list.files(outpath,full.names = T)[1]),"_old") ) )
dir.create( file.path(outpath, "p_1" ) )
# not needed - creates new folder
# file.copy( file.path(outpath,paste0(basename(list.files(outpath,full.names = T)[1]),"_old") ),
#            file.path(outpath, "p_1"  ), recursive = T )
#
dempath = "data/results/dem_shed/pshed.tif" 
outpath = "data/results/ipoints/p_1"
# delineate_ws_iter(dempath,outpath) # not set 

#
# save new snapped point
st_write(newsnp1, file.path(outpath,"snappedpp.shp") )
st_write(snp1, file.path(outpath,"point1.shp") )

# 4. delineate watersheds
wbt_watershed(d8_pntr = file.path(outpath,"D8pointer.tif"),
              pour_pts = file.path(outpath,"snappedpp.shp"),
              output = file.path(outpath,"watershed.tif"))
# observe
watershed_hillshade(outpath, shedoutpath, Zen = 40, Asp = 270, 
                    clipit=TRUE, multipts=id)

