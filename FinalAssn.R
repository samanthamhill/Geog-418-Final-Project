#DATA PREP
#Libraries
library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)
library(maptools)
library(sp)
library(rgeos)
library(gtable)
library(gridExtra)
library(grid)
library(ggplot2)
library("bcmaps")
library("bcmapsdata")

#Set working directory
dir <- "/Users/samanthahill/Desktop/Geog 418/Final Project"  
setwd(dir)

#Reading in elevation dataset
elev <- readOGR(dsn = "./VRI", layer = "ElevSample") #Read in data
elev <- spTransform(elev, CRS("+init=epsg:26910"))

#Reading in VRI data
VRI <- readOGR(dsn = "./VRI", layer = "WatershedVRI") #Read in shapefile
VRI <- spTransform(VRI, CRS("+init=epsg:26910"))
head(VRI@data)

vriCleanCols <- c("FID_VEG_CO", "POLYGON_ID", "PROJ_AGE_1",
                  "SITE_INDEX", "SPECIES__4", "SPECIES__5",
                  "PROJ_HEI_1", "SPECIES_PC", "SPECIES__6",
                  "VRI_LIVE_S", "BASAL_AREA", "WHOLE_STEM",
                  "CROWN_CL_1") 



#Create new file name "vriClean" and put 
#into it only the attributes in "vriCleanCols"
vriClean <- VRI[,vriCleanCols]

#Meta Data (https://www2.gov.bc.ca/assets/gov/farming-natural-resources-and-industry/forestry/stewardship/forest-analysis-inventory/data-management/standards/vegcomp_poly_rank1_data_dictionaryv5_2019.pdf)
# FID = Field ID
# PolyID = VRI Polygon ID
# Stand_Age = Estimated stand age projected to 2020 from estimated establishment date
# Site_Index = A value to estimate site quality. This describes the height that the stand could grow to by age 50 in meters.
# CoDom_Sp = The species code for the co-dominant tree species. Full list of codes: https://www.for.gov.bc.ca/hfp/publications/00026/fs708-14-appendix_d.htm
# Dom_Sp = The species code for the dominant tree species. Full list of codes: https://www.for.gov.bc.ca/hfp/publications/00026/fs708-14-appendix_d.htm
# Stand_HT = The estimated height for the stand
# DomSP_Perc = The estimated percentage of the dominent species
# CDomSP_Perc = The estimated percentage of the co-dominent species
# Stand_Dens = Estimated density of stand (Stems per hectare)
# Stand_BA = Estimated Basal area of the stand (square meters)
# Stand_StemBio = Estimated stand level stem biomass (tonnes per hectare)
# Stand_CrownCl = The percentage of ground area covered by tree crowns

newNames <- c("FID", "PolyID", "Stand_Age", "Site_Index",
              "CoDom_Sp", "Dom_Sp", "Stand_HT", "DomSP_Perc", 
              "CDomSP_Perc", "Stand_Dens", "Stand_BA", "Stand_StemBio", "Stand_CrownCl")

#set "vriClean" data column names to the "newNames" vector
colnames(vriClean@data) <- newNames


#Choose two variables, stand height and stand age
#Remove any null values from each of these attributes
vriClean <- vriClean[!is.na(vriClean@data$Stand_HT), ]



#Create choropleth map of height
map_Height <- tm_shape(vriClean) +
  tm_polygons(col = "Stand_HT",
              title = "Tree Stand Height (meters)",
              style = "jenks",
              palette = "viridis", n = 6,
              border.alpha = 0) +
  tm_legend(legend.position = c("LEFT", "BOTTOM"))

map_Height

#Study Area Map
bc <- as_Spatial(bc_neighbours()) #Get shp of BC bounds
bc <- spTransform(bc, CRS("+init=epsg:4326")) #project to WGS84 geographic (Lat/Long)

bc <- bc[which(bc$name == "British Columbia" ),] #Extract just the BC province


######################################################################################################################################################
#code retrieved from: https://www.jla-data.net/eng/adjusting-bounding-box-of-a-tmap-map/
bbox_new <- st_bbox(vriClean) # current bounding box

xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values

bbox_new[1] <- bbox_new[1] - (4.0 * xrange) # xmin - left
bbox_new[3] <- bbox_new[3] + (4.0 * xrange) # xmax - right
bbox_new[2] <- bbox_new[2] - (0.5 * yrange) # ymin - bottom
bbox_new[4] <- bbox_new[4] + (10.0 * yrange) # ymax - top

bbox_new <- bbox_new %>%  # take the bounding box ...
  st_as_sfc() # ... and make it a sf polygon
######################################################################################################################################################


studyArea <- tm_shape(bc, name = "BC", bbox = bbox_new)+
  tm_fill(col = "grey")+
  tm_shape(vriClean) +
  tm_fill(col = "red", borders.alpha = 0)+
  tm_legend(legend.position = c("LEFT", "BOTTOM"))

png("StudyArea.png")
studyArea
dev.off()

BCMap <- tm_shape(bc, name = "BC")+
  tm_fill(col = "grey") +
  tm_add_legend(type = c("text"), text = "Map created by Samantha Hill\nCoordinate Reference System: WGS84", size = 2, col = "black")+
  tm_compass(size = 3, position = c("RIGHT", "TOP"))+
  tm_scale_bar(size = 0.4)

png("BC_Area.png")
BCMap
dev.off()


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

#DESCRIPTIVE STATS
View(vriClean@data)

meanHeight <- mean(vriClean$Stand_HT, na.rm = TRUE)
meanHeight

modeHeight <- as.numeric(names(sort(table(vriClean$Stand_HT), decreasing = TRUE))[1])
modeHeight

medianHeight <- median(vriClean$Stand_HT, na.rm =TRUE)
medianHeight

sdHeight <- sd(vriClean$Stand_HT, na.rm = TRUE)
sdHeight

skewHeight <- skewness(vriClean$Stand_HT, na.rm = TRUE)
skewHeight

kurtHeight <- kurtosis(vriClean$Stand_HT, na.rm = TRUE)
kurtHeight

CoVHeight <- (sdHeight / meanHeight) * 100
CoVHeight

rangeHeight <- range(vriClean$Stand_HT)
rangeHeight

png("HeightHist.png")
hist(vriClean$Stand_HT, breaks = 30, main = "Histogram of Stand Height for the GVWSA", xlab = "Tree Height (meters)") #Base R style
dev.off()


meanElev <- mean(elev$grid_code, na.rm = TRUE)
meanElev

modeElev <- as.numeric(names(sort(table(elev$grid_code), decreasing = TRUE))[1])
modeElev

medianElev <- median(elev$grid_code)
medianElev

sdElev <- sd(elev$grid_code, na.rm = TRUE)
sdElev

skewElev <- skewness(elev$grid_code, na.rm = TRUE)
skewElev

kurtElev<- kurtosis(elev$grid_code, na.rm = TRUE)
kurtElev

CoVElev <- (sdElev / meanElev) * 100
CoVElev

rangeElev <- range(elev$grid_code)
rangeElev

png("ElevHist.png")
hist(elev$grid_code, breaks = 30, main = "Histogram of Elevations for the GVWSA", xlab = "Elevation (meters") #Base R style
dev.off()

#Make a table with all the descriptive stats
Samples = c("Tree Height", "Elevation") #Create an object for the labels
Means = c(round(meanHeight, 3), round(meanElev,3)) #Create an object for the means
S.D = c(round(sdHeight, 3), round(sdElev, 3)) #Create an object for the standard deviations
Median = c(round(medianHeight, 3), round(medianElev, 3)) #Create an object for the medians
Mode <- c(round(modeHeight, 3), round(modeElev, 3)) #Create an object for the modes
Skewness <- c(round(skewHeight, 3), round(skewElev, 3)) #Create an object for the skewness
Kurtosis <- c(round(kurtHeight, 3), round(kurtElev, 3)) #Create an object for the kurtosis
CoefficientOfVariation <- c(round(CoVHeight, 3), round(CoVElev, 3)) #Create an object for the CoV

Range <- c(round(rangeHeight, 3), round(rangeElev, 3))
temp1 <- c(Range[1], " - ", Range[2])
temp1
temp1 <- as.character(temp1)
temp1 <- toString(temp1)
temp1 <- sub(",", "", temp1)
temp1 <- sub(",", "", temp1)

temp2 <- c(Range[3], " - ",Range[4])
temp2
temp2 <- as.character(temp2)
temp2 <- toString(temp2)
temp2 <- sub(",", "", temp2)
temp2 <- sub(",", "", temp2)
Range<- c(temp1, temp2)
Range



data.for.table1 = data.frame(Samples, Means, Mode, Median, S.D, CoefficientOfVariation, Skewness, Kurtosis, Range)

#Make table 1
table1 <- tableGrob(data.for.table1, rows = c("","")) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 1: Descriptive Statistics for Tree Height and Elevation", gp = gpar(fontsize = 09))
padding <- unit(20, "mm")
table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1)

grid.arrange(table1, newpage = TRUE)

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

#Step 1: Determine if tree height is positively spatially autocorrelated

#1.2 Select neighborhoods

  #Default is Queens Weight
  vri.nb <- poly2nb(vriClean)

  #converts to line/network graph
  vri.net <- nb2lines(vri.nb, coords=coordinates(vriClean))

  #apply a coordinate system
  crs(vri.net) <- crs(vriClean)

  #weight matrix
  vri.lw <- nb2listw(vri.nb, zero.policy = TRUE, style = "W")

  print.listw(vri.lw, zero.policy = TRUE)

#1.2 Global Moran's I
  mi <- moran.test(vriClean$Stand_HT, vri.lw, zero.policy = TRUE)
  mi

  #1.2.1 Extract the results
  #Create varibles for moran's I, expected value, and variance
    mI <- mi$estimate[[1]]
    eI <- mi$estimate[[2]]
    var <- mi$estimate[[3]]
    
    #Z-test
    z <- ((mI-eI)/(sqrt(var)))
    z  
  
    mI <- formatC(mI, format = "e", digits = 2)
    eI <- formatC(eI, format = "e", digits = 2)
    var <- formatC(var, format = "e", digits = 2)
    z <- round(z,3)
    
    data.for.table3 = data.frame(mI, eI, var, z)
    table3 <- tableGrob(data.for.table3, rows = c("")) #make a table "Graphical Object" (GrOb) 
    t3Caption <- textGrob("Table 3: Results of Global Moran's I", gp = gpar(fontsize = 09))
    padding <- unit(20, "mm")
    table3 <- gtable_add_rows(table3, 
                              heights = grobHeight(t3Caption) + padding, 
                              pos = 0)
    
    table3 <- gtable_add_grob(table3,
                              t3Caption, t = 1, l = 2, r = ncol(data.for.table3) + 1)
    png("GlobMoranI_result_1.png")
    grid.arrange(table3, newpage = TRUE)
    dev.off()

    
    #1.2.2 Moran's I Range
    #Moran's I range
      moran.range <- function(lw) {
        wmat <- listw2mat(lw)
        return(range(eigen((wmat + t(wmat))/2)$values))
      }
      moran.range(vri.lw)
  
    
#1.3  Local Moran's I
  #1.3.1 Test
      lisa.test <- localmoran(vriClean$Stand_HT, vri.lw, zero.policy = TRUE)
      
      lisa.test
      
      #sticking the values in the polygons for mapping
      vriClean$Ii <- lisa.test[,1]
      vriClean$E.Ii<- lisa.test[,2]
      vriClean$Var.Ii<- lisa.test[,3]
      vriClean$Z.Ii<- lisa.test[,4]
      vriClean$P<- lisa.test[,5]
      
  #1.3.2 Map Local Moran's I
      map_LISA <- tm_shape(vriClean) + 
        tm_polygons(col = "Z.Ii", 
                    title = "Local Moran's I Z-score for Stand Height",
                    border.alpha = 0,
                    palette = "RdYlGn", n = 5, contrast = c(0,1),
                    style = "fixed",
                    breaks =  c(min(vriClean$Z.Ii, na.rm = TRUE),-1.96, 1.96, max(vriClean$Z.Ii, na.rm = TRUE)),
                    legend.show = TRUE,
                    labels = c("Negative Spatial Autocorrelation (Dispersed)", 
                               "Random",
                               "Positive Spatial Autocorrelation (Clustered)"),
                    legend.size = 15)+
        tm_compass(position = c("right","TOP"))+
        tm_scale_bar(position = "right")+
        tm_layout(title.size = 18, compass.type = "arrow", legend.width = 0.9, inner.margins = 0.1)
      
      png("LocalMoran_result_1.png")
      map_LISA
      dev.off()
  
  #1.3.3 Scatter Plot
      moran.plot(vriClean$Stand_HT, vri.lw, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Forest Stand Height", 
                 ylab="Forest Stand Height of Neighbours", quiet=NULL)

      
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
      
#Step 2: Create an interpolated surface for elevation
    #2.1.1 Create raster grid
      #Create a grid called grd to use in your interpolation
      # Create an empty grid where n is the total number of cells
      grd <- as.data.frame(spsample(elev, "regular", n=50000))
      names(grd)       <- c("X", "Y")
      coordinates(grd) <- c("X", "Y")
      gridded(grd)     <- TRUE  # Create SpatialPixel object
      fullgrid(grd)    <- TRUE  # Create SpatialGrid object
      proj4string(grd) <- proj4string(elev)

    #2.1.2 Omit null values
      #View(elev@data)
      #plot(elev)
      elev = na.omit(elev)
      
      
    #2.1.3 IDW
      P.idw <- gstat::idw(grid_code ~ 1, elev, newdata=grd, idp=3)
      r       <- raster(P.idw)
      r.m     <- mask(r, vriClean) #clips to area
      
      #2.1.3.1 Map the IDW surface
        IDW_map <- tm_shape(r.m) + 
        tm_raster(n=10,palette = "viridis",
                  title="Predicted Elevation\nin meters") + 
        tm_shape(elev) + tm_dots(size=0.1) +
        tm_layout(inner.margins = 0.1, title.size = 0.7, title.position = c("left", "TOP"))+
        tm_legend(legend.outside=TRUE, legend.outside.position = c("right"), position = c("LEFT","TOP"))+
        tm_compass(position = c("RIGHT","BOTTOM"), size = 2)+
        tm_scale_bar(position = c("RIGHT","BOTTOM"))
        
        png("IDW_surface.png")
        IDW_map    
        dev.off()
    
    #2.1.4 Cross Validation Leave One Out
        IDW.out <- vector(length = length(elev))
        for (i in 1:length(elev)) {
          IDW.out[i] <- idw(grid_code ~ 1, elev[-i,], elev[i,], idp=4)$var1.pred
        }
        
        # Plot the differences
        png("LOOCV_plot.png")
        OP <- par(pty="s", mar=c(4,3,0,0))
        plot(IDW.out ~ elev$grid_code, asp=1, xlab="Observed", ylab="Predicted", pch=16,
             col=rgb(0,0,0,0.5))
        abline(lm(IDW.out ~ elev$grid_code), col="red", lw=2,lty=2)
        abline(0,1)
        par(OP)
        dev.off()
        
        # if good estimator you would expect the red line to be the same at the black line.
        
        
        #root mean square error
        RMSE_3 <- sqrt(sum((IDW.out - elev$grid_code)^2) / length(elev))
        RMSE_2 <- sqrt(sum((IDW.out - elev$grid_code)^2) / length(elev))
        RMSE_4 <- sqrt(sum((IDW.out - elev$grid_code)^2) / length(elev))
        p_equal_2 <- round(RMSE_2, 3)
        p_equal_3 <- round(RMSE_3,3)
        p_equal_4 <- round(RMSE_4, 3)
        
        data.for.table4 = data.frame(p_equal_2, p_equal_3, p_equal_4)
        table4 <- tableGrob(data.for.table4, rows = c("")) #make a table "Graphical Object" (GrOb) 
        t4Caption <- textGrob("Table 2: RMSE for IDW Interpolation", gp = gpar(fontsize = 09))
        padding <- unit(20, "mm")
        table4 <- gtable_add_rows(table4, 
                                  heights = grobHeight(t4Caption) + padding, 
                                  pos = 0)
        
        table4 <- gtable_add_grob(table4,
                                  t4Caption, t = 1, l = 2, r = ncol(data.for.table4) + 1)
        png("RMSE_table.png")
        grid.arrange(table4, newpage = TRUE)
        dev.off()
        
      
      #2.1.5 Cross Validation Jack Knifing
        #we calculate how confident you are in each prediction
        img <- gstat::idw(grid_code ~ 1, elev, newdata=grd, idp=3) # check IDP value
        n   <- length(elev)
        Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)
        
        # Remove a point then interpolate (do this n times for each point)
        #same as above except stack the layers up.
        st <- stack()
        for (i in 1:n){
          Z1 <- gstat::idw(grid_code ~ 1, elev[-i,], newdata=grd, idp=3)
          st <- addLayer(st,raster(Z1,layer=1))
          # Calculated pseudo-value Z at j
          Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
        }
        
        # Get a confidence interval;
        # Jackknife estimator of parameter Z at location j
        Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )
        
        # Compute (Zi* - Zj)^2
        c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
        c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference
        
        # Compute the confidence interval
        CI <- sqrt( 1/(n*(n-1)) * c1)
        
        # Create (CI / interpolated value) raster
        img.sig   <- img
        img.sig$v <- CI /img$var1.pred 
        
        #Clip to area
        # Clip the confidence raster to GVWSA
        r.conf <- raster(img.sig, layer="v")
        r.m.conf <- mask(r.conf, vriClean)
        
        Jack_map <- tm_shape(r.m.conf) + 
          tm_raster(n=10,palette = "viridis",
                    title="Confidence Predicted Elevation\nin meters") + 
          tm_shape(elev) + tm_dots(size=0.1) +
          tm_layout(inner.margins = 0.1, title.size = 0.7, title.position = c("left", "TOP"))+
          tm_legend(legend.outside=TRUE, legend.outside.position = c("right"), position = c("LEFT","TOP"))+
          tm_compass(position = c("RIGHT","BOTTOM"), size = 2)+
          tm_scale_bar(position = c("RIGHT","BOTTOM"))
        
        png("jackknife_map.png")
        Jack_map  
        dev.off()
        
    #2.2 Combining and extracting the mean
        r <- raster(P.idw)
        sufaceMap <- tm_shape(r) + 
          tm_raster(n=5,palette = "viridis",
                    title="Elev (m)") +
          tm_shape(elev) + tm_dots(size=0.2)
        sufaceMap
        
        #If you have too many cells, 
        #you can reduce the number by aggregating values
        #agg <- aggregate(yourRasterFromKriging, fact=??, fun=mean)
        
        #Extract average elev for each polygon
        vriClean$Elev <- extract(r, vriClean, fun = mean)[,1]
        
        
        map_MeanElev <- tm_shape(vriClean) +
          tm_polygons(col = "Elev",
                      title = "Elevation (meters)",
                      style = "jenks",
                      palette = "viridis", n = 9,
                      border.alpha = 0) +
          tm_legend(legend.position = c("LEFT", "BOTTOM"))
        
        map_MeanElev
      
        View(vriClean@data)   

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
        
#Step 3 Regression Analysis
  #Let's say your dataset with both Elev and Height are stored in a dataset called VRI.
  #Plot Height and Elev from the VRI dataset you created
        
  # independent = elevation ; dependent = stand height
  plot(vriClean$Stand_HT ~ vriClean$Elev)
        
  #Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
  VRI.no0 <-  vriClean[which(vriClean$Stand_HT > 0), ]
  VRI.no0 <-  VRI.no0[which(VRI.no0$Elev > 0), ]
        
  #Now plot the data again
  plot(VRI.no0$Stand_HT ~ VRI.no0$Elev)
        
  #Perform a linear regression on the two variables. You should decide which one is dependent.
  lm.model <- lm(VRI.no0$Stand_HT ~ VRI.no0$Elev)
        
  #Add the regression model to the plot you created
  Tree_Height <- VRI.no0$Stand_HT
  Elevation <- VRI.no0$Elev
  
  png("Linear_Regression.png")
  plot(Tree_Height ~ Elevation)
  abline(lm.model, col = "red")
  dev.off()
        
  #Get the summary of the results
  summary(lm.model)
        
  #add the fitted values to your spatialpolygon dataframe
  VRI.no0$predictlm <- lm.model$fitted.values
        
  #You want to determine if the model residuals are spatially clustered. 
  #add the residuals to your spatialpolygon dataframe
  VRI.no0$residuals <- residuals.lm(lm.model)
        
  #Observe the result to make sure it looks correct
  head(VRI.no0@data)
  plot(VRI.no0$residuals)
        
  #Now, create choropleth map of residuals
  map_resid <- tm_shape(VRI.no0) +
  tm_polygons(col = "residuals",
    title = "Stand Height Residuals",
    style = "jenks",
    palette = "RdYlGn", n = 6, border.alpha = 0)+
    tm_compass(position = c("right", "top"))
      
  png("residualMap.png")
  map_resid
  dev.off()

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

#Step 4  Determine if we need to do GWR
  #2.4.1 Weight Matrix
  #Default is Queens Weight
  vri.nb <- poly2nb(VRI.no0)
  
  #converts to line/network graph
  vri.net <- nb2lines(vri.nb, coords=coordinates(VRI.no0))
  
  #apply a coordinate system
  crs(vri.net) <- crs(VRI.no0)
  
  #weight matrix
  vri.lw <- nb2listw(vri.nb, zero.policy = TRUE, style = "W")
  
  print.listw(vri.lw, zero.policy = TRUE)
  
  #4.2 Global Moran's I
  mi <- moran.test(VRI.no0$residuals, vri.lw, zero.policy = TRUE)
  mi
  
  #4.2.1 Extract the results
  #Create varibles for moran's I, expected value, and variance
  mI <- mi$estimate[[1]]
  eI <- mi$estimate[[2]]
  var <- mi$estimate[[3]]
  
  #Z-test
  z <- ((mI-eI)/(sqrt(var)))
  z  
  
  #4.2.2 Moran's I Range
  #Moran's I range
  moran.range <- function(lw) {
    wmat <- listw2mat(lw)
    return(range(eigen((wmat + t(wmat))/2)$values))
  }
  moran.range(vri.lw)

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
  
#Step 5 GWR
  #Let's say you are continuing with 
  #your data from the regression analysis. 
  #The first thing you need to do is to add the 
  #polygon coordinates to the spatialpolygondataframe.
  #You can obtain the coordinates using the 
  #"coordinates" function from the sp library
  VRI.no0.coords <- sp::coordinates(VRI.no0)
  #Observe the result:
  head(VRI.no0.coords)
  #Now add the coordinates back to the spatialpolygondataframe
  VRI.no0$X <- VRI.no0.coords[,1]
  VRI.no0$Y <- VRI.no0.coords[,2]
  
  ###Determine the bandwidth for GWR: this will take a while
  GWRbandwidth <- gwr.sel(VRI.no0$Stand_HT ~ VRI.no0$Elev, 
                          data=VRI.no0, coords=cbind(VRI.no0$X,VRI.no0$Y),adapt=T) 
  
  ###Perform GWR on the two variables with the bandwidth determined above
  ###This will take a looooooong while
  gwr.model = gwr(VRI.no0$Stand_HT ~ VRI.no0$Elev, 
                  data=VRI.no0, coords=cbind(VRI.no0$X,VRI.no0$Y), 
                  adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 
  
  #Print the results of the model
  gwr.model
  
  #Look at the results in detail
  results<-as.data.frame(gwr.model$SDF)
  head(results)
  
  #Now for the magic. Let's add our local r-square values to the map
  VRI.no0$localr <- results$localR2
  
  #VRI.no0.1 <-  VRI.no0[which(VRI.no0$localr >= -1.043), ]
  
  #Create choropleth map of r-square values
  map_r2 <- tm_shape(VRI.no0) +
    tm_polygons(col = "localr",
                title = "R2 values",
                style = "jenks",
                palette = "viridis", n = 6, border.alpha = 0)+
    tm_compass(position = c("right","top"))+
    tm_layout(title = "Map of Local R-squared Values", inner.margins = 0.1)+
    tm_scale_bar()
  
  png("R_squared.png")
  map_r2
  dev.off()
  
  #hist(VRI.no0$localr)
  
  #Time for more magic. Let's map the coefficients
  VRI.no0$coeff <- results$VRI.no0.Elev
  #Create choropleth map of the coefficients
  map_coef <- tm_shape(VRI.no0) +
    tm_polygons(col = "coeff",
                title = "Coefficients",
                style = "fisher",
                palette = "RdYlGn", n = 5, border.alpha = 0)+
    tm_compass(position = c("right","top"))+
    tm_layout(title = "Map of Coefficient Values", inner.margins = 0.1)+
    tm_scale_bar()
  
  png("Ceoff.png")
  map_coef
  dev.off()

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
  #Step 6 PPA to see if its clustered or not
  #6.1 Create point pattern analysis object
      kma <- elev
      kma$x <- coordinates(kma)[,1]
      kma$y <- coordinates(kma)[,2]
  
    #check for and remove duplicated points
    #first, finds zero distance among points to see if there are any duplicates
      zd <- zerodist(kma)
      zd
  
    #if there are duplicates, remove them
      kma <- remove.duplicates(kma)
  
    #create an "extent" object which can be used to create the observation window for spatstat
      kma.ext <- as.matrix(extent(vriClean))
      
  
    #observation window
      window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,]))
  
  
    #create ppp oject from spatstat
      kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)
  
  
  #6.2 Nearest Neighbour Distance
      nearestNeighbour <- nndist(kma.ppp)
  
    #Convert the nearestNeighbor object into a dataframe. BIG LIST OF DISTANCES
      nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
    #Change the column name to "Distance"
      colnames(nearestNeighbour) = "Distance"
  
    #determine number of elevation points
      View(elev@data)
      n = 365
    #Calculate the nearest neighbor statistic to test for a random spatial distribution.
    #mean nearest neighbour
      nnd = ((sum(nearestNeighbour$Distance))/(n))
      nnd  
  
    #study area calculation
      area <- gArea(vriClean, byid = FALSE)
      area
  
      studyArea <-  area  
      pointDensity <- n/area
  
      studyArea
      pointDensity
  
      #random NND
        r.nnd = 1/((2)*(sqrt(pointDensity)))
  
      #dispersed NND
        d.nnd = ((1.07453)/(sqrt(pointDensity)))
  
      #coefficient of variation
        R = nnd/r.nnd
  
      #sigma NND
        SE.NND <- ((0.26136)/(sqrt((n*pointDensity))))
  
      #calculate z score
        z = ((nnd - r.nnd)/SE.NND)
  
  #MAKE A TABLE
  NND <- c(round(nnd, 2))
  NNDr <- c(round(r.nnd, 2))
  NNDd <- c(round(d.nnd, 2))
  Sigma <- c(round(SE.NND,2))
  Z_score <- c(round(z,2))
  
  data.for.table5 = data.frame(NND, NNDr, NNDd, Sigma, Z_score)
  table5 <- tableGrob(data.for.table5, rows = c("")) #make a table "Graphical Object" (GrOb) 
  t5Caption <- textGrob("Table 3: Results of NND PPA", gp = gpar(fontsize = 09))
  padding <- unit(20, "mm")
  table5 <- gtable_add_rows(table5, 
                            heights = grobHeight(t5Caption) + padding, 
                            pos = 0)
  
  table5 <- gtable_add_grob(table5,
                            t5Caption, t = 1, l = 2, r = ncol(data.for.table5) + 1)
  
  png("NND_results.png")
  grid.arrange(table5, newpage = TRUE)
  dev.off()
  
  #6.3 KDE 
  kde.50 <- density(kma.ppp, sigma = 50, at = "pixels", eps = c(500, 500))
  kde.SG <- as(kde.50, "SpatialGridDataFrame")
  kde.100 <- density(kma.ppp, sigma = 100, at = "pixels", eps = c(500, 500))
  kde.SG <- cbind(kde.SG, as(kde.100, "SpatialGridDataFrame"))
  kde.250 <- density(kma.ppp, sigma = 250, at = "pixels", eps = c(500, 500))
  kde.SG <- cbind(kde.SG, as(kde.250, "SpatialGridDataFrame"))
  kde.500 <- density(kma.ppp, sigma = 500, at = "pixels", eps = c(500, 500))
  kde.SG <- cbind(kde.SG, as(kde.500, "SpatialGridDataFrame"))

  names(kde.SG) <- c("Size50", "Size100", "Size250", "Size500")  
  spplot(kde.SG, main = "KDE Result for 4 Bandwidths for Elevation")  

  #optimized bandwidth
  bw.d <- bw.diggle(kma.ppp)        
  plot(bw.d, ylim=c(-10, 10), main= "Cross-Validation Function for KDE of Elevation") 
  
  kde.bwo <- density(kma.ppp, sigma = bw.d, at = "pixels", eps = c(100, 100))
  
  png("KDE.png")
  plot(kde.bwo, main = "KDE Result with Cross-Validation for Elevation")  
  dev.off()