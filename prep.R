library(dplyr)
library(raster)
library(rgeos)
library(maptools)
#library(ggplot2)
library(rgdal)
library(gdalUtils)

setwd("~/Work/katie_myrt")

rm(list = ls())

source("R/globals.R")
source("R/range.R")

### get Myrtle Rust rasters

mrCurrent <- raster("data/myrtle_rust/raster/PA_Current.gri") %>%
  projectRaster(crs = aea, method = "ngb")

mrFuture <- raster("data/myrtle_rust/raster/PA_Future.gri") %>%
  projectRaster(crs = aea, method = "ngb")

# change raster values to NA if 0 (we only want values 1 and NA)
mrCurrent[mrCurrent == 0] <- NA
mrFuture[mrFuture == 0] <- NA

# write output rasters
writeRaster(mrCurrent, "output/myrtle_rust/mr_current.tif", overwrite = TRUE)
writeRaster(mrFuture, "output/myrtle_rust/mr_future.tif", overwrite = TRUE)


### get occurrence data

data <- read.csv("data/occurrence/New_Myrt_Aust_Complete_270916.csv", stringsAsFactors = FALSE) %>%
  setNames(c("id", "species", "latitude", "longitude")) %>%
  filter(!is.na(longitude) & !is.na(latitude))
write.csv(data, "output/occurrence/occurrence.csv", row.names = FALSE)

# get Australia boundary
aus <- readOGR("data/gis/australia/SHAPEFILE.shp", layer = "SHAPEFILE") %>%
  spTransform(aea)


### store taxon spatial data

# get list of taxa
taxa <- unique(data$species) %>% sort

mcp <- vector("list")

for (taxon in taxa) {
  message(taxon)
  
  # get occurrence data for taxon
  d <- data %>%
    filter(species == taxon) %>%
    dplyr::select(longitude, latitude) %>%
    SpatialPoints(CRS("+proj=longlat +datum=WGS84")) %>%
    spTransform(aea)
  
  # compute MCP
  mcpFull <- mcpPolygon(d@coords)

  # get number of vertices in MCP
  numVertices <- polygonToDataFrame(mcpFull) %>% unique %>% nrow  
  
  # intersect MCP with aus border; can only do intersection if 4 or more vertices in MCP
  if (numVertices >= 4) {
    mcpCrop <- gIntersection(mcpFull, aus)
  } else {
    mcpCrop <- mcpFull
  }  
  
  # create data frame of cropped MCP vertices
  mcpDF <- polygonToDataFrame(mcpCrop)  
  
  # write to list
  mcp[[taxon]]$mcp_full <- mcpFull
  mcp[[taxon]]$num_vertices <- numVertices
  mcp[[taxon]]$mcp_cropped <- mcpCrop
  mcp[[taxon]]$mcp_cropped_df <- mcpDF

}

saveRDS(mcp, "output/mcp.rds")

###

mcp <- readRDS("output/mcp.rds")

mrCurrent <- raster("output/myrtle_rust/mr_current.tif")
mrFuture <- raster("output/myrtle_rust/mr_future.tif")

tempShp <- "temp/temp.shp"

for (taxon in taxa[20]) {
  message(taxon)
  
  # generate base raster for MCP and 
  # only need to have raster extent of MCP polygon, no need to have full Aus extent  ===== OBSOLETE???
  e <- alignExtent(extent(mcp[[taxon]]$mcp_cropped), mrCurrent)
  baseRaster <- raster(ext = e, resolution = res(mrCurrent), crs = aea)
  
  mcp[[taxon]]$mcp_cropped %>%
    as("SpatialPolygonsDataFrame") %>%
    writeOGR(dirname(tempShp), gsub("[.]shp", "", basename(tempShp)), driver = "ESRI Shapefile", overwrite = TRUE)
  
  writeRaster(baseRaster, "temp/temp.tif", overwrite = TRUE)
  
  mcpRast <- gdal_rasterize(b = 1, burn = 1, l = "temp", src_datasource = tempShp, dst_filename = "temp/temp.tif", output_Raster = TRUE) %>%
    raster(layer = 1)
  
  # method using MUCH slower raster::rasterize()
  #mcpRast <- rasterize(mcp[[taxon]]$mcp_cropped, baseRaster)
  
  # write cropped raster to file
  writeRaster(mcpRast, paste0("output/mcp/", gsub(" ", "_", taxon), ".tif"), overwrite = TRUE)
  
  # compute raster of overlap between MCP and MR
  x <- crop(mrCurrent, mcpRast)
  out <- raster::mask(mcpRast, x, maskvalue = NA)

  writeRaster(out, paste0("output/mcp_mr/", gsub(" ", "_", taxon), ".tif"), overwrite = TRUE)
}



