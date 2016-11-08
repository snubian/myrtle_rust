library(dplyr)
library(raster)
library(rgeos)
library(maptools)
library(rgdal)
library(gdalUtils)

setwd("~/Work/katie_myrt")

rm(list = ls())

source("R/globals.R")
source("R/range.R")

############################################################################
### prepare Myrtle Rust rasters

# re-project to equal area projection, use nearest neighbour to retain 1 values
mrCurrent <- raster("data/myrtle_rust/raster/PA_Current.gri") %>%
  projectRaster(crs = aea, method = "ngb")

mrFuture <- raster("data/myrtle_rust/raster/PA_Future.gri") %>%
  projectRaster(crs = aea, method = "ngb")

# change raster values to NA if 0 (we only want values 1 and NA)
mrCurrent[mrCurrent == 0] <- NA
mrFuture[mrFuture == 0] <- NA

# write output MR rasters
writeRaster(mrCurrent, "output/myrtle_rust/mr_current.tif", overwrite = TRUE)
writeRaster(mrFuture, "output/myrtle_rust/mr_future.tif", overwrite = TRUE)

############################################################################
### prepare occurrence data

data <- read.csv("data/occurrence/New_Myrt_Aust_Complete_270916.csv", stringsAsFactors = FALSE) %>%
  setNames(c("id", "species", "latitude", "longitude")) %>%
  filter(!is.na(longitude) & !is.na(latitude))
write.csv(data, "output/occurrence/occurrence.csv", row.names = FALSE)

############################################################################
### prepare EOO MCPs for each species

# get Australia boundary
aus <- readOGR("data/gis/australia/SHAPEFILE.shp", layer = "SHAPEFILE") %>%
  spTransform(aea)

# get list of taxa
taxa <- unique(data$species) %>% sort

# create output vector
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
  
  # write to output vector
  mcp[[taxon]]$mcp_full <- mcpFull
  mcp[[taxon]]$num_vertices <- numVertices
  mcp[[taxon]]$mcp_cropped <- mcpCrop
  mcp[[taxon]]$mcp_cropped_df <- mcpDF

}

# write mcp data to file
saveRDS(mcp, "output/mcp.rds")

##############################################################################
### prepare species rasters for EOO (cropped MCP) and intersection of EOO & MR

# get mcp data
mcp <- readRDS("output/mcp.rds")

# get list of taxa
taxa <- unique(data$species) %>% sort

# get MR rasters
mrCurrent <- raster("output/myrtle_rust/mr_current.tif")
mrFuture <- raster("output/myrtle_rust/mr_future.tif")

for (taxon in taxa) {
  message(taxon)
  
  # generate base raster for MCP; only need to have raster extent of MCP polygon, no need to have full Aus extent
  baseRaster <- mcp[[taxon]]$mcp_cropped %>%
    extent %>%
    alignExtent(mrCurrent) %>%
    raster(ext = ., resolution = res(mrCurrent), crs = aea)
  
  # write base raster to temp file; this is used as base for output of rgdal_rasterize
  writeRaster(baseRaster, "temp/temp.tif", overwrite = TRUE)  
  
  # write MCP polygon to temp file
  mcp[[taxon]]$mcp_cropped %>%
    as("SpatialPolygonsDataFrame") %>%
    writeOGR("temp", "temp", driver = "ESRI Shapefile", overwrite = TRUE)
    
  # rasterize MCP
  mcpRast <- gdal_rasterize(b = 1, burn = 1, l = "temp", src_datasource = "temp/temp.shp", dst_filename = "temp/temp.tif", output_Raster = TRUE) %>%
    raster(layer = 1)
  
  # write cropped raster to file
  writeRaster(mcpRast, paste0("output/mcp/", gsub(" ", "_", taxon), ".tif"), overwrite = TRUE)
  
  # compute raster of overlap between MCP and MR; need to first crop MR raster to MCP raster extent
  mrCurrent %>%
    crop(mcpRast) %>%
    raster::mask(mcpRast, ., maskvalue = NA) %>%
    writeRaster(paste0("output/mcp_mr/current/", gsub(" ", "_", taxon), ".tif"), overwrite = TRUE)

  mrFuture %>%
    crop(mcpRast) %>%
    raster::mask(mcpRast, ., maskvalue = NA) %>%
    writeRaster(paste0("output/mcp_mr/future/", gsub(" ", "_", taxon), ".tif"), overwrite = TRUE)

}



