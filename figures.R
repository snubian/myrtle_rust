library(dplyr)
library(raster)
library(rgeos)
#library(maptools)
library(ggplot2)
library(rgdal)

setwd("~/Work/katie_myrt")

rm(list = ls())

source("R/globals.R")
source("R/range.R")

# get occurrence data
data <- read.csv("output/occurrence/occurrence.csv", stringsAsFactors = FALSE)

# get Australia boundary
aus <- readOGR("data/gis/australia/SHAPEFILE.shp", layer = "SHAPEFILE") %>%
  spTransform(aea)

# get Mrytle Rust rasters
mrCurrent <- raster("output/myrtle_rust/mr_current.tif")
mrFuture <- raster("output/myrtle_rust/mr_future.tif")

# compute raster cell area
cellAreaKm2 <- (mrCurrent %>% res %>% prod) / 1000000

# read MCP data from file
mcp <- readRDS("output/mcp.rds")

# get list of taxa
taxa <- unique(data$species) %>% sort

y <- mrCurrent %>%
  as("SpatialPixelsDataFrame") %>%
  as.data.frame

for (taxon in taxa) {
  
  # get occurrence data for taxon
  d <- data %>%
    filter(species == taxon) %>%
    dplyr::select(longitude, latitude) %>%
    SpatialPoints(CRS("+proj=longlat +datum=WGS84")) %>%
    spTransform(aea)
    
  # plot MR map
  plotTaxonMR(data.frame(d@coords), aus, mcp[[taxon]]$mcp_cropped, y)  
  ggsave(filename = paste0("figs/range_mr_maps/", gsub(" ", "_", taxon), ".jpg"), width = 8.01, height = 5.82)

  # plot EOO maps of range
  plotTaxonEOO(data.frame(d@coords), aus, mcp[[taxon]]$mcp_cropped, zoomToOccurrences = TRUE)
  ggsave(filename = paste0("~/Work/katie_myrt/figs/range_maps/", gsub(" ", "_", taxon), ".jpg"), width = 8.01, height = 5.82)
  
}