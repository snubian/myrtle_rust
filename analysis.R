library(dplyr)
library(raster)
library(rgeos)
library(maptools)
library(ggplot2)
library(rgdal)

setwd("~/Work/katie_myrt")

rm(list = ls())

source("R/globals.R")
source("R/range.R")

# get occurrence data
data <- read.csv("data/occurrence/New_Myrt_Aust_Complete_270916.csv", stringsAsFactors = FALSE) %>%
  setNames(c("id", "species", "latitude", "longitude")) %>%
  filter(!is.na(longitude) & !is.na(latitude))

# get Australia boundary
aus <- readOGR("data/gis/australia/SHAPEFILE.shp", layer = "SHAPEFILE") %>%
  spTransform(aea)

# get Mrytle Rust rasters
mrCurrent <- raster("output/myrtle_rust/mr_current.tif")
mrFuture <- raster("output/myrtle_rust/mr_future.tif")

y <- mrCurrent %>%
  as("SpatialPixelsDataFrame") %>%
  as.data.frame

cellAreaKm2 <- (mrCurrent %>% res %>% prod) / 1000000

# get list of taxa
taxa <- unique(data$species) %>% sort

for (taxon in taxa[1:20]) {
  message(taxon)
  
  # get occurrence data for taxon
  d <- data %>%
    filter(species == taxon) %>%
    dplyr::select(longitude, latitude) %>%
    SpatialPoints(CRS("+proj=longlat +datum=WGS84")) %>%
    spTransform(aea)
  
  # compute MCP
  mcp <- mcpPolygon(d@coords)

  # get number of vertices in MCP
  numVertices <- polygonToDataFrame(mcp) %>% unique %>% nrow
  
  # intersect MCP with aus border
  # can only do intersection if 4 or more vertices in MCP
  if (numVertices >= 4) {
    mcpCrop <- gIntersection(mcp, aus)
  } else {
    mcpCrop <- mcp
  }
  
  # create data frame of cropped MCP vertices
  mcpDF <- polygonToDataFrame(mcpCrop)
  
  # convert cropped MCP to raster, using MR raster as template (same extent and resolution)
  mcpRast <- rasterize(mcpCrop, mrCurrent)
  
  # get intersection between two rasters
  out <- raster::mask(mcpRast, mrCurrent, maskvalue = NA)
  
  # get number of cells in MCP; can simply add all cell values as they are 1
  numCellsInRange <- sum(mcpRast[], na.rm = TRUE)  
  
  # get number of cells overlapping Myrtle Rust area
  numCellsInMR <- sum(out[], na.rm = TRUE)
  
  # get proportion of range in Myrtle Rust area
  proportionMR <- numCellsInMR / numCellsInRange 
  
  p <- ggplot() +
    geom_polygon(data = aus, aes(x = long, y = lat, group = group), fill = NA, colour = "grey") +
    geom_tile(data = y, aes(x = x, y = y, fill = layer), colour = "lightgoldenrod2")
  
  if (length(out[out == 1]) > 0) {
    z <- out %>%
      as("SpatialPixelsDataFrame") %>%
      as.data.frame  
    
    p <- p + geom_tile(data = z, aes(x = x, y = y, fill = layer), colour = "tomato")
  }
  
  p <- p +
    geom_polygon(data = mcpCrop, aes(x = long, y = lat, group = group), fill = NA, colour = "springgreen4") +
    xlab("") +
    ylab("") +
    ggtitle(paste0(taxon, " [", round(proportionMR, 4), "]")) +
    guides(fill = FALSE) +
    theme_bw() +
    coord_equal() +
    theme(plot.title = element_text(face = "bold", size = 14, margin = margin(0, 0, 30, 0)),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 14),
          panel.grid = element_blank(),  ## remove grid lines
          axis.text = element_blank(),  ## remove y-axis label text
          axis.ticks = element_blank(),  ## remove y-axis tick marks
          panel.border = element_blank())

  p

  ggsave(filename = paste0("figs/range_mr_maps/", gsub(" ", "_", taxon), ".jpg"), width = 8.01, height = 5.82)
  
  # plot EOO maps of range
  #plotTaxonEOO(d, ausPoly, mcpCrop, zoomToOccurrences = TRUE)
  #ggsave(filename = paste0("~/Work/katie_myrt/figs/range_maps/", gsub(" ", "_", taxon), ".jpg"), width = 8.01, height = 5.82)
}
