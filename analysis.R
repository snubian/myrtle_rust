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

# setup output data frame
output <- data.frame(matrix(ncol = 4, nrow = length(taxa)))
colnames(output) <- c("species", "cells_in_eoo", "cells_in_mr_current", "cells_in_mr_future")

for (i in 1:length(taxa)) {
  taxon <- taxa[i]
  
  message(taxon)
  
  # get occurrence data for taxon
  d <- data %>%
    filter(species == taxon) %>%
    dplyr::select(longitude, latitude) %>%
    SpatialPoints(CRS("+proj=longlat +datum=WGS84")) %>%
    spTransform(aea)
  
  # get EOO cropped MCP
  mcpRast <- raster(paste0("output/mcp/", gsub(" ", "_", taxon), ".tif"))
  
  # get rasters of intersection of MCP and MR (for both current and future)
  mcpMrCurrent <- raster(paste0("output/mcp_mr/current/", gsub(" ", "_", taxon), ".tif"))
  mcpMrFuture <- raster(paste0("output/mcp_mr/future/", gsub(" ", "_", taxon), ".tif"))

  # write results to output data frame
  output[i, "species"] <- taxon
  # get number of cells in EOO MCP; can simply add all cell values as they are 1
  output[i, "cells_in_eoo"] <- sum(mcpRast[], na.rm = TRUE)

  # get number of cells overlapping Myrtle Rust area; both current and future
  output[i, "cells_in_mr_current"] <- sum(mcpMrCurrent[], na.rm = TRUE)
  output[i, "cells_in_mr_future"] <- sum(mcpMrFuture[], na.rm = TRUE)

}

# calculate proportion of MR cells in EOO
output <- output %>%
  mutate(proportion_of_mr_current = cells_in_mr_current / cells_in_eoo,
         proportion_of_mr_future = cells_in_mr_future / cells_in_eoo)

write.csv(output, "output/mr_proportion.csv", row.names = FALSE)
