plotTaxonMR <- function(d, outline, mcp, y) {
  
  mcpDF <- polygonToDataFrame(mcp)  
  
  p <- ggplot() +
    geom_polygon(data = outline, aes(x = long, y = lat, group = group), fill = NA, colour = "grey") +
    geom_tile(data = y, aes(x = x, y = y, fill = mr_current), colour = "lightgoldenrod2")
  
  #if (length(out[out == 1]) > 0) {
  #  z <- out %>%
  #    as("SpatialPixelsDataFrame") %>%
  #    as.data.frame  
  #  
  #  p <- p + geom_tile(data = z, aes(x = x, y = y, fill = layer), colour = "tomato")
  #}
  
  p <- p +
    geom_polygon(data = mcpDF, aes(x = longitude, y = latitude, group = id), fill = NA, colour = "springgreen4") +
    xlab("") +
    ylab("") +
    #ggtitle(paste0(taxon, " [", round(proportionMR, 4), "]")) +
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
  
}



plotTaxonEOO <- function(d, outline, mcp, zoomToOccurrences = TRUE) {
  
  xlimits <- c(min(d$longitude)*1.01, max(d$longitude)*0.99)
  ylimits <- c(min(d$latitude)*1.01, max(d$latitude)*0.99)
  
  mcpDF <- polygonToDataFrame(mcp)
  
  p <- ggplot() +
    geom_path(data = outline, aes(x = long, y = lat, group = group), colour = "antiquewhite") +
    geom_point(data = d, aes(x = longitude, y = latitude), size = 2, alpha = 0.3, colour = "red") +
    geom_polygon(data = mcpDF, aes(x = longitude, y = latitude, group = id), colour = "steelblue1", fill = "steelblue4", alpha = 0.05) +    
    xlab("") +
    ylab("") +
    ggtitle(taxon) +
    labs(colour = "") +
    theme_bw() +
    scale_colour_manual(values = c("red", "green"), drop = FALSE) +
    xlim(xlimits) +
    ylim(ylimits) +
    coord_equal() +
    theme(plot.title = element_text(face = "bold", size = 14, margin = margin(0, 0, 30, 0)),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 14))

  p
}

polygonToDataFrame <- function(p) {
  # create dataframe suitable for ggplot from polygon data  
  id <- 1
  for (i in 1:length(p@polygons)) {
    for (j in 1:length(p@polygons[[i]]@Polygons)) {
      d <- data.frame(p@polygons[[i]]@Polygons[[j]]@coords)
      d$id <- id
      colnames(d) <- c("longitude", "latitude", "id")
      if (id == 1) {
        out <- d
      } else {
        out <- rbind(out, d)
      }
      id <- id + 1
    }
  }
  return(out)  
}

mcpPolygon <- function(d) {
  d %>%
    chull %>%
    c(.[1]) %>%
    d[., ] %>%
    Polygon %>%
    list %>% Polygons(1) %>%
    list %>% SpatialPolygons(proj4string = CRS(aea))
}