plotTaxonEOO <- function(d, outline, mcp, zoomToOccurrences = TRUE) {
  
  xlimits <- c(min(d$longitude), max(d$longitude))
  ylimits <- c(min(d$latitude), max(d$latitude))
  
  mcpDF <- polygonToDataFrame(mcp)
  
  p <- ggplot() +
    geom_path(data = outline, aes(x = longitude, y = latitude, group = id), colour = "antiquewhite") +
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
    list %>% SpatialPolygons(proj4string = CRS("+proj=longlat +datum=WGS84"))
}