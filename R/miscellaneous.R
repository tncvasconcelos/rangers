#---------------------

# WWFload is taken from speciesgeocodeR; all credit goes to the original authors
WWFload <- function(x = NULL) {
  if (missing(x)) {
    x <- getwd()
  }
  download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
                destfile = file.path(x, "wwf_ecoregions.zip"), quiet=TRUE)
  unzip(file.path(x, "wwf_ecoregions.zip"), exdir = file.path(x, "WWF_ecoregions"))
  file.remove(file.path(x, "wwf_ecoregions.zip"))
  wwf <- maptools::readShapeSpatial(file.path(x, "WWF_ecoregions", "official",
                                              "wwf_terr_ecos.shp"))
  return(wwf)
}


localityToBiome <- function (points,lat="lat",lon="lon") {
  #colnames(points) <- c("acceptedScientificName","key","decimalLatitude","decimalLongitude","basisOfRecord","issues")
  cat("Getting biome from locality data...")
  points[,lat] <-  as.numeric(points[,lat])
  points[,lon] <-  as.numeric(points[,lon])
  locations.spatial <- sp::SpatialPointsDataFrame(coords=points[,c(lon, lat)], data=points)
  wwf <- WWFload(tempdir())
  mappedregions <- sp::over(locations.spatial, wwf)
  biomes <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga", "Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra", "Mediterranean Forests, Woodlands & Scrub", "Deserts & Xeric Shrublands", "Mangroves")
  points$eco_name <- mappedregions$ECO_NAME
  points$biome <- biomes[mappedregions$BIOME]
  return(points)
}


# getting biomes for each species
getBiomes <- function (points, species="species") {
  cat("Summarizing biome from locality data...")
  points <- as.data.frame(points) # not sure how to do it without transforming back to data.frame
  points <- subset(points, !is.na(points[,"biome"]))
  categories <- unique(points[,"biome"])
  taxa <- as.character(unique(points[,species]))
  result <- matrix(0, nrow=length(taxa), ncol=length(categories))
  rownames(result) <- taxa
  colnames(result) <- categories
  for (taxon_index in seq_along(taxa)) {
    for (category_index in seq_along(categories)) {
      x0 <- points[,species]==taxa[taxon_index]
      x1 <- points[,"biome"]==categories[category_index]
      result[taxon_index, category_index] <- length(which(x0 & x1))
    }
  }
  return(result)
}


#' Extract climate information from species points.
#' @description Function to extract climate information from species points.
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param res A number (2.5, 5 or 10) indicating resolution of climatic layers from where the climate data is being extracted.
#' @importFrom raster extract getData addLayer
#' @importFrom sp coordinates
#' @export
ClimateFromPoints <- function (points, layer, species="", lat = "lat", lon="lon") {
  tmp_points = points
  colnames(tmp_points)[which(colnames(tmp_points) == lon)] <- "lon"
  colnames(tmp_points)[which(colnames(tmp_points) == lat)] <- "lat"
  colnames(tmp_points)[which(colnames(tmp_points) == species)] <- "species"
  tmp_points <- tmp_points[,c("species","lat","lon")]
  cat("Extracting climatic information of", nrow(tmp_points), "points",  "\n")
  sp::coordinates(tmp_points) <- ~ lon + lat
  values <- raster::extract(layer, tmp_points)
  result <- cbind(tmp_points, values)
  return(as.data.frame(result))
}

#' Thinning distribution data to smooth sampling bias
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param n A number indicating how many points to keep in each cell after thinning
#' @importFrom raster extend extent
#' @importFrom dismo gridSample
#' @importFrom sp coordinates
#' @export
Thinning <- function(points, species="species", lat = "decimalLatitude", lon="decimalLongitude", n = 1) {
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  spp <- unique(tmp_points[,species])
  results <- list()
  for(species_index in 1:length(spp)) {
    coords <- tmp_points[tmp_points[,species]==spp[species_index],c("y","x")]
    coords <- coords[!duplicated(coords[,"x"]) & !duplicated(coords[,"y"]),]
    if(nrow(coords) > 1) {
      sp::coordinates(coords) <- ~ y + x
      raster::crs(coords) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      r0 <- raster::raster(coords)
      raster::res(r0) <- 1 # cell resolution
      r0 <- raster::extend(r0, raster::extent(r0) + 5) # expand the extent of the RasterLayer a little
      res <- cbind(spp[species_index], as.data.frame(dismo::gridSample(coords, r0, n))) # n = maximum number of points per cell
      colnames(res) <- c("species", "lat","lon")
      results[[species_index]] <- res
    } else {
      res <- cbind(spp[species_index],coords)
      colnames(res) <- c("species", "lat","lon")
      results[[species_index]] <- res
    }
  }
  results <- do.call(rbind, results)
  return(results)
}

#' Get summary statistics for climatic data of a list of species based on their distribution points
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param res A number (2.5, 5 or 10) indicating resolution of climatic layers from where the climate data is being extracted.
#' @importFrom raster getData extract
#' @importFrom dismo gridSample
#' @importFrom sp coordinates
#' @return A list with summaru statistics for the 19 WorldClim variables plut altitude
#' @export
GetClimateSummStats <- function (points, layer, species="species", lat = "lat", lon="lon") {
  tmp_points = points
  colnames(tmp_points)[which(colnames(tmp_points) == lon)] <- "lon"
  colnames(tmp_points)[which(colnames(tmp_points) == lat)] <- "lat"
  colnames(tmp_points)[which(colnames(tmp_points) == species)] <- "species"
  tmp_points <- tmp_points[,c("species","lon","lat")]
  #cat("Extracting climatic information of", nrow(points), "points",  "\n")
  spp <- unique(tmp_points$species)
  summ_stats <- as.data.frame(matrix(nrow=length(spp), ncol=5))
  for(species_index in 1:length(spp)){
    sp1 <- tmp_points[tmp_points$species==spp[species_index],]
    cat("\r", species_index, "   out of ", length(spp))
    sp::coordinates(sp1) <- ~ lon + lat
    values <- raster::extract(layer, sp1)
    values <- values[!is.na(values)]
    #if(length(values) > 2) {
    n0 <- as.numeric(length(values))
    mean0 <- as.numeric(round(mean(values), 2))
    sd0 <- as.numeric(round(stats::sd(values), 2))
    se0 <- as.numeric(round(sd0/ sqrt(n0), 2))
    tmp_summ_stats <- c(n0, mean0, sd0, se0)
    #} else {
    #  tmp_summ_stats <- rep("not_enough_points", 4)
    #}
    summ_stats[species_index,] <- c(spp[species_index], tmp_summ_stats)
    colnames(summ_stats) <- c("species","n", "mean", "sd", "se")
  }
  return(summ_stats)
}


#convert_to_f <- function(temp){
#  return(((temp/10) * 1.8) + 32)
#}


#' Get summary statistics for four variables (BIO1, BIO4, BIO12, BIO15) based on a list of ranges
#' Those variable represent: annual temperature, temperature seasonality, annual precipitation and
#' precipitation seasonality
#' @param points A list of ranges named as different species
#' @param res Resolution of the climatic layers to be used (a numeric: 2.5, 5 or 10)
#' @importFrom raster getData extract
#' @return A data.frame with means and standard deviations for each climatic variable
#' @export
GetClimate_current <- function (ranges, res=2.5) {
  bio <- raster::getData("worldclim", var="bio", res=res)
  temp <- bio[[1]]
  prec <- bio[[12]]
  seas_temp <- bio[[4]]
  seas_prec <- bio[[15]]
  results <- as.data.frame(matrix(ncol=9, nrow=length(ranges)))
  for(range_index in 1:length(ranges)) {
    cat("\r", "Now doing species number", range_index, "out of", length(ranges))
    points <- as.data.frame(raster::rasterToPoints(ranges[[range_index]])[,1:2])
    if(nrow(points) > 0 & ncol(points) > 1) {
      sp::coordinates(points) <- ~ x + y
      temp_values <- raster::extract(temp, points)
      prec_values <- raster::extract(prec, points)
      seas_temp_values <- raster::extract(seas_temp, points)
      seas_prec_values <- raster::extract(seas_prec, points)
      results[range_index,1] <- names(ranges)[range_index]
      results[range_index,2] <- round(mean(temp_values), 2)
      results[range_index,3] <- round(sd(temp_values), 2)
      results[range_index,4] <- round(mean(prec_values), 2)
      results[range_index,5] <- round(sd(prec_values), 2)
      results[range_index,6] <- round(mean(seas_temp_values), 2)
      results[range_index,7] <- round(sd(seas_temp_values), 2)
      results[range_index,8] <- round(mean(seas_prec_values), 2)
      results[range_index,9] <- round(sd(seas_prec_values), 2)
    }
  }
  colnames(results) <- c("species", "temp_mean", "temp_sd", "prec_mean","prec_sd","seas_temp_mean","seas_temp_sd","seas_prec_mean","seas_prec_sd")
  try(unlink("wc2-5", recursive = TRUE))
  try(unlink("wc5", recursive = TRUE))
  try(unlink("wc10", recursive = TRUE))
  return(results)
}

#' Get summary statistics for the temporal instability of four variables (BIO1, BIO4, BIO12, BIO15) based on a list of ranges
#' Those variable represent: annual temperature, temperature seasonality, annual precipitation and
#' precipitation seasonality. Means and standard deviations are taken from a layer that represent the
#' pairwise comparisons among five paleoclim layers.
#' @param points A list of ranges named as different species
#' @param res Resolution of the climatic layers to be used (a character: "2_5m,"5m" or "10m")
#' @importFrom raster getData extract stack
#' @importFrom rpaleoclim paleoclim
#' @return A data.frame with means and standard deviations for each climatic variable
#' @export
GetClimate_instability <- function (ranges, res="2_5m") {
  # set arguments for getting paleoclim layers
  resolution=res
  keep.vars <-  c("bio_1", "bio_4", "bio_12", "bio_15")
  time.slices <- c("cur","lgm","lig","mpwp","m2")
  # load four paleoclim layers: bio_1, bio_4, bio_12 and bio_15
  paleo_layers <- list()
  for(time_index in 1:length(time.slices)){
    tmp_layer <- rpaleoclim::paleoclim(time.slices[time_index], resolution=resolution)
    tmp_list <- list()
    for(var_index in 1:length(keep.vars)){
      tmp_list[[var_index]] <- tmp_layer[[which(names(tmp_layer) %in% keep.vars[var_index])]]
    }
    paleo_layers[[time_index]]  <- raster::stack(tmp_list)
    names(paleo_layers)[time_index] <- time.slices[time_index]
  }
  # now get instability layers (it can take a while)
  instability_layers <- GetClimateInstability_layers(paleo_layers) # this takes a while
  # get summary stats from those layers for all species
  temp <- instability_layers[[1]]
  prec <- instability_layers[[3]]
  seas_temp <- instability_layers[[2]]
  seas_prec <- instability_layers[[4]]
  results <- as.data.frame(matrix(ncol=9, nrow=length(ranges)))
  for(range_index in 1:length(ranges)) {
    cat("\r", "Now doing species number", range_index, "out of", length(ranges))
    points <- as.data.frame(rasterToPoints(ranges[[range_index]])[,1:2])
    sp::coordinates(points) <- ~ x + y
    temp_values <- raster::extract(temp, points)
    prec_values <- raster::extract(prec, points)
    seas_temp_values <- raster::extract(seas_temp, points)
    seas_prec_values <- raster::extract(seas_prec, points)
    results[range_index,1] <- names(ranges)[range_index]
    results[range_index,2] <- round(mean(temp_values[!is.na(temp_values)]), 2)
    results[range_index,3] <- round(sd(temp_values[!is.na(temp_values)]), 2)
    results[range_index,4] <- round(mean(prec_values[!is.na(prec_values)]), 2)
    results[range_index,5] <- round(sd(prec_values[!is.na(prec_values)]), 2)
    results[range_index,6] <- round(mean(seas_temp_values[!is.na(seas_temp_values)]), 2)
    results[range_index,7] <- round(sd(seas_temp_values[!is.na(seas_temp_values)]), 2)
    results[range_index,8] <- round(mean(seas_prec_values[!is.na(seas_prec_values)]), 2)
    results[range_index,9] <- round(sd(seas_prec_values[!is.na(seas_prec_values)]), 2)
  }
  colnames(results) <- c("species", "inst_temp_mean", "inst_temp_sd", "inst_prec_mean","inst_prec_sd","inst_seas_temp_mean","inst_seas_temp_sd","inst_seas_prec_mean","inst_seas_prec_sd")
  # unlink folders with paleoclim layers
  try(unlink("./wc2-5", recursive = TRUE))
  try(unlink("./wc5", recursive = TRUE))
  try(unlink("./wc10", recursive = TRUE))
  return(results)
}


#' @param all.rasters A list of ranges named as different species
#' @importFrom raster getData extract stack extent crop calc
#' @importFrom rdist rdist
GetClimateInstability_layers <- function(all.rasters) {
  if(var(unlist(lapply(all.rasters, nlayers))) != 0) {
    stop("elements of the list of layers don't have the same length")
  } else {
    xmin <- max(unlist(lapply(lapply(lapply(all.rasters, raster::extent), unextent), "[[", 1)))
    xmax <- min(unlist(lapply(lapply(lapply(all.rasters, raster::extent), unextent), "[[", 2)))
    ymin <- max(unlist(lapply(lapply(lapply(all.rasters, raster::extent), unextent), "[[", 3)))
    ymax <- min(unlist(lapply(lapply(lapply(all.rasters, raster::extent), unextent), "[[", 4)))
    to_crop <- c(xmin, xmax, ymin, ymax)
    instability_rasters <- list()
    #all.rasters <- lapply(all.rasters, crop, to_crop)
    for(bio_index in 1:lapply(all.rasters, nlayers)[[1]]){
      cat("\r", "Calculating instability layer...", bio_index, "in", lapply(all.rasters, nlayers)[[1]])
      r1 <- raster::crop(all.rasters[[1]][[bio_index]], to_crop)
      r2 <- raster::crop(all.rasters[[2]][[bio_index]], to_crop)
      tmp_pairwise <- raster::calc(stack(r1, r2), rdist::rdist)
      for(slice_index in 3:length(all.rasters)) {
        tmp_r <- crop(all.rasters[[slice_index]][[bio_index]], to_crop)
        tmp_pairwise <- raster::calc(stack(tmp_pairwise, tmp_r), rdist::rdist)
      }
      instability_rasters[[bio_index]] <- tmp_pairwise
      names(instability_rasters)[bio_index] <- names(r1)
      print(bio_index)
    }
  }
  cat("","\n")
  return(instability_rasters)
}

unextent <- function(x) {
  return(c(x[1], x[2], x[3], x[4]))
}

#' Get summary statistics for the spatial heterogeneity of four variables (BIO1, BIO4, BIO12, BIO15) based on a list of ranges
#' Those variable are: annual temperature, temperature seasonality, annual precipitation and
#' precipitation seasonality. Means and standard deviations are taken from a layer of slopes for each variable
#' @param points A list of ranges named as different species
#' @param res Resolution of the climatic layers to be used (a numeric: 2.5, 5 or 10)
#' @importFrom raster terrain extract getData
#' @return A data.frame with means and standard deviations for each climatic variable
#' @export
GetClimate_heterogeneity <- function (ranges, res=2.5) { #buffer=100000,
  bio <- raster::getData("worldclim", var="bio", res=res)
  temp_slope <- raster::terrain(bio[[1]])
  prec_slope <- raster::terrain(bio[[12]])
  seas_temp_slope <- raster::terrain(bio[[4]])
  seas_prec_slope <- raster::terrain(bio[[15]])
  results <- as.data.frame(matrix(ncol=9, nrow=length(ranges)))
  for(range_index in 1:length(ranges)) {
    cat("\r", "Now doing species number", range_index, "out of", length(ranges))
    points <- as.data.frame(rasterToPoints(ranges[[range_index]])[,1:2])
    if(nrow(points) > 0 & ncol(points) > 1) {
      sp::coordinates(points) <- ~ x + y
      #circle_around_point <- dismo::circles(points, d=buffer, lonlat=TRUE) # Radius of the circle in meters
      #circle_around_point <- sp::polygons(circle_around_point)
      circle_around_point <- points
      temp_values <- unlist(raster::extract(temp_slope, circle_around_point))
      prec_values <- unlist(raster::extract(prec_slope, circle_around_point))
      seas_temp_values <- unlist(raster::extract(seas_temp_slope, circle_around_point))
      seas_prec_values <- unlist(raster::extract(seas_prec_slope, circle_around_point))
      results[range_index,1] <- names(ranges)[range_index]
      results[range_index,2] <- round(mean(temp_values[!is.na(temp_values)]*100000), 2)
      results[range_index,3] <- round(sd(temp_values[!is.na(temp_values)]*100000), 2)
      results[range_index,4] <- round(mean(prec_values[!is.na(prec_values)]*100000), 2)
      results[range_index,5] <- round(sd(prec_values[!is.na(prec_values)]*100000), 2)
      results[range_index,6] <- round(mean(seas_temp_values[!is.na(seas_temp_values)]*100000), 2)
      results[range_index,7] <- round(sd(seas_temp_values[!is.na(seas_temp_values)]*100000), 2)
      results[range_index,8] <- round(mean(seas_prec_values[!is.na(seas_prec_values)]*100000), 2)
      results[range_index,9] <- round(sd(seas_prec_values[!is.na(seas_prec_values)]*100000), 2)
    }
  }
  colnames(results) <- c("species", "slope_temp_mean", "slope_temp_sd", "slope_prec_mean","slope_prec_sd","slope_seas_temp_mean","slope_seas_temp_sd","slope_seas_prec_mean","slope_seas_prec_sd")
  # unlink folders with paleoclim layers
  try(unlink("./wc2-5", recursive = TRUE))
  try(unlink("./wc5", recursive = TRUE))
  try(unlink("./wc10", recursive = TRUE))
  return(results)
}


#' Get proportion of cells in the range of a species that experience some frost and a lot of frost every year.
#' @param points A list of ranges named as different species
#' @param res Resolution of the climatic layers to be used (a numeric: 2.5, 5 or 10)
#' @importFrom raster extract getData
#' @return A data.frame with proportions of some frost and a lot of frost for each species
#' @export
Frost <- function(ranges, res=2.5) {
  bio <- raster::getData("worldclim", var="bio", res=res)
  some_frost <- bio[[6]]
  a_lot_of_frost <- bio[[11]]
  a_lot_of_frost[a_lot_of_frost[] < 0] <- 0
  a_lot_of_frost[a_lot_of_frost[] > 0] <- 1
  some_frost[some_frost[] < 0] <- 0
  some_frost[some_frost[] > 0] <- 1
  some_but_not_a_lot_of_frost <- calc(stack(a_lot_of_frost,some_frost), diff)
  some_but_not_a_lot_of_frost[] <- abs(some_but_not_a_lot_of_frost[])
  results <- as.data.frame(matrix(ncol=3, nrow=length(ranges)))
  for(range_index in 1:length(ranges)) {
    cat("\r", "Now doing species number", range_index, "out of", length(ranges))
    points <- as.data.frame(rasterToPoints(ranges[[range_index]])[,1:2])
    if(nrow(points) > 0 & ncol(points) > 1) {
      sp::coordinates(points) <- ~ x + y
      some_frost_values <- raster::extract(some_but_not_a_lot_of_frost, points)
      a_lot_of_frost_values <- raster::extract(a_lot_of_frost, points)
      some_frost_values <- table(some_frost_values)
      
      if(any(names(some_frost_values) == 1)){
        some_frost_values <- some_frost_values[which(names(some_frost_values) == 1)]
      } else {
        some_frost_values <- 0 }
      
      a_lot_of_frost_values <- table(a_lot_of_frost_values)
      if(any(names(a_lot_of_frost_values) == 0)){
        a_lot_of_frost_values <- a_lot_of_frost_values[which(names(a_lot_of_frost_values) == 0)]
      } else {
        a_lot_of_frost_values <- 0 }
      
      proportion_some_frost <- some_frost_values / length(points)
      proportion_a_lot_of_frost <- a_lot_of_frost_values / length(points)
      results[range_index,1] <- names(ranges)[range_index]
      results[range_index,2] <- round(proportion_some_frost, 2)
      results[range_index,3] <- round(proportion_a_lot_of_frost, 2)
    }
  }
  colnames(results) <- c("species", "some_frost", "a_lot_of_frost")
  # unlink folders with paleoclim layers
  try(unlink("./wc2-5", recursive = TRUE))
  try(unlink("./wc5", recursive = TRUE))
  try(unlink("./wc10", recursive = TRUE))
  return(results)
}


#' Because sometimes taxonomic changes are faster than the rate online databases can update them
#' And sometimes authors like to abreviate genera and species names in... weird ways.
sub_tiplabels <- function(phy, table) {
  new_tip_labels <- c()
  for(tip_index in 1:length(phy$tip.label)) {
    one_old_tip <- phy$tip.label[tip_index]
    new_tip_labels[tip_index] <- table[,2][which(table[,1] == one_old_tip)]
  }
   phy$tip.label <- new_tip_labels
   return(phy)
}

