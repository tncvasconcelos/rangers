#' A function to overlay distributions and create a map of species richness
#' @param list_of_ranges A list in the format of the output of GetRanges()
#' @return A raster of species richness
#' @importFrom raster resample calc stack mask
#' @export
GetSpRichness <- function (list_of_ranges) {
  ranges <- unlist(lapply(list_of_ranges, "[[", "range"))
  #template.map=NULL
  #if(is.null(template.map)) {
    template.map <- readRDS("data/template.map.Rdata")
    #template.map <- raster::getData("worldclim", var="bio", download=TRUE, res=10)[[1]]
    #template.map[!is.na(template.map)] <- 0
  #} else { template.map=template.map }
  tmp.raster.list <- list()
  for (i in 1:length(ranges)) {
    r1 <- ranges[[i]]
    r1 <- raster::resample(r1, template.map)
    r1[is.na(r1)] <- 0
    tmp.raster.list[[i]] <- raster::mask(r1, template.map)
    print(i)
  }
  names(tmp.raster.list) <- names(ranges)
  sprichness_map <- raster::calc(raster::stack(tmp.raster.list), sum)
  #saveRDS(tmp.raster.list, file=paste0("~/Desktop/MiSSEgradient/MiSSEGradient/regressions_plan/Data/1_rate_rasters/all_species_stack.Rdata"))
  return(sprichness_map)
}

#' A function to map distribution of traits (still with no PW)
#' @param list_of_ranges A list in the format of the output of GetRanges()
#' @param trait_data A data.frame where the first column contain species names and the second contains trait data
#' @param type "binary" or "continuous", depending on the kind of trait data
#' @return A raster with the traits mapped in space
#' @importFrom raster resample calc stack mask
#' @export
GetTraitDistribution <- function (list_of_ranges, trait_data, type=c("binary","continuous")) {
  if(length(list_of_ranges) != nrow(trait_data)) {
    stop("Length of list_of_ranges and trait_data do not match.")
  }
  if(any(!names(list_of_ranges) %in% trait_data[,1])) {
    stop("Some species in the list_of_ranges are not in the trait_data")
  }
  ranges <- lapply(list_of_ranges, function(x) x$range)
  template.map <- readRDS("R/template.map.Rdata")
  #template.map <- raster::getData("worldclim", var="bio", download=TRUE, res=10)[[1]]
  #template.map[!is.na(template.map)] <- 0
  if(type=="continuous") {
    tmp.raster.list_traits <- list()
    tmp.raster.list_sprich <- list()
    for (range_index in 1:length(ranges)) {
      species <- names(ranges)[[range_index]]
      r1 <- ranges[[range_index]]
      r1 <- raster::resample(r1, template.map)
      r1[is.na(r1)] <- 0
      tmp.raster.list_sprich[[range_index]] <- raster::mask(r1, template.map)
      r1[r1==1] <- trait_data[trait_data$species==species, 2]
      tmp.raster.list_traits[[range_index]] <- raster::mask(r1, template.map)
    }
    traits_sum <- raster::calc(raster::stack(tmp.raster.list_traits), sum)
    sp_sum <- raster::calc(raster::stack(tmp.raster.list_sprich), sum)
    result <- r1
    result[!is.na(result)] <- traits_sum / sp_sum
    return(mask(result, template.map))

  } else if (type=="binary") {
    trait_states_list <- list()
    states <- unique(trait_data[,2])
    if(length(states) != 2) {
      stop("Trait is not binary")
    }
    #trait_data[trait_data[,2]==states[1],2] <- 0
    #trait_data[trait_data[,2]==states[2],2] <- 1
    for(state_index in 1:length(states)) {
      t0 <- trait_data[which(trait_data[,2]==states[state_index]),]
      r0 <- ranges[which(names(ranges) %in% t0[,1])]
      tmp.raster.list <- list()
      for (range_index in 1:length(r0)) {
        r1 <- r0[[range_index]]
        r1 <- raster::resample(r1, template.map)
        r1[is.na(r1)] <- 0
        tmp.raster.list[[range_index]] <- raster::mask(r1, template.map)
      }
      trait_states_list[[state_index]] <- calc(stack(tmp.raster.list), sum)
      names(trait_states_list)[state_index] <- states[state_index]
    }
    total_sprich <- raster::calc(raster::stack(trait_states_list), sum)
    test <- total_sprich
    proportion_trait1 <- ((trait_states_list[[1]] * 100) / total_sprich) / 100
    proportion_trait2 <- ((trait_states_list[[2]] * 100) / total_sprich) / 100
    final <- c(proportion_trait1, proportion_trait2)
    names(final) <- names(trait_states_list)
    return(final)
  }
}

#remove_data_raster <- function(x) {
#  x[x[]==0] <- NA
#  x[!is.na(x)] <- 0
#  x
#}

#range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#' A function to map distribution of traits
#' @param raster1 A raster of mapped trait distribution or species-richness
#' @param raster2 A raster of mapped trait distribution or species-richness
#' @return A raster with mapped residuals of a linear regression between raster1 and raster2
#' @importFrom raster getValues crop extent
#' @importFrom stats lm na.exclude residuals.lm
#' @export
GetResDistribution <- function(raster1, raster2) {
  #if(is.null(template.map)) {
    template.map <- readRDS("R/template.map.Rdata")
    #template.map <- raster::getData("worldclim", var="bio", download=TRUE, res=10)[[1]]
    #template.map[!is.na(template.map)] <- 0
  #} else { template.map=template.map }
  template <- crop(template.map, raster::extent(-180, 180, -60, 90))
  # set pallete
  raster1[raster1[]==0] <- NA
  raster2[raster2[]==0] <- NA
  l.model <- stats::lm(raster::getValues(raster1) ~ raster::getValues(raster2), na.action = na.exclude)
  res.raster <- template
  res.raster[] <- as.numeric(stats::residuals.lm(l.model))
  return(res.raster)
}

