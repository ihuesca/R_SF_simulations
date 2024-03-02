#####################################################################
# Residence and Site Loyalty Indices
#####################################################################

# Note: The indices were adjusted to our case study.

# Ocurrency, Permanency, and Periodicity
# Indices based on Ballance (1990) 

# Ocurrency
Occurrence <- function(x) {
  output <- apply(x, 1, function(x) sum(x, na.rm = T))
  return(output)
}

# Permanency
Permanency <- function(x) {
  output <- apply(x, 1, function(x) {
    out <- 0
    if (sum(x, na.rm = T) > 0) {
      out <- max(which(x == 1)) - min(which(x == 1))
    }
    return(out)
  })
  return(output)
}

# Periodicity
Periodicity <- function(x, ...) {
  output <- apply(x, 1, function(x) {
    sighting <- which(x > 0)
    out <- 0
    if (length(sighting) > 1) {
      out <- (length(sighting) - 1)/sum(diff(sighting))
    }
    return(out)
  })
  return(output)
}

# Index based on Parra et al. (2006)

index_Parra <- function(x, year, ...) {
  output <- apply(x, 1, function(x) {
    MSR <- sum(x)/length(x)
    years <- rep(1:round(length(x)/year, 0), each = year, length.out = length(x))
    sighting.year <- aggregate(x, by = list(years), FUN = function(x) sum(x, na.rm = T))$x
    YSR <- length(sighting.year[sighting.year > 0])/length(sighting.year)
    out <- cbind(MSR = MSR, YSR = YSR)
    return(out)
  })
  output <- data.frame(MSR = output[1, ], YSR = output[2, ])
  return(output)
}

# Index based on Díaz-López (2012)

index_Diaz <- function(x, each, year, ...) {
  output <- apply(x, 1, function(x) {
    season <- rep(1:round(length(x)/each, 0), each = each, length.out = length(x))
    sighting.season <- aggregate(x, by = list(season), FUN = function(x) sum(x, na.rm = T))$x
    SSR <- length(sighting.season[sighting.season > 0])/length(sighting.season)
    years <- rep(1:round(length(x)/year, 0), each = year, length.out = length(x))
    sighting.year <- aggregate(x, by = list(years), FUN = function(x) sum(x, na.rm = T))$x
    YSR <- length(sighting.year[sighting.year > 0])/length(sighting.year)
    if (SSR >= 0.5 & YSR >= 0.5) {
      category <- "FRS"
    }
    if (SSR < 0.5 & YSR < 0.5 & SSR >= 0.25 & YSR >= 0.25) {
      category <- "FRV"
    }
    if (SSR < 0.5 & YSR >= 0.25) {
      category <- "FRV"
    }
    if (SSR < 0.25 & YSR > 0.25) {
      category <- "OCV"
    }
    if (SSR < 0.25 & YSR < 0.25) {
      category <- "SPV"
    }
    out <- cbind(SSR, YSR, category)
    return(out)
  })
  output <- data.frame(SSR = output[1, ], YSR = output[2, ], category = output[3, ])
  return(output)
}

# Index based on Quintana-Rizzo & Wells (2001), and  Culloch (2004)

index_Quintana <- function(x, each = 3, year, ...) {
  output1 <- apply(x, 1, function(x) {
    years <- rep(1:round(length(x)/year, 0), each = year, length.out = length(x))
    sighting.year <- aggregate(x, by = list(years), FUN = function(x) sum(x, na.rm = T))$x
    category <- rep(0, length(sighting.year))
    category[sighting.year >= 8] <- "Common"
    category[sighting.year > 5 & sighting.year < 8] <- "Frequent"
    category[sighting.year > 2 & sighting.year < 6] <- "Occasional"
    category[sighting.year < 3] <- "Rare"
    return(category)
  })
  output2 <- apply(t(output1), 1, function(x) {
    tab1 <- table(x)
    categories <- c("Rare", "Occasional", "Frequent", "Common")
    tab1 <- tab1[order(match(names(tab1), categories))]
    category <- names(which.max(tab1))
    return(category)
  })
  output <- data.frame(t(output1), General = output2)
  return(output)
}

# Index based on Balmer (2008)

index_Balmer <- function(x, ...) {
  sightings <- apply(x, 1, function(x) sum(x, na.rm = T))
  BINSIZE <- (2 * IQR(apply(x, 1, function(x) sum(x, na.rm = T))))/(nrow(x)^(1/3))
  if (BINSIZE > 0) {
    limits <- seq(from = min(sightings), to = max(sightings) + BINSIZE, by = BINSIZE)
    category <- sightings
    for (i in 1:(length(limits) - 1)) {
      category[category >= limits[i] & category < limits[i + 1]] <- i
    }
  } else {
    limits <- 0
    category <- rep(1, length(sightings))
  }
  output <- list(sightings = sightings, category = category, limits = limits, BINSIZE = BINSIZE)
  return(output)
}

# Index based on Tschopp et al. (2018)

index_Tschopp <- function(x, ...) {
  T.occasions <- ncol(x)
  F.time <- T.occasions - 1
  IO <- (Occurrence(x) - 1)/(T.occasions - 1)
  IT <- Permanency(x)/F.time
  It <- Periodicity(x)
  IA1 <- 1/3 * (IO + IT + It)
  IA2 <- 1/2 * (IO + IT)
  IA3 <- 1/2 * (IO + It)
  IA4 <- 1/2 * (IT + It)
  IH1 <- 3/((1/IO) + (1/It) + (1/IT))
  IH2 <- 2/((1/IO) + (1/IT))
  IH3 <- 2/((1/IO) + (1/It))
  IH4 <- 2/((1/It) + (1/IT))
  output <- data.frame(IA1, IA2, IA3, IA4, IH1, IH2, IH3, IH4)
  return(output)
}

# Index based on Ananias et al. (2008)

index_Ananias <- function(x, year, ...) {
  output1 <- apply(x, 1, function(x) {
    years <- rep(1:round(length(x)/year, 0), each = year, length.out = length(x))
    sighting.year <- aggregate(x, by = list(years), FUN = function(x) sum(x, na.rm = T))$x
    category <- rep(0, length(sighting.year))
    category[sighting.year > 5] <- "Year-round resident"
    category[sighting.year <= 5 & sighting.year >= 3] <- "Seasonal resident"
    category[sighting.year < 3] <- "Transient"
    return(category)
  })
  output2 <- apply(t(output1), 1, function(x) {
    tab1 <- table(x)
    categories <- c("Transient", "Seasonal resident", "Year-round resident")
    tab1 <- tab1[order(match(names(tab1), categories))]
    category <- names(which.max(tab1))
    return(category)
  })
  output <- data.frame(t(output1), General = output2)
  return(output)
}

# Index based on Chabanne et al. (2012)

index_Chabanne <- function(x, each, ...) {
  output <- apply(x, 1, function(x) {
    MSR <- sum(x, na.rm = T)/length(x)
    season <- rep(1:round(length(x)/each, 0), each = each, length.out = length(x))
    sighting.season <- aggregate(x, by = list(season), FUN = function(x) sum(x, na.rm = T))$x
    SSR <- length(sighting.season[sighting.season > 0])/length(sighting.season)
    if (MSR >= 0.1 & SSR >= 0.75) {
      category <- "Resident"
    }
    if (MSR >= 0.1 & SSR < 0.75) {
      category <- "OCV"
    }
    if (MSR < 0.1 & SSR < 0.75 & SSR >= 0.125) {
      category <- "OCV"
    }
    if (MSR < 0.1 & SSR < 0.125) {
      category <- "TRS"
    }
    out <- cbind(MSR = MSR, SSR = SSR, category = category)
    return(out)
  })
  output <- data.frame(MSR = output[1, ], SSR = output[2, ], category = output[3, ])
  return(output)
}

# Index based on Conway (2017)

index_Conway <- function(x, each = 3, year, ...) {
  output1 <- apply(x, 1, function(x) {
    years <- rep(1:round(length(x)/year, 0), each = year, length.out = length(x))
    sighting.year <- aggregate(x, by = list(years), FUN = function(x) sum(x, na.rm = T))$x
    n.sighting.year <- aggregate(x, by = list(years), FUN = length)$x
    p.sighting.year <- sighting.year/n.sighting.year
    season <- rep(1:round(length(x)/each, 0), each = each, length.out = length(x))
    sighting.season <- aggregate(x, by = list(season), FUN = function(x) sum(x, na.rm = T))$x
    sighting.season[sighting.season > 0] <- 1
    season.id <- rep(1:(length(x)/year), each = year/each)
    sighting.season.year <- aggregate(sighting.season, by = list(season.id), FUN = function(x) sum(x, na.rm = T))$x
    p.sighting.season <- sighting.season.year/(year/each)
    category <- rep("Non resident", length(x)/year)
    category[p.sighting.year >= 0.5] <- "Resident"
    category[p.sighting.season > 0.5] <- "Resident"
    return(category)
  })
  output2 <- apply(t(output1), 1, function(x) {
    tab1 <- table(x)
    categories <- c("Non resident", "Resident")
    tab1 <- tab1[order(match(names(tab1), categories))]
    category <- names(which.max(tab1))
    return(category)
  })
  output <- data.frame(t(output1), General = output2)
  return(output)
}

# Index based on Dinis et al. (2016)

index_Dinis <- function(x, each = 3, year, ...) {
  output <- apply(x, 1, function(x) {
    season <- rep(1:round(length(x)/each, 0), each = each, length.out = length(x))
    sighting.season <- aggregate(x, by = list(season), FUN = sum)$x
    sighting.season[sighting.season > 0] <- 1
    year.id <- rep(1:round(length(x)/year), each = (year/each), length.out = length(x)/each)
    sighting.year <- aggregate(sighting.season, by = list(year.id), FUN = function(x) sum(x, na.rm = T))$x
    category <- "Non resident"
    for (i in 1:(length(sighting.year) - 2)) {
      if (sighting.year[i] > 2 & sighting.year[i + 1] > 2 & sighting.year[i + 2] > 2) {
        category <- "Resident"
      }
    }
    return(category)
  })
  return(output)
}

# Index based on Martin & da-Silva (2004)

index_Martin <- function(x, year, ...) {
  output <- apply(x, 1, function(x) {
    years <- rep(1:round(length(x)/year, 0), each = year, length.out = length(x))
    sighting.year <- aggregate(x, by = list(years), FUN = function(x) sum(x, na.rm = T))$x
    if (any(sighting.year >= 7) == TRUE) {
      category <- "Partial resident"
    }
    if (all(sighting.year >= 7) == TRUE) {
      category <- "Permanent resident"
    }
    if (all(sighting.year < 7) == TRUE) {
      category <- "Non resident"
    }
    out <- category
    return(out)
  })
  return(output)
}

# Index based on Möller et al. (2002)

index_Moller <- function(x, each, ...) {
  output <- apply(x, 1, function(x) {
    prop <- sum(x, na.rm = T)/length(x)
    if (prop < 0.1) {
      sighting.rate <- "LSR"
    }
    if (prop >= 0.1 & prop <= 0.3) {
      sighting.rate <- "MSR"
    }
    if (prop > 0.3) {
      sighting.rate <- "HSR"
    }
    season <- rep(1:round(length(x)/each, 0), each = each, length.out = length(x))
    sighting.season <- aggregate(x, by = list(season), FUN = function(x) sum(x, na.rm = T))$x
    length.sighting <- sighting.season[sighting.season > 0]
    category <- NULL
    if ((sighting.rate == "HSR" | sighting.rate == "MSR") & length(length.sighting) > 1) {
      category <- "Resident"
    }
    if ((sighting.rate == "HSR" | sighting.rate == "MSR") & length(length.sighting) == 1) {
      category <- "TRS"
    }
    if (sighting.rate == "LSR" & length(length.sighting) == 1) {
      category <- "TRS"
    }
    if (sighting.rate == "LSR" & length(length.sighting) > 1) {
      category <- "OCV"
    }
    out <- cbind(sighting.rate = sighting.rate, category = category)
    return(out)
  })
  output <- data.frame(Sighting.rate = output[1, ], Category = output[2, ])
  return(output)
}

# Index based on Rosel et al. (2011)

index_Rosel <- function(x, year, ...) {
  output1 <- apply(x, 1, function(x) {
    years <- rep(1:round(length(x)/year, 0), each = year, length.out = length(x))
    sighting.year <- aggregate(x, by = list(years), FUN = function(x) sum(x, na.rm = T))$x
    n.sighting.year <- aggregate(x, by = list(years), FUN = length)$x
    p.sighting.year <- sighting.year/n.sighting.year
    out <- rep("Non resident", length(sighting.year))
    out[p.sighting.year > 0.5] <- "Resident"
    return(out)
  })
  output2 <- apply(t(output1), 1, function(x) {
    tab1 <- table(x)
    categories <- c("Non resident", "Resident")
    tab1 <- tab1[order(match(names(tab1), categories))]
    category <- names(which.max(tab1))
    return(category)
  })
  output <- data.frame(t(output1), General = output2)
  return(output)
}

# Index based on Zolman (2002)

index_Zolman <- function(x, each = 3, year, ...) {
  output <- apply(x, 1, function(x) {
    year.id <- rep(1:round(length(x)/year), each = (year/each), length.out = (year/each) * each)
    season <- rep(1:round(length(x)/each, 0), each = each, length.out = length(x))
    sighting.season <- aggregate(x, by = list(season), FUN = function(x) sum(x, na.rm = T))$x
    season.id <- rep(1:(year/each), times = round(length(x)/year), length.out = (year/each) * each)
    category <- NA
    for (i in 1:(length(x)/year)) {
      sight.year <- sighting.season[year.id == i]
      if (all(sight.year > 0) == TRUE) {
        category <- "Resident"
      }
    }
    for (i in 1:round(year/each)) {
      sight.season <- sighting.season[season.id == i]
      if (all(sight.season > 0) == TRUE & all(sighting.season[season.id > i] == 0) == TRUE) {
        category <- "Seasonal resident"
      }
    }
    if (is.na(category) == TRUE) {
      category <- "Transient"
    }
    return(category)
  })
  return(output)
}

# Bibliography

# Ballance, L. T. (1990). Residency patterns, group organization, and surfacing associations of bottlenose dolphins 
# in Kino Bay, Gulf of California, Mexico. En The Bottlenose Dolphin (pp. 267–283). Academic Press. 
# https://doi.org/10.1016/B978-0-12-440280-5.50017-2
#
# Parra, G. J., Corkeron, P. J., & Marsh, H. (2006). Population sizes, site fidelity and residency patterns of 
# Australian snubfin and Indo-Pacific humpback dolphins: Implications for conservation. Biological Conservation, 
# 129(2), 167–180. https://doi.org/10.1016/j.biocon.2005.10.031
#
# Díaz-López, B. (2012). Bottlenose dolphins and aquaculture: Interaction and site fidelity on the North-Eastern 
# Coast of Sardinia (Italy). Marine Biology, 159(10), 2161–2172. https://doi.org/10.1007/s00227-012-2002-x
#
# Quintana-Rizzo, E., & Wells, R. S. (2001). Resighting and association patterns of bottlenose dolphins (Tursiops 
# truncatus) in the Cedar Keys, Florida: Insights into social organization. Canadian Journal of Zoology, 79(3), 
# 447–456. https://doi.org/10.1139/cjz-79-3-447
#
# Balmer, B. C., Wells, R. S., Nowacek, S. M., Nowacek, D. P., Schwacke, L. H., McLellan, W. A., Scharf, F. S., 
# Rowles, T. K., Hansen, L. J., Spradlin, T. R., & Pabst, D. A. (2008). Seasonal abundance and distribution patterns 
# of common bottlenose dolphins (Tursiops truncatus) near St. Joseph Bay, Florida, USA. Cetacean Res. Manage, 10(2), 
# 157–167.
#
# Tschopp, A., Ferrari, M. A., Crespo, E. A., & Coscarella, M. A. (2018). Development of a site fidelity index based 
# on population capture-recapture data. PeerJ, 6, e4782. https://doi.org/10.7717/peerj.4782
#
# Ananias, S. M. A., Jesus, A. H., & Yamamoto, M. E. (2008). Recorrência e fidelidade espacial do boto-cinza Sotalia 
# guianensis na enseada do Curral, Pipa/RN, avaliada através da fotoidentificação. En A. H. Jesus, P. I. A. P. 
# Medeiros, & F. J. L. Silva, Boto-cinza Sotalia guianensis: Pesquisa e conservação no nordeste do Brasil (pp. 
# 61–77). Edições UERN.
#
# Chabanne, D., Finn, H., Salgado-Kent, C., & Bedjer, L. (2012). Identification of a resident community of bottlenose 
# dolphins (Tursiops aduncus) in the SwanCanning Riverpark, Western Australia, using behavioural information. Pacific 
# Conservation Biology, 18(4), 247–262. https://doi.org/10.1071/PC120247
#
# Conway, J. N. (2017). Estimating density and residency of bottlenose dolphins (Tursiops truncatus) in three estuarine 
# sites in South Carolina [Electronic Theses and Dissertations]. Coastal Carolina University.
#
# Dinis, A., Alves, F., Nicolau, C., Ribeiro, C., Kaufmann, M., Cañadas, A., & Freitas, L. (2016). Bottlenose dolphin 
# Tursiops truncatus group dynamics, site fidelity, residency and movement patterns in the Madeira Archipelago (North-
# East Atlantic). African Journal of Marine Science, 38(2), 151–160. https://doi.org/10.2989/1814232X.2016.1167780
#
# Martin, A. R., & Silva, V. M. F. (2004). Number, seasonal movements, and residency characteristics of river dolphins 
# in an Amazonian floodplain lake system. Canadian Journal of Zoology, 82(8), 1307–1315. https://doi.org/10.1139/z04-109
#
# Möller, L. M., Allen, S. J., & Harcourt, R. G. (2002). Group characteristics, site fidelity and seasonal abundance of 
# bottlenosed dolphins (Tursiops aduncus) in Jervis Bay and Port Stephens, South-Eastern Australia. Australian Mammalogy,
# 24, 11–21.3
#
# Rosel, P. E., Mullin, K. D., Garrison, L., Schwacke, L., Adams, J., Balmer, B., Conn, P., Conroy, M. J., Eguchi, T., 
# Gorgone, A., Hohn, A., Mazzoil, M., Schwarz, C., Sinclair, C., Speakman, T., Urian, K., Vollmer, N., Wade, P., Wells, 
# R., & Zolman, E. (2011). Photo-identification capture-mark-recapture techniques for estimating abundance of bay, sound 
# and estuary populations of bottlenose dolphins along the U.S. East Coast and Gulf of Mexico: A workshop report. NOAA 
# Technical Memorandum NMFS-SEFSC-621.
#
# Zolman, E. S. (2002). Residency patterns of bottlenose dolphins (Tursiops truncatus) in the Stono River Estuary, 
# Charleston County, South Carolina, U.S.A. Marine Mammal Science, 18(4), 879–892. 
# https://doi.org/10.1111/j.1748-7692.2002.tb01079.x

