
##################### PAPER 3 #############################
#################### Turkey vultures - RASO ######################
#################### Spatial Network ######################

# Pos doc project - Spatial Networks
# Date: 11/08/26
# Author: Daniela P. Coelho
# Aim: Base code for performing the Paper 3 analyses - spatial networks of Turkey vultures from Raso da Catarina.

## defining the working directory
setwd("path")

# used packages to import movebank data - online
library(RCurl) # Data download directly from the internet
library(bitops) # RCurl extension
library(move) # Movement data / to work with class move object

# Login to movebank
loginStored <- movebankLogin(username="username", password="password")

# Informations about the study
getMovebankStudy(study="Cathartid ecology - Dénes - Brazil", login=loginStored)

# Getting data from movebank for the group of individuals from Raso and removing duplicated timestamps (n = 10 individuals)
LM <- getMovebankData(study="Cathartid ecology - Dénes - Brazil", animalName=c("8833","8835","8838","8841","8842","8848","8849","8850","8851","8852","8853","8907","8908","8909","8910","8911","8912","8913","8914","8916","8917","8918","8919","8920","8940"), sensorID="GPS", login=loginStored, removeDuplicatedTimestamps=TRUE)

# Packages for manipulating movement data
library(move2) # New version of the move package / Improved integration with the sf package
library(sf) # Movement data / maps trajectories and corrects coordinate system (CRS)
library(lubridate) # Working with dates
library(dplyr) # Working with data and time / 
library(ggplot2)

############################### STEP 1 ####################################
####################### Manipulating the dataset ##########################

# Transforming the move object data into a dataframe
df <- as.data.frame(LM)

# Selecting the specific columns
df <- df %>%
  dplyr::select(location_lat, location_long, timestamp, tag_local_identifier)

# Creating a move2 object (necessary to apply my_filter_per_interval function)
df <- df %>%
  rename(x = location_long, y = location_lat, time = timestamp, id =  tag_local_identifier)

move2_df <- mt_as_move2(df, coords = c("x", "y"), time_column = "time", track_id_column = "id", crs = 4326)

# Filtering 30' intervals by coordinate and selecting the columns of interest
move2_df <- mt_filter_per_interval(move2_df, unit = "30 minute", criterion = "first") # hour or minute 

# Extracting the coordinates (lat/long) from the geometry column of the move2_df Transformig the move2 object into a dataframe. 

coordinates_ll <- st_coordinates(move2_df)

df_fit <- move2_df %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  mutate(
    location_long = coordinates_ll[, 1],
    location_lat  = coordinates_ll[, 2],
    
    # Keeping times in UTC to calculate time intervals
    time_utc = as.POSIXct(time, tz = "UTC"),
    
    # Local time for classifying month and season
    time_local = with_tz(time_utc, tzone = "America/Bahia"),
    
    month = month(time_local),
    
    season = case_when(
      month %in% 6:12 ~ "Dry",
      month %in% 1:5  ~ "Rainy"
    )
  ) %>%
  arrange(id, time_utc)


gc() # Garbage collection (Identifies objects that are no longer in use; Frees up the memory they occupy; Returns a memory usage report)


############################### STEP 2 ####################################
########## Transforming spatial data into an adjacency matrix #############
# 1 →  Converting the geographic projection (CRS) from lat/long (degrees) to UTM (meters);
# 2 → Create an object of the ltraj class needed to transform the movement data into an adjacency matrix;
# 3 → Creating a spatial grid that will be superimposed on the movement route where each cell of the grid will be equivalent to a node of the spatial network with a size of 2,000 m on each side = 4km2 in area.
# 4 → Convert the movement routes into the adjacency matrix.


## Converting the dataframe (df_fit) into a ltraj class (necessary to converts a trajectory object into an adjacency matrix)
library(adehabitatLT) # working with movement paths / Creates and manipulates ltraj (trajectory) objects. The function as.ltraj() transforms raw data into a structured trajectory; However, the problem is that ltraj objects do not store CRS (coordinate system).
library(terra) # For converting objects of the sf class into raster objects
library(moveNT) # Transforming trajectories into spatial networks (movement networks)

# Projecting all data to UTM 24S

df_sf_all <- st_as_sf(
  df_fit,
  coords = c("location_long", "location_lat"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(32724)

coordinates_utm <- st_coordinates(df_sf_all)

df_utm <- df_sf_all %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  mutate(
    x_utm = coordinates_utm[, 1],
    y_utm = coordinates_utm[, 2]
  ) %>%
  arrange(id, time_utc)


# Creating a single 2×2 km base raster

grid_spacing <- 2000

# The boundaries are aligned to exact multiples of 2,000 m. An additional cell is added in each direction to ensure the raster is larger than the extent of the points.

xmin <- floor(min(df_utm$x_utm, na.rm = TRUE) /
                grid_spacing) * grid_spacing - grid_spacing

xmax <- ceiling(max(df_utm$x_utm, na.rm = TRUE) /
                  grid_spacing) * grid_spacing + grid_spacing

ymin <- floor(min(df_utm$y_utm, na.rm = TRUE) /
                grid_spacing) * grid_spacing - grid_spacing

ymax <- ceiling(max(df_utm$y_utm, na.rm = TRUE) /
                  grid_spacing) * grid_spacing + grid_spacing


# Exact number of rows and columns
n_cols <- as.integer((xmax - xmin) / grid_spacing)
n_rows <- as.integer((ymax - ymin) / grid_spacing)


# Raster base
grid_raster <- raster::raster(
  xmn = xmin,
  xmx = xmax,
  ymn = ymin,
  ymx = ymax,
  ncols = n_cols,
  nrows = n_rows
)

# Assigning the CRS using WKT
raster::crs(grid_raster) <- sf::st_crs(32724)$wkt

# ach cell receives its ID in the raster.
raster::values(grid_raster) <-
  seq_len(raster::ncell(grid_raster))


# Creating the vector version of the same grid

grid_bbox <- st_bbox(
  c(
    xmin = xmin,
    ymin = ymin,
    xmax = xmax,
    ymax = ymax
  ),
  crs = st_crs(32724)
)

# grid_bbox
#  xmin    ymin     xmax    ymax 
# 432000  8752000  664000  9046000 


grid_geometry <- st_make_grid(
  st_as_sfc(grid_bbox),
  cellsize = grid_spacing,
  offset = c(xmin, ymin),
  square = TRUE
)

grid_sf <- st_sf(geometry = grid_geometry)


# The cell number is obtained from the raster, not from the order returned by st_make_grid().

grid_centroids <- st_coordinates(st_centroid(grid_sf))

grid_sf$cell_id <- raster::cellFromXY(
  grid_raster,
  grid_centroids
)

grid_sf <- grid_sf %>%
  arrange(cell_id)

# Checking: a geometry for each raster cell
stopifnot(nrow(grid_sf) == raster::ncell(grid_raster))


# Assigning each ID to a permanent grid cell

df_utm$cell_id <- raster::cellFromXY(
  grid_raster,
  as.matrix(df_utm[, c("x_utm", "y_utm")])
)

# No point should be outside the raster
stopifnot(!any(is.na(df_utm$cell_id)))


# Split the seasons

dry_season <- df_utm %>%
  filter(season == "Dry") %>%
  arrange(id, time_utc)

rainy_season <- df_utm %>%
  filter(season == "Rainy") %>%
  arrange(id, time_utc)


# Defining independent time segments: A new burst begins when the interval exceeds 300 minutes (5h). This value should be verified against the distribution of actual intervals. For data aggregated into 30-minute intervals, a 300-minute threshold allows for up to two missing windows without linking widely separated periods.

max_gap_minutes <- 300


# DRY SEASON

dry_season <- dry_season %>%
  group_by(id) %>%
  arrange(time_utc, .by_group = TRUE) %>%
  mutate(
    dt_minutes = as.numeric(
      difftime(
        time_utc,
        lag(time_utc),
        units = "mins"
      )
    ),
    
    new_burst = is.na(dt_minutes) |
      dt_minutes > max_gap_minutes,
    
    burst_number = cumsum(new_burst),
    
    burst_id = paste(
      id,
      "Dry",
      sprintf("%03d", burst_number),
      sep = "_"
    )
  ) %>%
  ungroup()



# RAINY SEASON

rainy_season <- rainy_season %>%
  group_by(id) %>%
  arrange(time_utc, .by_group = TRUE) %>%
  mutate(
    dt_minutes = as.numeric(
      difftime(
        time_utc,
        lag(time_utc),
        units = "mins"
      )
    ),
    
    new_burst = is.na(dt_minutes) |
      dt_minutes > max_gap_minutes,
    
    burst_number = cumsum(new_burst),
    
    burst_id = paste(
      id,
      "Rainy",
      sprintf("%03d", burst_number),
      sep = "_"
    )
  ) %>%
  ungroup()


# Creating ltraj objects

ltraj_df_dry <- as.ltraj(
  xy = as.matrix(dry_season[, c("x_utm", "y_utm")]),
  date = dry_season$time_utc,
  id = dry_season$id,
  burst = dry_season$burst_id,
  typeII = TRUE,
  proj4string = sp::CRS(SRS_string = "EPSG:32724")
)

ltraj_df_rainy <- as.ltraj(
  xy = as.matrix(rainy_season[, c("x_utm", "y_utm")]),
  date = rainy_season$time_utc,
  id = rainy_season$id,
  burst = rainy_season$burst_id,
  typeII = TRUE,
  proj4string = sp::CRS(SRS_string = "EPSG:32724")
)


# Creating dry-season transitions
# cell_from is the cell occupied at time t.
# cell_to is the cell occupied in the next record for the same individual and the same burst.

edges_dry <- dry_season %>%
  group_by(id, burst_id) %>%
  arrange(time_utc, .by_group = TRUE) %>%
  mutate(
    cell_to = lead(cell_id)
  ) %>%
  ungroup() %>%
  filter(!is.na(cell_to)) %>%
  transmute(
    cell_from = cell_id,
    cell_to = cell_to
  ) %>%
  
  # Adds equal transitions between all individuals.
  count(
    cell_from,
    cell_to,
    name = "weight"
  )


# Dry season node lookup

used_cells_dry <- sort(unique(dry_season$cell_id))

node_lookup_dry <- data.frame(
  node_id = seq_along(used_cells_dry),
  cell_id = used_cells_dry
)

node_coordinates_dry <- raster::xyFromCell(
  grid_raster,
  node_lookup_dry$cell_id
)

node_lookup_dry <- node_lookup_dry %>%
  mutate(
    x = node_coordinates_dry[, 1],
    y = node_coordinates_dry[, 2],
    season = "Dry"
  )


write.csv2(node_lookup_dry,"path/Nodes_coord_dry_turkey.csv", row.names = FALSE)

# Relacionar cell_id permanente com a posição compacta utilizada na matriz.

edges_dry_nodes <- edges_dry %>%
  left_join(
    node_lookup_dry %>%
      select(cell_id, node_id),
    by = c("cell_from" = "cell_id")
  ) %>%
  rename(node_from = node_id) %>%
  
  left_join(
    node_lookup_dry %>%
      select(cell_id, node_id),
    by = c("cell_to" = "cell_id")
  ) %>%
  rename(node_to = node_id)



# Creating dry season adjacency matrix
mat_dry <- matrix(
  0,
  nrow = nrow(node_lookup_dry),
  ncol = nrow(node_lookup_dry),
  dimnames = list(
    as.character(node_lookup_dry$cell_id),
    as.character(node_lookup_dry$cell_id)
  )
)

mat_dry[
  cbind(
    edges_dry_nodes$node_from,
    edges_dry_nodes$node_to
  )
] <- edges_dry_nodes$weight


# Raster relating matrix position to geographic position
patch_dry <- raster::raster(grid_raster)
raster::values(patch_dry) <- NA_integer_

raster::values(patch_dry)[node_lookup_dry$cell_id] <-
  node_lookup_dry$node_id


# Object with the structure expected by moveNT
adj_dry <- list(
  mat_dry,
  grid_raster,
  patch_dry
)

class(adj_dry) <- "adjmov"


# Creating the rainy season transition 

edges_rainy <- rainy_season %>%
  group_by(id, burst_id) %>%
  arrange(time_utc, .by_group = TRUE) %>%
  mutate(
    cell_to = lead(cell_id)
  ) %>%
  ungroup() %>%
  filter(!is.na(cell_to)) %>%
  transmute(
    cell_from = cell_id,
    cell_to = cell_to
  ) %>%
  count(
    cell_from,
    cell_to,
    name = "weight"
  )


# Rainy season node lookup

used_cells_rainy <- sort(unique(rainy_season$cell_id))

node_lookup_rainy <- data.frame(
  node_id = seq_along(used_cells_rainy),
  cell_id = used_cells_rainy
)

node_coordinates_rainy <- raster::xyFromCell(
  grid_raster,
  node_lookup_rainy$cell_id
)

node_lookup_rainy <- node_lookup_rainy %>%
  mutate(
    x = node_coordinates_rainy[, 1],
    y = node_coordinates_rainy[, 2],
    season = "Rainy"
  )


write.csv2(node_lookup_rainy,"path/0.Nodes_coord_rainy_turkey.csv", row.names = FALSE)


edges_rainy_nodes <- edges_rainy %>%
  left_join(
    node_lookup_rainy %>%
      select(cell_id, node_id),
    by = c("cell_from" = "cell_id")
  ) %>%
  rename(node_from = node_id) %>%
  
  left_join(
    node_lookup_rainy %>%
      select(cell_id, node_id),
    by = c("cell_to" = "cell_id")
  ) %>%
  rename(node_to = node_id)


# Rainy season adjacency matrix

mat_rainy <- matrix(
  0,
  nrow = nrow(node_lookup_rainy),
  ncol = nrow(node_lookup_rainy),
  dimnames = list(
    as.character(node_lookup_rainy$cell_id),
    as.character(node_lookup_rainy$cell_id)
  )
)

mat_rainy[
  cbind(
    edges_rainy_nodes$node_from,
    edges_rainy_nodes$node_to
  )
] <- edges_rainy_nodes$weight


patch_rainy <- raster::raster(grid_raster)
raster::values(patch_rainy) <- NA_integer_

raster::values(patch_rainy)[node_lookup_rainy$cell_id] <-
  node_lookup_rainy$node_id


adj_rainy <- list(
  mat_rainy,
  grid_raster,
  patch_rainy
)

class(adj_rainy) <- "adjmov"



#### Checking

# Do the two objects use exactly the same raster?
raster::compareRaster(
  adj_dry[[2]],
  adj_rainy[[2]],
  extent = TRUE,
  rowcol = TRUE,
  crs = TRUE,
  res = TRUE,
  orig = TRUE,
  stopiffalse = TRUE
)


# Resolution: 2.000 × 2.000 m
raster::res(grid_raster)


# Number of nodes per season
nrow(adj_dry[[1]]) # 532 nodes

nrow(adj_rainy[[1]]) # 1686 nodes

# Total number of valid transitions expected:
# number of records minus number of bursts.

expected_transitions_dry <-
  nrow(dry_season) -
  dplyr::n_distinct(dry_season$burst_id)

expected_transitions_rainy <-
  nrow(rainy_season) -
  dplyr::n_distinct(rainy_season$burst_id)


# The sum of the matrix weights must equal the number of valid transitions.

sum(adj_dry[[1]]) # 13828
expected_transitions_dry # 13828

sum(adj_rainy[[1]]) # 40951
expected_transitions_rainy # 40951

stopifnot(
  sum(adj_dry[[1]]) == expected_transitions_dry
)

stopifnot(
  sum(adj_rainy[[1]]) == expected_transitions_rainy
)


# Cells used in both seasons
common_cells <- intersect(
  node_lookup_dry$cell_id,
  node_lookup_rainy$cell_id
)

length(common_cells) # 399

# Table showing the position occupied by the same cell
# in each of the matrices.

node_correspondence <- node_lookup_dry %>%
  select(
    cell_id,
    node_id_dry = node_id,
    x,
    y
  ) %>%
  full_join(
    node_lookup_rainy %>%
      select(
        cell_id,
        node_id_rainy = node_id
      ),
    by = "cell_id"
  ) %>%
  mutate(
    used_dry = !is.na(node_id_dry),
    used_rainy = !is.na(node_id_rainy)
  )

write.csv2(node_correspondence,"path/0.Nodes_correspondence_turkey.csv", 
           row.names = FALSE)


# Visual checking

points_dry_sf <- st_as_sf(
  dry_season,
  coords = c("x_utm", "y_utm"),
  crs = 32724
)

points_rainy_sf <- st_as_sf(
  rainy_season,
  coords = c("x_utm", "y_utm"),
  crs = 32724
)


ggplot() +
  geom_sf(
    data = grid_sf,
    fill = NA,
    color = "grey70",
    linewidth = 0.1
  ) +
  geom_sf(
    data = points_dry_sf,
    color = "orange",
    size = 0.2
  ) +
  geom_sf(
    data = points_rainy_sf,
    color = "lightblue",
    size = 0.2
  ) +
  coord_sf(datum = NA) +
  theme_bw()



#### OBS:
# cell_id: This is the permanent spatial identifier:
# node_id: This is simply the position of the occupied cell within that station's compact matrix.


gc() #garbage collection



#---------------------------#
### Choosing the burst limit
#---------------------------#

# Checking the distribution of time intervals:

interval_distribution <- df_utm %>%
  
  dplyr::group_by(
    id,
    season
  ) %>%
  
  dplyr::arrange(
    time_utc,
    .by_group = TRUE
  ) %>%
  
  dplyr::mutate(
    dt_minutes = as.numeric(
      difftime(
        time_utc,
        dplyr::lag(time_utc),
        units = "mins"
      )
    )
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::filter(
    !is.na(dt_minutes)
  )


summary(
  interval_distribution$dt_minutes
)


quantile(
  interval_distribution$dt_minutes,
  probs = c(
    0.50,
    0.75,
    0.90,
    0.95,
    0.99
  ),
  na.rm = TRUE
)


# Checking the frequencies:
interval_distribution %>%
  
  dplyr::count(
    dt_minutes
  ) %>%
  
  dplyr::arrange(
    dt_minutes
  )


ggplot(
  interval_distribution,
  aes(x = dt_minutes)
) +
  
  geom_histogram(
    binwidth = 30,
    fill = "grey70",
    color = "black"
  ) +
  
  geom_vline(
    xintercept = 300,
    linetype = "dashed",
    linewidth = 1
  ) +
  
  coord_cartesian(
    xlim = c(0, 600)
  ) +
  
  labs(
    x = "Time interval between consecutive locations (min)",
    y = "Frequency"
  ) +
  
  theme_bw()









