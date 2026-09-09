
##################### PAPER 3 #############################
#################### Turkey vultures - RASO ######################
#################### Spatial Network ######################

# Pos doc project - Spatial Networks
# Date: 11/08/26
# Author: Daniela P. Coelho
# Aim: Performing movement analysis and plotting figures


################ RUN "0.Base_code.R" #########################


# Checking the distribution of locations among months by seasons
(m_dry <- as.data.frame(table(dry_season$month)))
(m_dry1 <- m_dry %>% rename(month = Var1))


# Number of locations by month
(g1 <- ggplot(m_dry1, aes(x = month, y = Freq)) +
    geom_bar(stat = "identity", fill = "orange") +
    labs(title = "",
         x = "Month", y = "") +
    ylim(0,16000) +
    geom_text(aes(label = Freq), vjust = -0.5, size = 2.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1)) +
    theme(plot.title = element_text(hjust=0.5)))


m_rainy <- table(rainy_season$month)
(m_rainy <- as.data.frame(table(rainy_season$month)))
(m_rainy1 <- m_rainy %>% rename(month = Var1))

(g2 <- ggplot(m_rainy1, aes(x = month, y = Freq)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "", 
         x = "Month", y = "Number of locations") +
    geom_text(aes(label = Freq), vjust = -0.5, size = 2.5) +
    ylim(0,16000) +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1)) +
    theme(plot.title = element_text(hjust=0.5)))


#### Organize multiple charts into a single figure
library(ggpubr)

ggarrange(g2, g1, nrow =1, ncol = 2, legend = "right", labels = c("a", "b"),  common.legend = TRUE) 


# Counting the number of points per individual for the dry season
dry_counts <- dry_season %>%
  count(id, name = "Fixes") %>%
  arrange(desc(Fixes))


(g3 <- ggplot(dry_counts, aes(x = reorder(as.character(id), -Fixes), y = Fixes)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(title = "", x = "Individual ID", y = "") +
    geom_text(aes(label = Fixes), vjust = -0.5, size = 2.5) +
    ylim(0,6500) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust=0.5)))

write.csv2(dry_counts, "path/dry_individual_nlocs_turkey.csv", row.names = F)

# Rainy season
rainy_counts <- rainy_season %>%
  count(id, name = "Fixes") %>%
  arrange(desc(Fixes))

(g4 <- ggplot(rainy_counts, aes(x = reorder(as.character(id), -Fixes), y = Fixes)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "", x = "Individual ID", y = "Number of locations") +
    geom_text(aes(label = Fixes), vjust = -0.5, size = 2.5) +
    ylim(0, 6500) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust=0.5)))

write.csv2(rainy_counts, "path/rainy_individual_nlocs_turkey.csv", row.names = F)

ggarrange(g4, g3, nrow =1, ncol = 2, legend = "right", labels = c("a", "b"),  common.legend = TRUE) 


## Combining the 4 figures into a single panel
ggarrange(g2, g1, g4, g3, nrow =2, ncol = 2, legend = "right", labels = c("a", "b", "c", "d"),  common.legend = TRUE) 



### Analyzing the distances traveled by individuals at each station
library(geosphere) # Calculating step-lengths

# Auxiliary function to calculate step lengths per individual
calc_distances <- function(df) {
  df <- df %>%
    arrange(id, time) %>%
    group_by(id) %>%
    mutate(
      lead_lon = lead(location_long),
      lead_lat = lead(location_lat),
      step_length = distGeo(matrix(c(location_long, location_lat), ncol = 2),matrix(c(lead_lon, lead_lat), ncol = 2))
    ) %>%
    ungroup()
  return(df)
}

# Applying to the seasons
dry_season_dist <- calc_distances(dry_season)
rainy_season_dist <- calc_distances(rainy_season)

# Abstract function
resumo_dist <- function(df) {
  df %>%
    group_by(id) %>%
    summarise(
      n_pontos = n(),
      min_dist = min(step_length, na.rm = TRUE),
      max_dist = max(step_length, na.rm = TRUE),
      mean_dist = mean(step_length, na.rm = TRUE),
      median_dist = median(step_length, na.rm = TRUE)
    ) %>%
    arrange(desc(mean_dist))
}

# Applying to the seasons
dry_summary <- resumo_dist(dry_season_dist)
rainy_summary <- resumo_dist(rainy_season_dist)

write.csv2(dry_summary, "path/dry_individual_step_lenght_turkey.csv", row.names = F)

write.csv2(rainy_summary, "path/rainy_individual_step_lenght_turkey.csv", row.names = F)


# Maximum distances per id
(g6 <- ggplot(dry_summary, aes(x = reorder(as.character(id), -max_dist), y = max_dist)) +
          geom_bar(stat = "identity", fill = "orange") +
          labs(title = "", 
               x = "Individual ID", y = "") +
          ylim(0, 80000) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(plot.title = element_text(hjust=0.5)))

(g5 <- ggplot(rainy_summary, aes(x = reorder(as.character(id), -max_dist), y = max_dist)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "", x = "Individual ID", y = "Maximum step-length (m)") +
    ylim(0, 136000) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust=0.5)))


ggarrange(g5, g6, nrow =1, ncol = 2, legend = "right", labels = c("a", "b"),  common.legend = TRUE) 


# Step length
(g7 <- ggplot(dry_season_dist, aes(x = as.character(id), y = step_length)) +
    geom_boxplot(outlier.colour = "orange", outlier.shape = 16) +
    labs(title = "", x = "Individual ID", y = "") +
    #ylim(0, 200000) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust=0.5)))


(g8 <- ggplot(rainy_season_dist, aes(x = as.character(id), y = step_length)) +
    geom_boxplot(outlier.colour = "skyblue", outlier.shape = 16) +
    labs(title = "", x = "Individual ID", y = "Step-length (m)") +
    #ylim(0, 50000) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust=0.5)))


ggarrange(g8, g7, nrow =1, ncol = 2, legend = "right", labels = c("a", "b"),  common.legend = TRUE) 



gc() #garbage collection























#############################################################################################################################################
##### PLOT locations per individual

# Using the `transform_coord_rainy` object from '0.Base_code.R'
transform_coord_rainy <- st_coordinates(df_sf_all)
transform_coord_dry <- st_coordinates(df_sf_all)

# Transforming IDs into characters - used to color individuals in the plot
df_sf_rainy1 <- as.character(df_sf_all$id)
df_sf_dry1 <- as.character(df_sf_all$id)

# Tracking plots (30') - UTM - RAINY
(g9 <- ggplot(transform_coord_rainy, aes(x = X, y = Y, color = df_sf_rainy1)) + 
    geom_point(size = 0.8) +  
    #scale_color_manual(values=c("5568" = "#dd3497", "5570" = "#4daf4a", "9025" = "#ff7f00", "8956" = "#386cb0", "9026" = "#252525", "9250" = "#c994c7")) +
    labs(title = "", 
                x = "Longitude", 
                y = "Latitude",
                colour = "ID Individual") + 
    #ylim(8700000, 9100000) + 
    #xlim(440000, 660000) +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5)))


# Tracking plots (30') - UTM - DRY
(g10 <- ggplot(transform_coord_dry, aes(x = X, y = Y, color = df_sf_dry1)) + 
    geom_point(size = 0.8) +
    #scale_color_manual(values=c("5568" = "#dd3497", "5570" = "#4daf4a", "9025" = "#ff7f00", "9024" = "#00441b", "5574" = "#4d004b", "5573" = "#e31a1c", "9252" = "#0570b0")) +
    labs(title = "", 
         x = "Longitude", 
         y = "",
         colour = "ID Individual") + 
    #ylim(8700000, 9100000) + 
    #xlim(440000, 660000) +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5)))

ggarrange(g9, g10, nrow =1, ncol = 2, legend = "bottom", labels = c("a", "b")) 



############################################################
####################### NETWORK ANALYSIS ###################
############################################################

##### Using moveNT package 

# Using the `ltraj_df_dry and rainy` from '0.Base_code.R'

summary_dry <- summary(ltraj_df_dry) # Summary of trajectories
df_summary_dry <- as.data.frame(summary_dry)

write.csv2(df_summary_dry, "C:/Users/danip/Desktop/Pos doc/Projeto/Results_macaw/coord_interval_30'/Raso/Seasonality/Four/Dry_summary.csv", row.names = F)

summary_rainy <- summary(ltraj_df_rainy)
df_summary_rainy <- as.data.frame(summary_rainy)

write.csv2(df_summary_rainy, "C:/Users/danip/Desktop/Pos doc/Projeto/Results_macaw/coord_interval_30'/Raso/Seasonality/Four/Rainy_summary.csv", row.names = F)


# Plotting the movement route in ltraj with coordinates in meters
plot(ltraj_df_dry, xlab='Longitude', ylab='Latitude')
plot(ltraj_df_rainy, xlab='Longitude', ylab='Latitude')

par(mar=c(2,2,2,2))
par(mfrow = c(2,2))
plot(ltraj_df_dry[1], xlab='Longitude', ylab='Latitude', main = "ID 5568 - Dry")
plot(ltraj_df_dry[2], xlab='Longitude', ylab='', main = "ID 5570 - Dry")
plot(ltraj_df_rainy[1], xlab='Longitude', ylab='Latitude', main = "ID 5568 - Rainy")
plot(ltraj_df_rainy[2], xlab='Longitude', ylab='', main = "ID 5570 - Rainy")

plot(ltraj_df_dry[7], xlab='Longitude', ylab='Latitude', main = "ID 9025 - Dry")
plot(ltraj_df_dry[5], xlab='Longitude', ylab='Latitude', main = "ID  9018 - Dry")
plot(ltraj_df_rainy[5], xlab='Longitude', ylab='Latitude', main = "ID 9025 - Rainy")
plot(ltraj_df_rainy[4], xlab='Longitude', ylab='Latitude', main = "ID 9018 - Rainy")

plot(ltraj_df_dry[3], xlab='Longitude', ylab='Latitude', main = "ID 5573 - Dry")
plot(ltraj_df_dry[4], xlab='Longitude', ylab='Latitude', main = "ID 5574 - Dry")
plot(ltraj_df_dry[6], xlab='Longitude', ylab='Latitude', main = "ID 9024 - Dry")
plot(ltraj_df_dry[8], xlab='Longitude', ylab='Latitude', main = "ID 9252 - Dry")


plot(ltraj_df_rainy[3], xlab='Longitude', ylab='Latitude', main = "ID 8956 - Rainy")
plot(ltraj_df_rainy[6], xlab='Longitude', ylab='Latitude', main = "ID 9026 - Rainy")
plot(ltraj_df_rainy[7], xlab='Longitude', ylab='Latitude', main = "ID 9250 - Rainy")

dev.off()

gc() #garbage collection
