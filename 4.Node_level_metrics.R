

##################### PAPER 3 #############################
#################### Turkey vultures - RASO ######################
#################### Spatial Network ######################

# Pos doc project - Spatial Networks
# Date: 11/08/26
# Author: Daniela P. Coelho
# Aim: Calculating node-level metrics (self loop, degree and betweenness centrality, and global clustering (transitivity), according to seasonality (dry season and rainy season) for Turkey vultures (n = 25) from RASO using node size  = 4 km2 (2000 x 2000 m).


########################## RUN "0.Base_code.R" #########################

library (moveNT) # spatial networks
library(ggplot2) # plots
library(ggpubr) # union of figures
library(terra) # raster manipulation
library(adehabitatLT) # ltraj
library(dplyr) # data manipulation
library(raster)
library(igraph)
library(tidyr)
library(lubridate)


################################################################
# OBTENDO OS RESULTADOS ORIGINAIS CALCULADOS PELO PACOTE MOVE NT
################################################################

stck_dry <- moveNT::adj2stack(
  adj_dry,
  mode = "directed",
  weighted = TRUE
)

stck_rainy <- moveNT::adj2stack(
  adj_rainy,
  mode = "directed",
  weighted = TRUE
)


# Conferir os nomes das camadas
names(stck_dry)
names(stck_rainy)



################################################################
# EXTRAINDO OS RESULTADOS pelo MOVENT PARA CONFERÊNCIA POSTERIOR
################################################################

metrics_moveNT_dry <- node_lookup_dry %>%
  transmute(
    node_id,
    cell_id,
    x,
    y,
    season = "Dry",
    
    weight_moveNT =
      raster::getValues(stck_dry[[2]])[cell_id],
    
    self_loop_moveNT =
      raster::getValues(stck_dry[[3]])[cell_id],
    
    degree_moveNT =
      raster::getValues(stck_dry[[4]])[cell_id],
    
    betweenness_moveNT =
      raster::getValues(stck_dry[[5]])[cell_id],
    
    transitivity_moveNT =
      raster::getValues(stck_dry[[6]])[cell_id]
  )


metrics_moveNT_rainy <- node_lookup_rainy %>%
  transmute(
    node_id,
    cell_id,
    x,
    y,
    season = "Rainy",
    
    weight_moveNT =
      raster::getValues(stck_rainy[[2]])[cell_id],
    
    self_loop_moveNT =
      raster::getValues(stck_rainy[[3]])[cell_id],
    
    degree_moveNT =
      raster::getValues(stck_rainy[[4]])[cell_id],
    
    betweenness_moveNT =
      raster::getValues(stck_rainy[[5]])[cell_id],
    
    transitivity_moveNT =
      raster::getValues(stck_rainy[[6]])[cell_id]
  )


############################################################
# MÉTRICAS — DRY SEASON
############################################################

mat_dry <- adj_dry[[1]]

n_nodes_dry <- nrow(mat_dry)
total_transitions_dry <- sum(mat_dry)


# ----------------------------------------------------------
# Matriz sem self-loops
# ----------------------------------------------------------

mat_dry_no_loop <- mat_dry

diag(mat_dry_no_loop) <- 0


# ----------------------------------------------------------
# Matriz binária
#
# Uma ligação recebe 1 independentemente da frequência.
# ----------------------------------------------------------

mat_dry_binary <- 1L * (mat_dry_no_loop > 0)


# ----------------------------------------------------------
# Grafo binário direcionado
# ----------------------------------------------------------

g_dry_binary <- igraph::graph_from_adjacency_matrix(
  mat_dry_binary,
  mode = "directed",
  weighted = NULL,
  diag = FALSE
)


# ----------------------------------------------------------
# Grafo ponderado pela frequência das transições
# ----------------------------------------------------------

g_dry_frequency <- igraph::graph_from_adjacency_matrix(
  mat_dry_no_loop,
  mode = "directed",
  weighted = TRUE,
  diag = FALSE
)


# Frequências são forças de conexão.
# Para betweenness, precisamos convertê-las em distância (é uma distância ou resistência topológica derivada da intensidade de uso).
# maior frequência = menor distância.

if (igraph::ecount(g_dry_frequency) > 0) {
  
  igraph::E(g_dry_frequency)$distance <-
    1 / igraph::E(g_dry_frequency)$weight
}


# ----------------------------------------------------------
# Grafo não direcionado para clustering
# ----------------------------------------------------------

g_dry_undirected <- igraph::as_undirected(
  g_dry_binary,
  mode = "collapse"
)

g_dry_undirected <- igraph::simplify(
  g_dry_undirected,
  remove.multiple = TRUE,
  remove.loops = TRUE
)


# ----------------------------------------------------------
# Self-loop
# ----------------------------------------------------------

self_loop_count_dry <- diag(mat_dry)

self_loop_prop_total_dry <-
  self_loop_count_dry / total_transitions_dry

departures_dry <- rowSums(mat_dry)

self_loop_rate_dry <- ifelse(
  departures_dry > 0,
  self_loop_count_dry / departures_dry,
  NA_real_
)


# ----------------------------------------------------------
# Degree
# ----------------------------------------------------------

degree_in_dry <- igraph::degree(
  g_dry_binary,
  mode = "in",
  loops = FALSE
)

degree_out_dry <- igraph::degree(
  g_dry_binary,
  mode = "out",
  loops = FALSE
)

degree_total_dry <-
  degree_in_dry + degree_out_dry


# Em uma rede direcionada, o maior grau total possível sem loops é 2 × (n - 1).

degree_total_norm_dry <- if (n_nodes_dry > 1) {
  
  degree_total_dry /
    (2 * (n_nodes_dry - 1))
  
} else {
  
  rep(NA_real_, n_nodes_dry)
}


# ----------------------------------------------------------
# Betweenness não ponderado
# ----------------------------------------------------------

# Um nó possui alta betweenness quando funciona como uma ponte entre diferentes regiões da rede

betweenness_unweighted_dry <-
  igraph::betweenness(
    g_dry_binary,
    directed = TRUE,
    weights = NA,
    normalized = TRUE
  )


# ----------------------------------------------------------
# Betweenness ponderado pelo inverso da frequência
# ----------------------------------------------------------

if (igraph::ecount(g_dry_frequency) > 0) {
  
  betweenness_inverse_frequency_dry <-
    igraph::betweenness(
      g_dry_frequency,
      directed = TRUE,
      weights =
        igraph::E(g_dry_frequency)$distance,
      normalized = TRUE
    )
  
} else {
  
  betweenness_inverse_frequency_dry <-
    rep(0, n_nodes_dry)
}


# ----------------------------------------------------------
# Coeficiente de agrupamento local
# ----------------------------------------------------------

clustering_dry <- igraph::transitivity(
  g_dry_undirected,
  type = "local",
  isolates = "zero"
)


# ----------------------------------------------------------
# Tabela final da estação seca
# ----------------------------------------------------------

metrics_dry <- node_lookup_dry %>%
  transmute(
    node_id,
    cell_id,
    x,
    y,
    season = "Dry",
    
    transitions_from_node = departures_dry,
    
    self_loop_count =
      self_loop_count_dry,
    
    self_loop_prop_total =
      self_loop_prop_total_dry,
    
    self_loop_rate =
      self_loop_rate_dry,
    
    degree_in =
      degree_in_dry,
    
    degree_out =
      degree_out_dry,
    
    degree_total =
      degree_total_dry,
    
    degree_total_norm =
      degree_total_norm_dry,
    
    betweenness_unweighted_norm =
      betweenness_unweighted_dry,
    
    betweenness_inverse_frequency_norm =
      betweenness_inverse_frequency_dry,
    
    clustering =
      clustering_dry
  ) %>%
  
  left_join(
    metrics_moveNT_dry %>%
      select(
        cell_id,
        weight_moveNT,
        self_loop_moveNT,
        degree_moveNT,
        betweenness_moveNT,
        transitivity_moveNT
      ),
    by = "cell_id"
  )


write.csv2(metrics_dry,"path/Metrics_dry_turkey.csv", row.names = FALSE)


############################################################
# MÉTRICAS — RAINY SEASON
############################################################

mat_rainy <- adj_rainy[[1]]

n_nodes_rainy <- nrow(mat_rainy)
total_transitions_rainy <- sum(mat_rainy)


# ----------------------------------------------------------
# Matriz sem self-loops
# ----------------------------------------------------------

mat_rainy_no_loop <- mat_rainy

diag(mat_rainy_no_loop) <- 0


# ----------------------------------------------------------
# Matriz binária
# ----------------------------------------------------------

mat_rainy_binary <- 1L * (mat_rainy_no_loop > 0)


# ----------------------------------------------------------
# Grafo binário direcionado
# ----------------------------------------------------------

g_rainy_binary <- igraph::graph_from_adjacency_matrix(
  mat_rainy_binary,
  mode = "directed",
  weighted = NULL,
  diag = FALSE
)


# ----------------------------------------------------------
# Grafo ponderado pela frequência
# ----------------------------------------------------------

g_rainy_frequency <- igraph::graph_from_adjacency_matrix(
  mat_rainy_no_loop,
  mode = "directed",
  weighted = TRUE,
  diag = FALSE
)


if (igraph::ecount(g_rainy_frequency) > 0) {
  
  igraph::E(g_rainy_frequency)$distance <-
    1 / igraph::E(g_rainy_frequency)$weight
}


# ----------------------------------------------------------
# Grafo não direcionado para clustering
# ----------------------------------------------------------

g_rainy_undirected <- igraph::as_undirected(
  g_rainy_binary,
  mode = "collapse"
)

g_rainy_undirected <- igraph::simplify(
  g_rainy_undirected,
  remove.multiple = TRUE,
  remove.loops = TRUE
)


# ----------------------------------------------------------
# Self-loop
# ----------------------------------------------------------

self_loop_count_rainy <- diag(mat_rainy)

self_loop_prop_total_rainy <-
  self_loop_count_rainy / total_transitions_rainy

departures_rainy <- rowSums(mat_rainy)

self_loop_rate_rainy <- ifelse(
  departures_rainy > 0,
  self_loop_count_rainy / departures_rainy,
  NA_real_
)


# ----------------------------------------------------------
# Degree
# ----------------------------------------------------------

degree_in_rainy <- igraph::degree(
  g_rainy_binary,
  mode = "in",
  loops = FALSE
)

degree_out_rainy <- igraph::degree(
  g_rainy_binary,
  mode = "out",
  loops = FALSE
)

degree_total_rainy <-
  degree_in_rainy + degree_out_rainy


degree_total_norm_rainy <- if (n_nodes_rainy > 1) {
  
  degree_total_rainy /
    (2 * (n_nodes_rainy - 1))
  
} else {
  
  rep(NA_real_, n_nodes_rainy)
}


# ----------------------------------------------------------
# Betweenness não ponderado
# ----------------------------------------------------------

betweenness_unweighted_rainy <-
  igraph::betweenness(
    g_rainy_binary,
    directed = TRUE,
    weights = NA,
    normalized = TRUE
  )


# ----------------------------------------------------------
# Betweenness ponderado pelo inverso da frequência
# ----------------------------------------------------------

if (igraph::ecount(g_rainy_frequency) > 0) {
  
  betweenness_inverse_frequency_rainy <-
    igraph::betweenness(
      g_rainy_frequency,
      directed = TRUE,
      weights =
        igraph::E(g_rainy_frequency)$distance,
      normalized = TRUE
    )
  
} else {
  
  betweenness_inverse_frequency_rainy <-
    rep(0, n_nodes_rainy)
}


# ----------------------------------------------------------
# Coeficiente de agrupamento
# ----------------------------------------------------------

clustering_rainy <- igraph::transitivity(
  g_rainy_undirected,
  type = "local",
  isolates = "zero"
)


# ----------------------------------------------------------
# Tabela final da estação chuvosa
# ----------------------------------------------------------

metrics_rainy <- node_lookup_rainy %>%
  transmute(
    node_id,
    cell_id,
    x,
    y,
    season = "Rainy",
    
    transitions_from_node =
      departures_rainy,
    
    self_loop_count = # Número de permanências na mesma célula
      self_loop_count_rainy,
    
    self_loop_prop_total = # Contribuição da célula ao total de transições. Usada para identificar as principais manchas de permanência em toda a rede
      self_loop_prop_total_rainy,
    
    self_loop_rate = # Proporção das transições iniciadas na célula que permaneceram nela. Usada para identificar manchas com elevada tendência local de retenção
      self_loop_rate_rainy,
    
    degree_in =
      degree_in_rainy,
    
    degree_out =
      degree_out_rainy,
    
    degree_total =
      degree_total_rainy,
    
    degree_total_norm = # Número de manchas distintas conectadas
      degree_total_norm_rainy,
    
    betweenness_unweighted_norm = # Caminhos mínimos considerando apenas a existência das conexões. Identifica pontes estruturais da rede: manchas necessárias para conectar diferentes regiões usando o menor número de etapas.
      betweenness_unweighted_rainy,
    
    betweenness_inverse_frequency_norm = # Caminhos preferenciais considerando frequência das transições. identifica manchas localizadas em rotas de movimento intensamente utilizadas, e não apenas em rotas com poucas etapas.
      betweenness_inverse_frequency_rainy,
    
    clustering = # Formação de triângulos entre manchas vizinhas
      clustering_rainy
  ) %>%
  
  left_join(
    metrics_moveNT_rainy %>%
      select(
        cell_id,
        weight_moveNT,
        self_loop_moveNT,
        degree_moveNT,
        betweenness_moveNT,
        transitivity_moveNT
      ),
    by = "cell_id"
  )


write.csv2(metrics_rainy,"paths/Metrics_rainy_turkey.csv", row.names = FALSE)


############################################################
# TRANSIÇÕES COM INFORMAÇÕES TEMPORAIS
## Este bloco utiliza os burst_id já criados no código-base
############################################################

# ----------------------------------------------------------
# DRY SEASON
# ----------------------------------------------------------

transitions_dry <- dry_season %>%
  group_by(id, burst_id) %>%
  arrange(time_utc, .by_group = TRUE) %>%
  mutate(
    cell_to = lead(cell_id),
    time_to = lead(time_utc),
    time_local_to = lead(time_local)
  ) %>%
  ungroup() %>%
  filter(
    !is.na(cell_to),
    !is.na(time_to)
  ) %>%
  transmute(
    id,
    burst_id,
    season = "Dry",
    
    time_from = time_utc,
    time_to = time_to,
    
    time_local_from = time_local,
    time_local_to = time_local_to,
    
    cell_from = cell_id,
    cell_to = cell_to
  )


# ----------------------------------------------------------
# RAINY SEASON
# ----------------------------------------------------------

transitions_rainy <- rainy_season %>%
  group_by(id, burst_id) %>%
  arrange(time_utc, .by_group = TRUE) %>%
  mutate(
    cell_to = lead(cell_id),
    time_to = lead(time_utc),
    time_local_to = lead(time_local)
  ) %>%
  ungroup() %>%
  filter(
    !is.na(cell_to),
    !is.na(time_to)
  ) %>%
  transmute(
    id,
    burst_id,
    season = "Rainy",
    
    time_from = time_utc,
    time_to = time_to,
    
    time_local_from = time_local,
    time_local_to = time_local_to,
    
    cell_from = cell_id,
    cell_to = cell_to
  )


############################################################
# CLASSIFICAR OS HORÁRIOS DAS TRANSIÇÕES
############################################################

## Como agora usamos time_local, os horários devem ser os horários locais da Bahia:
# Roosting: 04:00–06:00 e 18:00–20:30;
# Feeding/Resting: 06:00–18:00;
# 1 h = 60 min

transitions_all <- bind_rows(
  transitions_dry,
  transitions_rainy
) %>%
  mutate(
# Converter o horário local em minutos após a meia-noite
    minute_from =
      hour(time_local_from) * 60 +
      minute(time_local_from),
    
    minute_to =
      hour(time_local_to) * 60 +
      minute(time_local_to),
    
    # Classificação do ponto inicial da transição    
    activity_from = case_when(
      # 04:00–06:00 ou 18:00–20:30  
      between(minute_from, 240, 360) |
        between(minute_from, 1080, 1230) ~
        "Roosting",
      # 06:00–18:00
      between(minute_from, 360, 1080) ~
        "Feeding",
      
      TRUE ~ "Other"
    ),
    
    # Classificação do ponto final da transição  
    activity_to = case_when(
      # 04:00–06:00 ou 18:00–20:30
      between(minute_to, 240, 360) |
        between(minute_to, 1080, 1230) ~
        "Roosting",
      # 06:00–18:00
      between(minute_to, 360, 1080) ~
        "Feeding",
      
      TRUE ~ "Other"
    ),
    
    
    # A transição só recebe uma classe quando os dois
    # registros pertencem ao mesmo período de atividade.
    activity = dplyr::if_else(
      activity_from == activity_to,
      activity_from,
      NA_character_
    )
  )


## Uma transição só recebe uma atividade quando os dois pontos pertencem ao mesmo período, ou seja, Uma transição iniciada às 17:50 e terminada às 18:20, por exemplo, não será classificada nem como feeding nem como roosting.


############################################################
# SELF-LOOPS DE ROOSTING E FEEDING
############################################################

activity_transitions <- transitions_all %>%
  filter(
    activity %in% c("Roosting", "Feeding")
  )


# ----------------------------------------------------------
# Células envolvidas em cada atividade
# ----------------------------------------------------------

activity_cells <- bind_rows(
  
  activity_transitions %>%
    transmute(
      season,
      activity,
      cell_id = cell_from
    ),
  
  activity_transitions %>%
    transmute(
      season,
      activity,
      cell_id = cell_to
    )
  
) %>%
  distinct()


# ----------------------------------------------------------
# Número de transições iniciadas em cada célula
# ----------------------------------------------------------

activity_departures <- activity_transitions %>%
  count(
    season,
    activity,
    cell_id = cell_from,
    name = "n_departures"
  )


# ----------------------------------------------------------
# Número de self-loops em cada célula
# ----------------------------------------------------------

activity_self_loops <- activity_transitions %>%
  filter(cell_from == cell_to) %>%
  count(
    season,
    activity,
    cell_id = cell_from,
    name = "self_loop_count"
  )


# ----------------------------------------------------------
# Total de transições por estação e atividade
# ----------------------------------------------------------

activity_totals <- activity_transitions %>%
  count(
    season,
    activity,
    name = "n_transitions_activity"
  )


# ----------------------------------------------------------
# Coordenadas permanentes das células
# ----------------------------------------------------------

all_activity_cells <- sort(
  unique(activity_cells$cell_id)
)

activity_xy <- raster::xyFromCell(
  grid_raster,
  all_activity_cells
)

activity_coordinates <- data.frame(
  cell_id = all_activity_cells,
  x = activity_xy[, 1],
  y = activity_xy[, 2]
)


# ----------------------------------------------------------
# Tabela final
# ----------------------------------------------------------

self_loop_activity <- activity_cells %>%
  
  left_join(
    activity_departures,
    by = c(
      "season",
      "activity",
      "cell_id"
    )
  ) %>%
  
  left_join(
    activity_self_loops,
    by = c(
      "season",
      "activity",
      "cell_id"
    )
  ) %>%
  
  left_join(
    activity_totals,
    by = c(
      "season",
      "activity"
    )
  ) %>%
  
  left_join(
    activity_coordinates,
    by = "cell_id"
  ) %>%
  
  dplyr::mutate(
    
    n_departures =
      coalesce(n_departures, 0L),
    
    self_loop_count =
      coalesce(self_loop_count, 0L),
    
    
    # Equivalente conceitual ao self-loop do moveNT,
    # mas considerando somente aquela atividade.
    
    self_loop_prop_activity = # Representa a contribuição de cada célula ao total de self-loops daquele período de atividade.
      self_loop_count /
      n_transitions_activity,
    
    
    # Probabilidade de permanecer na mesma célula,
    # condicionada ao uso daquela célula.
    
    self_loop_rate = dplyr::if_else(
      n_departures > 0,
      self_loop_count / n_departures,
      NA_real_
    )
  )


## Separar os quatro data frames

df_self_dry_roost <- self_loop_activity %>%
  filter(
    season == "Dry",
    activity == "Roosting"
  )

write.csv2(df_self_dry_roost,"C:/Users/danip/Desktop/Pos doc/Projeto/Paper_3/Results/Self_dry_roost_turkey.csv", row.names = FALSE)


df_self_rainy_roost <- self_loop_activity %>%
  filter(
    season == "Rainy",
    activity == "Roosting"
  )

write.csv2(df_self_rainy_roost,"path/4.Self_rainy_roost_turkey.csv", row.names = FALSE)


df_self_dry_feeding <- self_loop_activity %>%
  filter(
    season == "Dry",
    activity == "Feeding"
  )

write.csv2(df_self_dry_feeding,"path/4.Self_dry_feeding_turkey.csv", row.names = FALSE)


df_self_rainy_feeding <- self_loop_activity %>%
  filter(
    season == "Rainy",
    activity == "Feeding"
  )


write.csv2(df_self_rainy_feeding,"path/4.Self_rainy_feeding_turkey.csv", row.names = FALSE)


############################################################
# COMPARAÇÃO ESPACIAL ENTRE ESTAÇÕES
############################################################

## Criar uma tabela para comparar as mesmas manchas

metrics_comparison <- metrics_dry %>%
  transmute(
    cell_id,
    x,
    y,
    
    used_dry = TRUE,
    
    self_loop_dry =
      self_loop_prop_total,
    
    self_loop_rate_dry =
      self_loop_rate,
    
    degree_dry =
      degree_total_norm,
    
    betweenness_dry =
      betweenness_unweighted_norm,
    
    clustering_dry =
      clustering
  ) %>%
  
  full_join(
    
    metrics_rainy %>%
      transmute(
        cell_id,
        x,
        y,
        
        used_rainy = TRUE,
        
        self_loop_rainy =
          self_loop_prop_total,
        
        self_loop_rate_rainy =
          self_loop_rate,
        
        degree_rainy =
          degree_total_norm,
        
        betweenness_rainy =
          betweenness_unweighted_norm,
        
        clustering_rainy =
          clustering
      ),
    
    by = c(
      "cell_id",
      "x",
      "y"
    )
  ) %>%
  
  mutate(
    
    used_dry =
      coalesce(used_dry, FALSE),
    
    used_rainy =
      coalesce(used_rainy, FALSE),
    
    
    use_category = case_when(
      
      used_dry & used_rainy ~
        "Both seasons",
      
      used_dry & !used_rainy ~
        "Dry only",
      
      !used_dry & used_rainy ~
        "Rainy only"
    ),
    
    
    self_loop_change =
      self_loop_rainy -
      self_loop_dry,
    
    degree_change =
      degree_rainy -
      degree_dry,
    
    betweenness_change =
      betweenness_rainy -
      betweenness_dry,
    
    clustering_change =
      clustering_rainy -
      clustering_dry
  )


write.csv2(metrics_comparison,"path/4.Comparacao_entre_as_metrics_por_no_turkey.csv", row.names = FALSE)


# Os valores não devem ser transformados em zero quando uma célula não foi usada em determinada estação. Nesse caso, a centralidade é NA, porque o nó não pertencia à rede daquela estação. A ausência de uso fica representada separadamente por:
  
# used_dry
# used_rainy
# use_category

# Agora já temos métricas comparáveis:
# self_loop_prop_total: proporção;
# self_loop_rate: entre 0 e 1;
# degree_total_norm: entre 0 e 1;
# betweenness_unweighted_norm: normalizado;
# clustering: entre 0 e 1.


## Mapas de grau, betweenness e clustering

#  grid_bbox
#  xmin    ymin     xmax    ymax 
# 432000  8752000  664000  9046000 

# PALETA

pal_blue <- colorRampPalette(
  c(
    "#ece7f2",
    "#7fcdbb",
    "#41b6c4",
    "#253494",
    "#081d58"
  )
)


# DEGREE

p_degree_dry <- ggplot(
  metrics_dry,
  aes(
    x = x,
    y = y,
    fill = degree_total_norm
  )
) +
  geom_tile(
    width = grid_spacing,
    height = grid_spacing
  ) +
  scale_fill_gradientn(
    colors = pal_blue(5),
    limits = c(0, 0.12)
  ) +
  coord_equal(
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    expand = FALSE
  ) +
  labs(
    x = "Longitude (m)",
    y = "",
    fill = "Normalized\ndegree"
  ) +
  theme_bw()


p_degree_rainy <- ggplot(
  metrics_rainy,
  aes(
    x = x,
    y = y,
    fill = degree_total_norm
  )
) +
  geom_tile(
    width = grid_spacing,
    height = grid_spacing
  ) +
  scale_fill_gradientn(
    colors = pal_blue(5),
    limits = c(0, 0.12)
  ) +
  coord_equal(
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    expand = FALSE
  ) +
  labs(
    x = "Longitude (m)",
    y = "Latitude (m)",
    fill = "Normalized\ndegree"
  ) +
  theme_bw()


ggarrange(
  p_degree_rainy,
  p_degree_dry,
  nrow = 1,
  ncol = 2,
  labels = c("a", "b"),
  legend = "right",
  common.legend = TRUE
)




# BETWEENNESS

p_bet_dry <- ggplot(
  metrics_dry,
  aes(
    x = x,
    y = y,
    fill = betweenness_unweighted_norm
  )
) +
  geom_tile(
    width = grid_spacing,
    height = grid_spacing
  ) +
  scale_fill_gradientn(
    colors = pal_blue(5),
    limits = c(0, 0.32)
  ) +
  coord_equal(
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    expand = FALSE
  ) +
  labs(
    x = "Longitude (m)",
    y = "",
    fill = "Normalized\nbetweenness"
  ) +
  theme_bw()


p_bet_rainy <- ggplot(
  metrics_rainy,
  aes(
    x = x,
    y = y,
    fill = betweenness_unweighted_norm
  )
) +
  geom_tile(
    width = grid_spacing,
    height = grid_spacing
  ) +
  scale_fill_gradientn(
    colors = pal_blue(5),
    limits = c(0, 0.32)
  ) +
  coord_equal(
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    expand = FALSE
  ) +
  labs(
    x = "Longitude (m)",
    y = "Latitude (m)",
    fill = "Normalized\nbetweenness"
  ) +
  theme_bw()


ggarrange(
  p_bet_rainy,
  p_bet_dry,
  nrow = 1,
  ncol = 2,
  labels = c("a", "b"),
  legend = "right",
  common.legend = TRUE
)



# CLUSTERING COEFFICIENT

pal_blue_c <- colorRampPalette(c( "#f7fbff", "#deebf7","#deebf7",
                                  "#081d58"))

p_clustering_dry <- ggplot(
  metrics_dry,
  aes(
    x = x,
    y = y,
    fill = clustering
  )
) +
  geom_tile(
    width = grid_spacing,
    height = grid_spacing
  ) +
  scale_fill_gradientn(
    colors = pal_blue_c(4),
    limits = c(0, 1)
  ) +
  coord_equal(
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    expand = FALSE
  ) +
  labs(
    x = "Longitude (m)",
    y = "",
    fill = "Clustering"
  ) +
  theme_bw()


p_clustering_rainy <- ggplot(
  metrics_rainy,
  aes(
    x = x,
    y = y,
    fill = clustering
  )
) +
  geom_tile(
    width = grid_spacing,
    height = grid_spacing
  ) +
  scale_fill_gradientn(
    colors = pal_blue_c(4),
    limits = c(0, 1)
  ) +
  coord_equal(
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    expand = FALSE
  ) +
  labs(
    x = "Longitude (m)",
    y = "Latitude (m)",
    fill = "Clustering"
  ) +
  theme_bw()


ggarrange(
  p_clustering_rainy,
  p_clustering_dry,
  nrow = 1,
  ncol = 2,
  labels = c("a", "b"),
  legend = "right",
  common.legend = TRUE
)



### Mapas de self-loop

#######################
# SELF-LOOP — ROOSTING
#######################

roosting_max <- max(
  self_loop_activity$self_loop_prop_activity[
    self_loop_activity$activity == "Roosting"
  ],
  na.rm = TRUE
)


p_roosting_dry <- ggplot(
  df_self_dry_roost,
  aes(
    x = x,
    y = y,
    fill = self_loop_prop_activity
  )
) +
  geom_tile(
    width = grid_spacing,
    height = grid_spacing
  ) +
  scale_fill_gradientn(
    colors = pal_blue(5),
    limits = c(0, roosting_max)
  ) +
  coord_equal(
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    expand = FALSE
  ) +
  labs(
    x = "",
    y = "",
    fill = "Roosting\nself-loop"
  ) +
  theme_bw()


p_roosting_rainy <- ggplot(
  df_self_rainy_roost,
  aes(
    x = x,
    y = y,
    fill = self_loop_prop_activity
  )
) +
  geom_tile(
    width = grid_spacing,
    height = grid_spacing
  ) +
  scale_fill_gradientn(
    colors = pal_blue(5),
    limits = c(0, roosting_max)
  ) +
  coord_equal(
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    expand = FALSE
  ) +
  labs(
    x = "",
    y = "Latitude (m)",
    fill = "Roosting\nself-loop"
  ) +
  theme_bw()


ggarrange(
  p_roosting_rainy,
  p_roosting_dry,
  nrow = 1,
  ncol = 2,
  labels = c("a", "b"),
  legend = "right",
  common.legend = TRUE
)



###############################
# SELF-LOOP — FEEDING/RESTING
###############################

feeding_max <- max(
  self_loop_activity$self_loop_prop_activity[
    self_loop_activity$activity == "Feeding"
  ],
  na.rm = TRUE
)


p_feeding_dry <- ggplot(
  df_self_dry_feeding,
  aes(
    x = x,
    y = y,
    fill = self_loop_prop_activity
  )
) +
  geom_tile(
    width = grid_spacing,
    height = grid_spacing
  ) +
  scale_fill_gradientn(
    colors = pal_blue(5),
    limits = c(0, feeding_max)
  ) +
  coord_equal(
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    expand = FALSE
  ) +
  labs(
    x = "Longitude (m)",
    y = "",
    fill = "Feeding\nself-loop"
  ) +
  theme_bw()


p_feeding_rainy <- ggplot(
  df_self_rainy_feeding,
  aes(
    x = x,
    y = y,
    fill = self_loop_prop_activity
  )
) +
  geom_tile(
    width = grid_spacing,
    height = grid_spacing
  ) +
  scale_fill_gradientn(
    colors = pal_blue(5),
    limits = c(0, feeding_max)
  ) +
  coord_equal(
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    expand = FALSE
  ) +
  labs(
    x = "Longitude (m)",
    y = "Latitude (m)",
    fill = "Feeding\nself-loop"
  ) +
  theme_bw()


ggarrange(p_feeding_rainy, p_feeding_dry,
  nrow = 1, ncol = 2, labels = c("c", "d"),
  legend = "right", common.legend = T)




### Identificar as manchas centrais
## Como primeiro critério descritivo, eu primeiramente defini como áreas chave aqueles nós situados nos 10% superiores de cada métrica dentro de cada estação. Esse critério permite verificar diretamente quais manchas são centrais apenas na seca, apenas na chuva e simultaneamente nas duas estações. 

############################################################
# IDENTIFICANDO NÓS CENTRAIS
############################################################

metrics_all <- bind_rows(
  metrics_dry,
  metrics_rainy
) %>%
  group_by(season) %>%
  mutate(
    
    self_loop_percentile =
      percent_rank(self_loop_prop_total),
    
    degree_percentile =
      percent_rank(degree_total_norm),
    
    betweenness_percentile =
      percent_rank(
        betweenness_unweighted_norm
      ),
    
    clustering_percentile =
      percent_rank(clustering),
    
    
    key_self_loop =
      self_loop_percentile >= 0.9, # 10%
    
    key_degree =
      degree_percentile >= 0.9,
    
    key_betweenness =
      betweenness_percentile >= 0.9,
    
    key_clustering =
      clustering_percentile >= 0.9 
  ) %>%
  ungroup()


write.csv2(metrics_all,"path/4.Central_nodes_by_metric_10%_turkey.csv", row.names = FALSE)



######################################################
##### Plotando os histogramas das métricas por estação
######################################################

cor_seca  <- "orange"
cor_chuva <- "lightblue"


## Preparando os vetores das métricas

### ESTAÇÃO SECA ####

self_dry <- metrics_dry %>%
  dplyr::select(self_loop_prop_total) %>%
  filter(!is.na(self_loop_prop_total))

degree_dry_df <- metrics_dry %>%
  dplyr::select(degree_total_norm) %>%
  filter(!is.na(degree_total_norm))

bet_dry_df <- metrics_dry %>%
  dplyr::select(betweenness_unweighted_norm) %>%
  filter(!is.na(betweenness_unweighted_norm))

clust_dry_df <- metrics_dry %>%
  dplyr::select(clustering) %>%
  filter(!is.na(clustering))

# Self-loop centrality
(p_self_dry <- ggplot(self_dry, aes(x = self_loop_prop_total)) +
  geom_histogram(fill = cor_seca, color = "black") +
  labs(
    x = "Self-loop",
    y = "") +
  #scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  #scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 10)) +
  theme_classic(base_size = 14))


# Degree centrality
(p_degree_dry <- ggplot(degree_dry_df, aes(x = degree_total_norm)) +
    geom_histogram(fill = cor_seca, color = "black") +
    labs(
      x = "Degree",
      y = "Habitat patches"
    ) +
    theme_classic(base_size = 14))

# Betweenness centrality
(p_bet_dry <- ggplot(bet_dry_df, aes(x = betweenness_unweighted_norm)) +
  geom_histogram(fill = cor_seca, color = "black") +
  labs(
    x = "Betweenness",
    y = ""
  ) +
  theme_classic(base_size = 14))


# Clustering coefficient 
(p_clust_dry <- ggplot(clust_dry_df, aes(x = clustering)) +
  geom_histogram(fill = cor_seca, color = "black") +
  labs(
    x = "Clustering coefficient",
    y = ""
  ) +
  #scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  #scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_classic(base_size = 14))



### Estação chuvosa ####

self_rainy <- metrics_rainy %>%
  dplyr::select(self_loop_prop_total) %>%
  filter(!is.na(self_loop_prop_total))

degree_rainy_df <- metrics_rainy %>%
  dplyr::select(degree_total_norm) %>%
  filter(!is.na(degree_total_norm))

bet_rainy_df <- metrics_rainy %>%
  dplyr::select(betweenness_unweighted_norm) %>%
  filter(!is.na(betweenness_unweighted_norm))

clust_rainy_df <- metrics_rainy %>%
  dplyr::select(clustering) %>%
  filter(!is.na(clustering))


# Self-loop centrality
(p_self_rainy <- ggplot(self_rainy, aes(x = self_loop_prop_total)) +
  geom_histogram(fill = cor_chuva, color = "black") +
  labs(
    x = "Self-loop",
    y = "Habitat patches"
  ) +
  #scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  #scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 10)) +
  theme_classic(base_size = 14))

# Degree centrality
(p_degree_rainy <- ggplot(degree_rainy_df, aes(x = degree_total_norm)) +
  geom_histogram(fill = cor_chuva, color = "black") +
  labs(
    x = "Degree",
    y = "Habitat patches"
  ) +
  theme_classic(base_size = 14))


# Betweenness centrality
(p_bet_rainy <- ggplot(bet_rainy_df, aes(x = betweenness_unweighted_norm)) +
  geom_histogram(fill = cor_chuva, color = "black") +
  labs(
    x = "Betweenness",
    y = "Habitat patches"
  ) +
  theme_classic(base_size = 14))


# Clustering
(p_clust_rainy <- ggplot(clust_rainy_df, aes(x = clustering)) +
  geom_histogram(fill = cor_chuva, color = "black") +
  labs(
    x = "Clustering coefficient",
    y = "Habitat patches"
  ) +
  #scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  #scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_classic(base_size = 14))


## Figures 
(fig1 <- ggarrange(
  p_degree_rainy, p_degree_dry, 
  ncol = 2, nrow = 1,
  labels = c("a", "b")
))


(fig2 <- ggarrange(
  p_self_rainy, p_self_dry, 
  p_bet_rainy, p_bet_dry, 
  p_clust_rainy, p_clust_dry, 
  ncol = 2, nrow = 3,
  labels = c("a", "b", "c", "d", "e", "f")))




### Para calcular a sobreposição das manchas centrais eu usei o índice de Jaccard, onde: J = número de manchas centrais nas duas estações / n° total de manchas centrais em pelo menos uma estação

# Sobreposição das manchas centrais entre estações

############################################################
# ÍNDICE DE JACCARD DOS NÓS CENTRAIS
############################################################

dry_self_loop_key <- metrics_all %>%
  filter(
    season == "Dry",
    key_self_loop
  ) %>%
  pull(cell_id)

rainy_self_loop_key <- metrics_all %>%
  filter(
    season == "Rainy",
    key_self_loop
  ) %>%
  pull(cell_id)


dry_degree_key <- metrics_all %>%
  filter(
    season == "Dry",
    key_degree
  ) %>%
  pull(cell_id)

rainy_degree_key <- metrics_all %>%
  filter(
    season == "Rainy",
    key_degree
  ) %>%
  pull(cell_id)


dry_betweenness_key <- metrics_all %>%
  filter(
    season == "Dry",
    key_betweenness
  ) %>%
  pull(cell_id)

rainy_betweenness_key <- metrics_all %>%
  filter(
    season == "Rainy",
    key_betweenness
  ) %>%
  pull(cell_id)


dry_clustering_key <- metrics_all %>%
  filter(
    season == "Dry",
    key_clustering
  ) %>%
  pull(cell_id)

rainy_clustering_key <- metrics_all %>%
  filter(
    season == "Rainy",
    key_clustering
  ) %>%
  pull(cell_id)


key_node_overlap <- data.frame(
  
  metric = c(
    "Self-loop",
    "Degree",
    "Betweenness",
    "Clustering"
  ),
  
  intersection = c(
    
    length(intersect(
      dry_self_loop_key,
      rainy_self_loop_key
    )),
    
    length(intersect(
      dry_degree_key,
      rainy_degree_key
    )),
    
    length(intersect(
      dry_betweenness_key,
      rainy_betweenness_key
    )),
    
    length(intersect(
      dry_clustering_key,
      rainy_clustering_key
    ))
  ),
  
  union = c(
    
    length(union(
      dry_self_loop_key,
      rainy_self_loop_key
    )),
    
    length(union(
      dry_degree_key,
      rainy_degree_key
    )),
    
    length(union(
      dry_betweenness_key,
      rainy_betweenness_key
    )),
    
    length(union(
      dry_clustering_key,
      rainy_clustering_key
    ))
  )
)


(key_node_overlap$jaccard <-
  key_node_overlap$intersection /
  key_node_overlap$union)


write.csv2(key_node_overlap,"path/4.Jaccard_key_areas_10%_turkey.csv", row.names = FALSE)

#    metric     intersection  union   jaccard (10% das métricas com maiores valores)
# 1 Self-loop       33         182   0.18131868
# 2 Degree          39         180   0.21666667
# 3 Betweenness     7          216   0.03240741
# 4 Clustering      29         193   0.15025907

# Interpretação:
# Valores de Jaccard próximos de:
# 1: as mesmas manchas são centrais nas duas estações;
# 0: as manchas centrais são completamente diferentes;
# valores intermediários: sobreposição parcial.

# Esse índice é um resumo descritivo inicial. A comparação estatística ainda precisará considerar a diferença de esforço amostral entre seca e chuva, especialmente para self-loop e medidas dependentes do número de transições.



##### JACCARD PRA SELF LOOP ROOSTING AND FEEDING SITES

############################################################
# IDENTIFICAR MANCHAS-CHAVE DE ROOSTING E FEEDING
############################################################

# As manchas-chave são definidas separadamente dentro de cada combinação de estação e atividade: season × activity
# Critérios: 10% dos valores superiores da distribuição de self-loop percentil >= 0,9

self_loop_activity_key <- self_loop_activity %>%
  
  dplyr::filter(
    activity %in%
      c(
        "Roosting",
        "Feeding"
      )
  ) %>%
  
  dplyr::group_by(
    season,
    activity
  ) %>%
  
  dplyr::mutate(
    
    ########################################################
    # Posição percentílica dentro de season × activity
    ########################################################
    
    self_loop_percentile_activity =
      dplyr::percent_rank(
        self_loop_prop_activity
      ),
    
    
    ########################################################
    # Percentil de corte específico de cada atividade
    ########################################################
    
    activity_percentile_cutoff =
      dplyr::case_when(
        
        activity == "Roosting" ~
          0.9,
        
        activity == "Feeding" ~
          0.9,
        
        TRUE ~
          NA_real_
      ),
    
    
    ########################################################
    # Classificação das manchas-chave
    ########################################################
    
    key_self_loop_activity =
      
      self_loop_percentile_activity >=
      activity_percentile_cutoff &
      
      self_loop_count > 0
  ) %>%
  
  dplyr::ungroup()


# Conferência dos cortes utilizados

self_loop_activity_key %>%
  
  dplyr::distinct(
    season,
    activity,
    activity_percentile_cutoff
  ) %>%
  
  dplyr::arrange(
    season,
    activity
  )

## O esperado:
# season  activity   activity_percentile_cutoff
# Dry     Feeding            0.9
# Dry     Roosting           0.9
# Rainy   Feeding            0.9
# Rainy   Roosting           0.9


############################################################
# CONTRIBUIÇÃO AO TOTAL DE SELF-LOOPS DA ATIVIDADE
############################################################

self_loop_activity_key <- self_loop_activity_key %>%
  
  dplyr::group_by(
    season,
    activity
  ) %>%
  
  dplyr::mutate(
    
    total_self_loops_activity =
      sum(
        self_loop_count,
        na.rm = TRUE
      ),
    
    
    self_loop_share_activity =
      dplyr::if_else(
        
        total_self_loops_activity > 0,
        
        self_loop_count /
          total_self_loops_activity,
        
        NA_real_
      )
  ) %>%
  
  dplyr::ungroup()


# self_loop_share_activity representa diretamente a proporção dos self-loops daquela atividade que ocorreu em cada célula. O conjunto dos 10% superiores será igual ao obtido com self_loop_prop_activity, porque ambas as variáveis diferem apenas por um denominador constante dentro de cada grupo.


############################################################
# RESUMO DAS MANCHAS-CHAVE
############################################################

# Conferir o número de manchas selecionadas
# Por causa de empates nos valores de self-loop, a quantidade selecionada pode não ser exatamente 10%. Por exemplo, se várias células tiverem o mesmo valor no limite do percentil, todas poderão ser classificadas juntas. Isso preserva o critério baseado no valor da métrica.

key_activity_summary <- self_loop_activity_key %>%
  
  dplyr::group_by(
    season,
    activity
  ) %>%
  
  dplyr::summarise(
    
    n_activity_cells =
      dplyr::n(),
    
    n_cells_with_self_loop =
      sum(
        self_loop_count > 0,
        na.rm = TRUE
      ),
    
    percentile_cutoff =
      dplyr::first(
        activity_percentile_cutoff
      ),
    
    n_key_cells =
      sum(
        key_self_loop_activity,
        na.rm = TRUE
      ),
    
    proportion_key_cells =
      n_key_cells /
      n_activity_cells,
    
    minimum_selected_percentile =
      if (
        any(
          key_self_loop_activity,
          na.rm = TRUE
        )
      ) {
        
        min(
          self_loop_percentile_activity[
            key_self_loop_activity
          ],
          na.rm = TRUE
        )
        
      } else {
        
        NA_real_
      },
    
    minimum_selected_self_loop =
      if (
        any(
          key_self_loop_activity,
          na.rm = TRUE
        )
      ) {
        
        min(
          self_loop_prop_activity[
            key_self_loop_activity
          ],
          na.rm = TRUE
        )
        
      } else {
        
        NA_real_
      },
    
    .groups = "drop"
  )


key_activity_summary

write.csv2(key_activity_summary,"path/4.key_activity_summary_self_loop_10%_turkey.csv",
  
  row.names = FALSE
)

############################################################
# TODAS AS MANCHAS OCUPADAS — CAMADA DE FUNDO
############################################################

# Criar uma tabela única que possui todas as células ocupadas em cada estação, independentemente de terem sido ou não classificadas como áreas-chave.

all_used_habitat_patches <- dplyr::bind_rows(
  
  metrics_dry %>%
    
    dplyr::transmute(
      season = "Dry",
      cell_id,
      x,
      y
    ),
  
  
  metrics_rainy %>%
    
    dplyr::transmute(
      season = "Rainy",
      cell_id,
      x,
      y
    )
) %>%
  
  dplyr::distinct(
    season,
    cell_id,
    .keep_all = TRUE
  )  


############################################################
# EXTRAIR OS CELL_ID DAS MANCHAS-CHAVE
############################################################

# ----------------------------------------------------------
# DRY — ROOSTING
# ----------------------------------------------------------

dry_roosting_key <- self_loop_activity_key %>%
  
  dplyr::filter(
    season == "Dry",
    activity == "Roosting",
    key_self_loop_activity
  ) %>%
  
  dplyr::pull(
    cell_id
  ) %>%
  
  unique()


# ----------------------------------------------------------
# DRY — FEEDING
# ----------------------------------------------------------

dry_feeding_key <- self_loop_activity_key %>%
  
  dplyr::filter(
    season == "Dry",
    activity == "Feeding",
    key_self_loop_activity
  ) %>%
  
  dplyr::pull(
    cell_id
  ) %>%
  
  unique()


# ----------------------------------------------------------
# RAINY — ROOSTING
# ----------------------------------------------------------

rainy_roosting_key <- self_loop_activity_key %>%
  
  dplyr::filter(
    season == "Rainy",
    activity == "Roosting",
    key_self_loop_activity
  ) %>%
  
  dplyr::pull(
    cell_id
  ) %>%
  
  unique()


# ----------------------------------------------------------
# RAINY — FEEDING
# ----------------------------------------------------------

rainy_feeding_key <- self_loop_activity_key %>%
  
  dplyr::filter(
    season == "Rainy",
    activity == "Feeding",
    key_self_loop_activity
  ) %>%
  
  dplyr::pull(
    cell_id
  ) %>%
  
  unique()

length(dry_roosting_key) # 9
length(dry_feeding_key) # 50
length(rainy_roosting_key) # 38
length(rainy_feeding_key) # 164

anyDuplicated(dry_roosting_key)
anyDuplicated(dry_feeding_key)
anyDuplicated(rainy_roosting_key)
anyDuplicated(rainy_feeding_key)



############################################################
# ÍNDICE DE JACCARD — ROOSTING VERSUS FEEDING
############################################################

jaccard_roosting_feeding <- data.frame(
  
  season = c(
    "Dry",
    "Rainy"
  ),
  
  
  n_roosting_key = c(
    length(
      dry_roosting_key
    ),
    
    length(
      rainy_roosting_key
    )
  ),
  
  
  n_feeding_key = c(
    length(
      dry_feeding_key
    ),
    
    length(
      rainy_feeding_key
    )
  ),
  
  
  intersection = c(
    
    length(
      intersect(
        dry_roosting_key,
        dry_feeding_key
      )
    ),
    
    length(
      intersect(
        rainy_roosting_key,
        rainy_feeding_key
      )
    )
  ),
  
  
  union = c(
    
    length(
      union(
        dry_roosting_key,
        dry_feeding_key
      )
    ),
    
    length(
      union(
        rainy_roosting_key,
        rainy_feeding_key
      )
    )
  )
)


## calculando o índice e algumas medidas complementares:

jaccard_roosting_feeding <- jaccard_roosting_feeding %>%
  
  dplyr::mutate(
    
    ########################################################
    # Jaccard
    ########################################################
    
    jaccard =
      dplyr::if_else(
        
        union > 0,
        
        intersection /
          union,
        
        NA_real_
      ),
    
    
    ########################################################
    # Não sobreposição ou turnover
    ########################################################
    
    turnover = 1 - jaccard,
    
    
    ########################################################
    # Proporção das manchas de roosting que também foram
    # classificadas como manchas-chave de feeding
    ########################################################
    
    proportion_roosting_also_feeding =
      dplyr::if_else(
        
        n_roosting_key > 0,
        
        intersection /
          n_roosting_key,
        
        NA_real_),
    
    
    ########################################################
    # Proporção das manchas de feeding que também foram
    # classificadas como manchas-chave de roosting
    ########################################################
    
    proportion_feeding_also_roosting =
      dplyr::if_else(
        
        n_feeding_key > 0,
        
        intersection /
          n_feeding_key,
        
        NA_real_
      ))


jaccard_roosting_feeding

write.csv2(jaccard_roosting_feeding, "path/4.Jaccard_roosting_feeding_by_season_self_loop_10%_turkey.csv",
           row.names = FALSE)


############################################################
# MANCHAS-CHAVE COMPARTILHADAS
############################################################

shared_dry_roosting_feeding <- intersect(
  dry_roosting_key,
  dry_feeding_key)


shared_rainy_roosting_feeding <- intersect(
  rainy_roosting_key,
  rainy_feeding_key)


shared_activity_key_cells <- dplyr::bind_rows(
  
  data.frame(
    season = "Dry",
    cell_id =
      shared_dry_roosting_feeding),
  
  data.frame(
    season = "Rainy",
    cell_id =
      shared_rainy_roosting_feeding))

shared_activity_key_cells <- shared_activity_key_cells %>%
  
  dplyr::left_join(
    activity_coordinates,
    by = "cell_id")


shared_activity_key_cells

write.csv2(shared_activity_key_cells, "path/4.shared_activity_key_cells_self_loop_10%_turkey.csv",row.names = FALSE)


## Classificar todas as manchas-chave por categoria. Este bloco cria uma tabela espacial mostrando quais manchas foram: 
# importantes apenas para roosting;
# importantes apenas para feeding;
# importantes para ambas as atividades.


############################################################
# CATEGORIA DE SOBREPOSIÇÃO — DRY
############################################################

dry_activity_overlap <- dplyr::full_join(
  
  data.frame(cell_id = dry_roosting_key,
    
    key_roosting = TRUE),
  
  data.frame(cell_id = dry_feeding_key, key_feeding = TRUE),by = "cell_id") %>%
  
  dplyr::mutate(season = "Dry", key_roosting =
      dplyr::coalesce(key_roosting, FALSE),
    key_feeding =
      dplyr::coalesce(
        key_feeding,
        FALSE),
    
    overlap_category =
      dplyr::case_when(
        
        key_roosting &
          key_feeding ~
          "Roosting and feeding",
        
        key_roosting &
          !key_feeding ~
          "Roosting only",
        
        !key_roosting &
          key_feeding ~
          "Feeding only"
      )
  ) %>%
  
  dplyr::left_join(
    activity_coordinates,
    by = "cell_id"
  )


############################################################
# CATEGORIA DE SOBREPOSIÇÃO — RAINY
############################################################

rainy_activity_overlap <- dplyr::full_join(
  
  data.frame(
    cell_id =
      rainy_roosting_key,
    
    key_roosting =
      TRUE
  ),
  
  data.frame(
    cell_id =
      rainy_feeding_key,
    
    key_feeding =
      TRUE
  ),
  
  by = "cell_id"
) %>%
  
  dplyr::mutate(
    
    season =
      "Rainy",
    
    key_roosting =
      dplyr::coalesce(
        key_roosting,
        FALSE
      ),
    
    key_feeding =
      dplyr::coalesce(
        key_feeding,
        FALSE
      ),
    
    overlap_category =
      dplyr::case_when(
        
        key_roosting &
          key_feeding ~
          "Roosting and feeding",
        
        key_roosting &
          !key_feeding ~
          "Roosting only",
        
        !key_roosting &
          key_feeding ~
          "Feeding only"
      )
  ) %>%
  
  dplyr::left_join(
    activity_coordinates,
    by = "cell_id"
  )


### Juntando
activity_overlap_cells <- dplyr::bind_rows(
  dry_activity_overlap,
  rainy_activity_overlap
)


activity_overlap_cells %>%
  
  dplyr::count(
    season,
    overlap_category
  )

write.csv2(activity_overlap_cells, "path/4.Jaccard_activity_overlap_cells_self_loop_10%_turkey.csv", row.names = FALSE)


############################################################
# MAPA — SOBREPOSIÇÃO ENTRE ROOSTING E FEEDING
############################################################

activity_overlap_cells$overlap_category <- factor(
  activity_overlap_cells$overlap_category,
  
  levels = c(
    "Roosting only",
    "Feeding only",
    "Roosting and feeding"
  )
)



############################################################
# MAPA — DRY SEASON
############################################################

p_overlap_activity_dry <- ggplot2::ggplot() +
  
##########################################################
# Camada de fundo: todas as manchas usadas na seca
##########################################################

ggplot2::geom_tile(
  
  data = all_used_habitat_patches %>%
    dplyr::filter(
      season == "Dry"
    ),
  
  ggplot2::aes(
    x = x,
    y = y
  ),
  
  width =
    grid_spacing,
  
  height =
    grid_spacing,
  
  fill =
    "#f1e2cc",
  
  color =
    "grey75",
  
  linewidth =
    0.08
) +
  
  
##########################################################
# Camada superior: manchas-chave
##########################################################

ggplot2::geom_tile(
  
  data = activity_overlap_cells %>%
    dplyr::filter(
      season == "Dry"
    ),
  
  ggplot2::aes(
    x = x,
    y = y,
    fill = overlap_category
  ),
  
  width =
    grid_spacing,
  
  height =
    grid_spacing,
  
  color =
    "black",
  
  linewidth =
    0.15
) +
  
  
##########################################################
# Cores das categorias
##########################################################

ggplot2::scale_fill_manual(
  
  values = c(
    
    "Roosting only" =
      "#e7298a",
    
    "Feeding only" =
      "#b2df8a",
    
    "Roosting and feeding" =
      "#ffff33"
  ),
  
  drop = FALSE
) +
  
  
##########################################################
# Extensão espacial
##########################################################

ggplot2::coord_equal(
  
  xlim = c(432000, 664000),
  ylim = c(8752000, 9046000),
  
  expand = FALSE
) +
  
  
  ggplot2::labs(
    x = "Longitude (m)",
    y = "",
    fill = "Key-area category"
  ) +
  
  
  ggplot2::theme_bw()


############################################################
# MAPA — RAINY SEASON
############################################################

p_overlap_activity_rainy <- ggplot2::ggplot() +
  
##########################################################
# Camada de fundo: todas as manchas usadas na chuva
##########################################################

ggplot2::geom_tile(
  
  data = all_used_habitat_patches %>%
    dplyr::filter(
      season == "Rainy"
    ),
  
  ggplot2::aes(
    x = x,
    y = y
  ),
  
  width =
    grid_spacing,
  
  height =
    grid_spacing,
  
  fill =
    "#f1e2cc",
  
  color =
    "grey75",
  
  linewidth =
    0.08
) +
  
  
##########################################################
# Camada superior: manchas-chave
##########################################################

ggplot2::geom_tile(
  
  data = activity_overlap_cells %>%
    dplyr::filter(
      season == "Rainy"
    ),
  
  ggplot2::aes(
    x = x,
    y = y,
    fill = overlap_category
  ),
  
  width =
    grid_spacing,
  
  height =
    grid_spacing,
  
  color =
    "black",
  
  linewidth =
    0.15
) +
  
  
##########################################################
# Cores das categorias
##########################################################

ggplot2::scale_fill_manual(
  
  values = c(
    
    "Roosting only" =
      "#e7298a",
    
    "Feeding only" =
      "#b2df8a",
    
    "Roosting and feeding" =
      "#ffff33"
  ),
  
  drop = FALSE
) +
  
  
##########################################################
# Extensão espacial
##########################################################

ggplot2::coord_equal(
  
  xlim = c(432000, 664000),
  ylim = c(8752000, 9046000),
  
  expand = FALSE
) +
  
  
  ggplot2::labs(
    x = "Longitude (m)",
    y = "Latitude (m)",
    fill = "Key-area category"
  ) +
  
  
  ggplot2::theme_bw()



figure_activity_overlap <- ggpubr::ggarrange(
  
  p_overlap_activity_rainy,
  p_overlap_activity_dry,
  
  nrow = 1,
  ncol = 2,
  
  labels = c(
    "a",
    "b"
  ),
  
  common.legend = TRUE,
  
  legend = "right"
)


figure_activity_overlap

############################################################
# JUNTAR OS MAPAS FINAIS
############################################################

figure_activity_overlap <- ggpubr::ggarrange(
  
  p_overlap_activity_rainy,
  p_overlap_activity_dry,
  
  nrow = 1,
  ncol = 2,
  
  labels = c(
    "a",
    "b"
  ),
  
  common.legend = TRUE,
  
  legend = "right"
)


figure_activity_overlap


############################################################
# TABELA COMPLETA DE CATEGORIAS PARA O MAPA
############################################################

activity_overlap_map <- all_used_habitat_patches %>%
  
  dplyr::left_join(
    
    activity_overlap_cells %>%
      dplyr::select(
        season,
        cell_id,
        overlap_category
      ),
    
    by = c(
      "season",
      "cell_id"
    )
  ) %>%
  
  dplyr::mutate(
    
    overlap_category =
      dplyr::coalesce(
        as.character(
          overlap_category
        ),
        "Other used patches"
      ),
    
    overlap_category =
      factor(
        overlap_category,
        
        levels = c(
          "Other used patches",
          "Roosting only",
          "Feeding only",
          "Roosting and feeding"
        )
      )
  )


## Mapa mais simples usando uma única camada

### DRY
p_overlap_activity_dry <- ggplot2::ggplot(
  
  activity_overlap_map %>%
    dplyr::filter(
      season == "Dry"
    ),
  
  ggplot2::aes(
    x = x,
    y = y,
    fill = overlap_category
  )
) +
  
  ggplot2::geom_tile(
    width = grid_spacing,
    height = grid_spacing,
    color = "grey65",
    linewidth = 0.08
  ) +
  
  ggplot2::scale_fill_manual(
    
    values = c(
      
      "Other used patches" =
        "#f1e2cc",
      
      "Roosting only" =
        "#e7298a",
      
      "Feeding only" =
        "#b2df8a",
      
      "Roosting and feeding" =
        "#ffff33"
    ),
    
    drop = FALSE
  ) +
  
  ggplot2::coord_equal(
    
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    
    expand = FALSE
  ) +
  
  ggplot2::labs(
    x = "Longitude (m)",
    y = "Latitude (m)",
    fill = "Habitat-patch category"
  ) +
  
  ggplot2::theme_bw()

### RAINY
p_overlap_activity_rainy <- ggplot2::ggplot(
  
  activity_overlap_map %>%
    dplyr::filter(
      season == "Rainy"
    ),
  
  ggplot2::aes(
    x = x,
    y = y,
    fill = overlap_category
  )
) +
  
  ggplot2::geom_tile(
    width = grid_spacing,
    height = grid_spacing,
    color = "grey65",
    linewidth = 0.08
  ) +
  
  ggplot2::scale_fill_manual(
    
    values = c(
      
      "Other used patches" =
        "#f1e2cc",
      
      "Roosting only" =
        "#e7298a",
      
      "Feeding only" =
        "#b2df8a",
      
      "Roosting and feeding" =
        "#ffff33"
    ),
    
    drop = FALSE
  ) +
  
  ggplot2::coord_equal(
    
    xlim = c(432000, 664000),
    ylim = c(8752000, 9046000),
    
    expand = FALSE
  ) +
  
  ggplot2::labs(
    x = "Longitude (m)",
    y = "",
    fill = "Habitat-patch category"
  ) +
  
  ggplot2::theme_bw()




############################################################
# JACCARD SAZONAL PARA CADA ATIVIDADE
############################################################

jaccard_activity_between_seasons <- data.frame(
  
  activity = c(
    "Roosting",
    "Feeding"
  ),
  
  
  n_dry_key = c(
    length(
      dry_roosting_key
    ),
    
    length(
      dry_feeding_key
    )
  ),
  
  
  n_rainy_key = c(
    length(
      rainy_roosting_key
    ),
    
    length(
      rainy_feeding_key
    )
  ),
  
  
  intersection = c(
    
    length(
      intersect(
        dry_roosting_key,
        rainy_roosting_key
      )
    ),
    
    length(
      intersect(
        dry_feeding_key,
        rainy_feeding_key
      )
    )
  ),
  
  
  union = c(
    
    length(
      union(
        dry_roosting_key,
        rainy_roosting_key
      )
    ),
    
    length(
      union(
        dry_feeding_key,
        rainy_feeding_key
      )
    )
  )
) %>%
  
  dplyr::mutate(
    
    jaccard =
      dplyr::if_else(
        union > 0,
        intersection / union,
        NA_real_
      ),
    
    turnover =
      1 - jaccard
  )


jaccard_activity_between_seasons
## Results

# jaccard_activity_between_seasons
#   activity  n_dry_key  n_rainy_key  intersection  union   jaccard    turnover
# 1 Roosting      9          38            6         41    0.1463415  0.8536585
# 2 Feeding       50         164           33        181   0.1823204  0.8176796



write.csv2(jaccard_activity_between_seasons, "path/4.Jaccard_activity_between_seasons_10%_turkey.csv", row.names = FALSE)


jaccard_roosting_feeding # Dentro de cada estação, qual é a sobreposição entre as principais manchas de roosting e feeding?

jaccard_activity_between_seasons # Para cada atividade, qual é a sobreposição das manchas-chave entre seca e chuva?




