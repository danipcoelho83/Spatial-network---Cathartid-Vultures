
##################### PAPER 3 #############################
#################### Turkey vultures - RASO ######################
#################### Spatial Network ######################

# Pos doc project - Spatial Networks
# Date: 11/08/26
# Author: Daniela P. Coelho
# Aim: Calculating graph level metrics (network size, connectance, global clustering (transitivity), and degree distribution according to seasonality (dry season and rainy season) for Turkey vultures (n = 25) from RASO using node size  = 4 km2 (2000x2000m).


########################## RUN "0.Base_code.R" #########################

library(moveNT)
library(igraph)
library(dplyr)
library(raster)


############################################################
# ESTRUTURA DA MATRIZ — DRY SEASON
############################################################

mat_dry <- adj_dry[[1]]
range(mat_dry) # 0 - 1376
dim(mat_dry) # 532 x 532
length(mat_dry) # 283024

# Número de entradas não nulas, incluindo a diagonal
(nonzero_entries_dry <- sum(mat_dry > 0)) # 2215

# Número de entradas iguais a zero
(zero_entries_dry <- sum(mat_dry == 0)) # 280809

# Número total de transições observadas
(total_transitions_dry <- sum(mat_dry)) # 13828

# Número total de eventos de self-loop
(self_loop_events_dry <- sum(diag(mat_dry))) # 7226

# Número de manchas com pelo menos um self-loop
(self_loop_nodes_dry <- sum(diag(mat_dry) > 0)) # 216

matrix_diagnostics_dry <- data.frame(
  season = "Dry",

  minimum_weight =
  min(mat_dry),
  
  maximum_weight =
    max(mat_dry),
  
  matrix_rows =
    nrow(mat_dry),
  
  matrix_columns =
    ncol(mat_dry),
  
  matrix_entries =
    length(mat_dry),
  
  nonzero_entries =
    nonzero_entries_dry,
  
  zero_entries =
    zero_entries_dry,
  
  total_transitions =
    total_transitions_dry,
  
  self_loop_events =
    self_loop_events_dry,
  
  nodes_with_self_loop =
    self_loop_nodes_dry
)

matrix_diagnostics_dry


############################################################
# ESTRUTURA DA MATRIZ — RAINY SEASON
############################################################

mat_rainy <- adj_rainy[[1]]
range(mat_rainy) # 0 - 2668
dim(mat_rainy) # 1686  x 1686 
length(mat_rainy) # 2842596


(nonzero_entries_rainy <- sum(mat_rainy > 0)) # 8555
(zero_entries_rainy <- sum(mat_rainy == 0)) # 2834041
(total_transitions_rainy <- sum(mat_rainy)) # 40951
(self_loop_events_rainy <- sum(diag(mat_rainy))) # 22479
(self_loop_nodes_rainy <- sum(diag(mat_rainy) > 0)) # 757


matrix_diagnostics_rainy <- data.frame(
  
  season = "Rainy",
  
  minimum_weight =
    min(mat_rainy),
  
  maximum_weight =
    max(mat_rainy),
  
  matrix_rows =
    nrow(mat_rainy),
  
  matrix_columns =
    ncol(mat_rainy),
  
  matrix_entries =
    length(mat_rainy),
  
  nonzero_entries =
    nonzero_entries_rainy,
  
  zero_entries =
    zero_entries_rainy,
  
  total_transitions =
    total_transitions_rainy,
  
  self_loop_events =
    self_loop_events_rainy,
  
  nodes_with_self_loop =
    self_loop_nodes_rainy
)

matrix_diagnostics_rainy


## Juntando os diagnósticos

matrix_diagnostics <- bind_rows(
  matrix_diagnostics_dry,
  matrix_diagnostics_rainy
)

matrix_diagnostics

write.csv2(matrix_diagnostics,"path/Adj_matrix_diagnostics_turkey.csv", row.names = FALSE)



############################################################
# MATRIZ BINÁRIA SEM SELF-LOOPS — DRY
############################################################

mat_dry_no_loop <- mat_dry

# Remover movimentos que permaneceram na mesma célula
diag(mat_dry_no_loop) <- 0

# Converter frequência em presença ou ausência de link
mat_dry_binary <- 1L * (mat_dry_no_loop > 0)


############################################################
# MATRIZ BINÁRIA SEM SELF-LOOPS — RAINY
############################################################

mat_rainy_no_loop <- mat_rainy

diag(mat_rainy_no_loop) <- 0

mat_rainy_binary <- 1L * (mat_rainy_no_loop > 0)


############################################################
# GRAFOS BINÁRIOS DIRECIONADOS
############################################################

g_dry_directed <- igraph::graph_from_adjacency_matrix(
  mat_dry_binary,
  mode = "directed",
  weighted = NULL,
  diag = FALSE
)


g_rainy_directed <- igraph::graph_from_adjacency_matrix(
  mat_rainy_binary,
  mode = "directed",
  weighted = NULL,
  diag = FALSE
)



## Criar grafos não direcionados para o clustering global

############################################################
# GRAFOS NÃO DIRECIONADOS PARA CLUSTERING
############################################################

g_dry_undirected <- igraph::as_undirected(
  g_dry_directed,
  mode = "collapse"
)

g_dry_undirected <- igraph::simplify(
  g_dry_undirected,
  remove.multiple = TRUE,
  remove.loops = TRUE
)


g_rainy_undirected <- igraph::as_undirected(
  g_rainy_directed,
  mode = "collapse"
)

g_rainy_undirected <- igraph::simplify(
  g_rainy_undirected,
  remove.multiple = TRUE,
  remove.loops = TRUE
)


############################################################
# NETWORK SIZE
############################################################

(network_size_dry <- igraph::vcount(
  g_dry_directed)) # 532

(network_size_rainy <- igraph::vcount(
  g_rainy_directed)) # 1686



############################################################
# CONNECTANCE — DRY SEASON
############################################################

# Número de links direcionados observados
(n_links_dry <- igraph::ecount(
  g_dry_directed)) # 1999


# Número máximo possível de links direcionados sem loops
(possible_links_dry <-
  network_size_dry *
  (network_size_dry - 1)) # 282492


(connectance_dry <-
  igraph::edge_density(
    g_dry_directed,
    loops = FALSE)) # 0.007076307



############################################################
# CONNECTANCE — RAINY SEASON
############################################################

(n_links_rainy <- igraph::ecount(
  g_rainy_directed)) # 7798

(possible_links_rainy <-
  network_size_rainy *
  (network_size_rainy - 1)) # 2840910


(connectance_rainy <-
  igraph::edge_density(
    g_rainy_directed,
    loops = FALSE)) # 0.002744895


############################################################
# GLOBAL CLUSTERING COEFFICIENT
############################################################

(global_clustering_dry <-
  igraph::transitivity(
    g_dry_undirected,
    type = "global" )) # 0.3919406


(global_clustering_rainy <-
  igraph::transitivity(
    g_rainy_undirected,
    type = "global" )) # 0.3914546




############################################################
# NÓS ISOLADOS — CONTROLE DE QUALIDADE
############################################################

## O código-base cria um nó para toda célula ocupada. Entretanto, uma célula pode ter sido representada somente por um registro pertencente a um burst com um único ponto e, portanto, não possuir nenhuma transição válida. É útil conferir:

(isolated_nodes_dry <- sum(
  igraph::degree(
    g_dry_directed,
    mode = "all",
    loops = FALSE) == 0)) # 24


(isolated_nodes_rainy <- sum(
  igraph::degree(
    g_rainy_directed,
    mode = "all",
    loops = FALSE) == 0)) # 5


## Calculando o número de nós que efetivamente participam de pelo menos um link na estação chuvosa:

(connected_nodes_dry <-
    network_size_dry - isolated_nodes_dry) # 508

(connected_nodes_rainy <-
  network_size_rainy - isolated_nodes_rainy) # 1681



### Embora não esteja entre as três métricas principais, esta verificação ajuda a avaliar se a rede está muito fragmentada.

############################################################
# GIANT COMPONENT — DIAGNÓSTICO COMPLEMENTAR
############################################################

(components_dry <- igraph::components(
  g_dry_undirected)) # 37
# 1 Comp = 479 nodes
# 2 Comps = 4 nodes
# 1 Comp = 3 nodes
# 9 Comps = 2 nodes
# 24 Comps = 1 node

(components_rainy <- igraph::components(
  g_rainy_undirected)) # 7
# 1 Comp = 1679 nodes
# 1 Comp = 2 nodes
# 5 Comps 3 = 1 node

(giant_component_dry <- max(
  components_dry$csize)) # 479

(giant_component_rainy <- max(
  components_rainy$csize)) # 1679


(giant_component_prop_dry <-
  giant_component_dry /
  network_size_dry) # 0.9003759

(giant_component_prop_rainy <-
  giant_component_rainy /
  network_size_rainy) # 0.9958482



############################################################
# FINAL NETWORK-LEVEL METRICS - TABLE
############################################################

(network_metrics <- data.frame(
  
  season = c(
    "Dry",
    "Rainy"
  ),
  
  network_size_nodes = c(
    network_size_dry,
    network_size_rainy
  ),
  
  directed_links = c(
    n_links_dry,
    n_links_rainy
  ),
  
  possible_directed_links = c(
    possible_links_dry,
    possible_links_rainy
  ),
  
  connectance = c(
    connectance_dry,
    connectance_rainy
  ),
  
  global_clustering = c(
    global_clustering_dry,
    global_clustering_rainy
  ),
  
  isolated_nodes = c(
    isolated_nodes_dry,
    isolated_nodes_rainy
  ),
  
  connected_nodes = c(
    connected_nodes_dry,
    connected_nodes_rainy
  ),
  
  giant_component_nodes = c(
    giant_component_dry,
    giant_component_rainy
  ),
  
  giant_component_proportion = c(
    giant_component_prop_dry,
    giant_component_prop_rainy
  ),
  
  total_transitions = c(
    total_transitions_dry,
    total_transitions_rainy
  )
))

write.csv2(network_metrics,"path/Network_level_metrics_turkey.csv", row.names = FALSE)



### Plotting Degree distribution
library(ggplot2)
library(ggpubr)

############################################################
# CALCULAR O GRAU — DRY SEASON
############################################################

degree_in_dry <- igraph::degree(
  g_dry_directed,
  mode = "in",
  loops = FALSE
)

degree_out_dry <- igraph::degree(
  g_dry_directed,
  mode = "out",
  loops = FALSE
)

degree_total_dry <- igraph::degree(
  g_dry_directed,
  mode = "all",
  loops = FALSE
)


############################################################
# CALCULAR O GRAU — RAINY SEASON
############################################################

degree_in_rainy <- igraph::degree(
  g_rainy_directed,
  mode = "in",
  loops = FALSE
)

degree_out_rainy <- igraph::degree(
  g_rainy_directed,
  mode = "out",
  loops = FALSE
)

degree_total_rainy <- igraph::degree(
  g_rainy_directed,
  mode = "all",
  loops = FALSE
)

## Conferência
all(degree_total_dry == degree_in_dry + degree_out_dry)

all(degree_total_rainy == degree_in_rainy + degree_out_rainy)


############################################################
# DATA FRAME — DRY SEASON
############################################################

df_degree_dry <- data.frame(
  node_id = seq_along(degree_total_dry),
  
  In_degree = as.numeric(degree_in_dry),
  
  Out_degree = as.numeric(degree_out_dry),
  
  Degree = as.numeric(degree_total_dry)
)


############################################################
# DATA FRAME — RAINY SEASON
############################################################

df_degree_rainy <- data.frame(
  node_id = seq_along(degree_total_rainy),
  
  In_degree = as.numeric(degree_in_rainy),
  
  Out_degree = as.numeric(degree_out_rainy),
  
  Degree = as.numeric(degree_total_rainy)
)


############################################################
# RESUMO DO GRAU
############################################################

range(df_degree_dry$Degree) # 0-70

range(df_degree_rainy$Degree) # 0-121


# Número de nós com grau zero
sum(df_degree_dry$Degree == 0) # 24

sum(df_degree_rainy$Degree == 0) # 5


# Grau médio
mean(df_degree_dry$Degree) # 7.515038
 
mean(df_degree_rainy$Degree) # 9.250297

# Grau mediano
median(df_degree_dry$Degree) # 4

median(df_degree_rainy$Degree) # 2


############################################################
# LIMITES COMUNS DOS EIXOS
############################################################

maximum_degree <- max(
  df_degree_dry$Degree,
  df_degree_rainy$Degree,
  na.rm = TRUE
)


maximum_frequency <- max(
  as.numeric(table(df_degree_dry$Degree)),
  as.numeric(table(df_degree_rainy$Degree)),
  na.rm = TRUE
)


y_upper_limit <- ceiling(
  maximum_frequency * 1.05
)


############################################################
# DEGREE DISTRIBUTION — DRY SEASON
############################################################

(G1 <- ggplot(
  df_degree_dry,
  aes(x = Degree)
) +
  
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    fill = "orange",
    color = "black",
    linewidth = 0.4
  ) +
  
  labs(
    title = "",
    x = "Degree",
    y = ""
  ) +
  
  scale_x_continuous(
    limits = c(-0.5, maximum_degree + 0.5),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  
  coord_cartesian(
    ylim = c(0, y_upper_limit)
  ) +
  
  theme_minimal(base_size = 16) +
  
  theme(
    plot.title = element_text(hjust = 0.5),
    
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(size = 16),
    
    axis.text = element_text(size = 14)
  ))



############################################################
# DEGREE DISTRIBUTION — RAINY SEASON
############################################################

(G2 <- ggplot(
  df_degree_rainy,
  aes(x = Degree)
) +
  
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    fill = "skyblue",
    color = "black",
    linewidth = 0.4
  ) +
  
  labs(
    title = "",
    x = "Degree",
    y = "Habitat patches"
  ) +
  
  scale_x_continuous(
    limits = c(-0.5, maximum_degree + 0.5),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  
  coord_cartesian(
    ylim = c(0, y_upper_limit)
  ) +
  
  theme_minimal(base_size = 16) +
  
  theme(
    plot.title = element_text(hjust = 0.5),
    
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(size = 16),
    
    axis.text = element_text(size = 14)
  ))



ggarrange(G2, G1, nrow = 1, ncol = 2, legend = "right", labels = c("a", "b"))









