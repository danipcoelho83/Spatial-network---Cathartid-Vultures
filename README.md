# Spatial interaction networks of the Turkey Vultures (*Cathartes aura*) on the Caatinga biome

This repository contains the R scripts used to investigate seasonal variation in habitat use by individual Turkey Vultures inhabiting a seasonally dry tropical forest (Caatinga biome) in northern Bahia state, Brazil. The results of these analyses are being used to prepare a manuscript titled "Climate seasonality shapes the use of space by a species of vulture, promoting modifications in habitat connectivity in a seasonally dry tropical forest."

The study employs an interaction network approach to investigate the effects of seasonality on habitat use by individuals monitored via GPS devices (biologgers). We used network-level metrics (size, connectance, and modularity) to investigate the spatial organization of habitat use between the dry and rainy seasons. Additionally, we used land-cover data from MapBiomas to investigate how the vultures promote landscape-scale connectivity among different habitat types. Finally, we used node-level metrics (degree centrality, self-loop centrality, betweenness centrality, and clustering coefficient) to examine variations in preferentially used areas across seasons.


## Project status

**This repository is currently under development.**

Data analyses for the associated scientific manuscript are ongoing. Some scripts contain completed data-processing and descriptive analytical workflows, whereas other analyses are exploratory or await additional procedures to account for differences in sampling effort between seasons.

Results and analytical decisions may therefore change as the study develops.

The repository should not yet be considered the final reproducible archive of the associated manuscript.

## Study overview

GPS tracking data are used to construct seasonal spatial movement networks in which:

* **nodes** represent 2 × 2 km spatial cells (4 km² habitat patches);
* **edges** represent observed movements between consecutive occupied cells;
* **edge direction** preserves the observed direction of movement;
* **edge weights** represent the number of observed transitions between pairs of cells; and
* **self-loops** represent consecutive observations occurring within the same spatial cell.

A single common spatial grid is used for both seasonal networks. Therefore, the permanent spatial identifier `cell_id` always represents the same geographic location in the dry- and rainy-season networks.

The biologgers capture geographic coordinates at intervals ranging from 5 to 10 minutes per individual. To standardize data collection, we filtered the dataset to obtain one coordinate every 30 minutes.

Seasonal classification is based on local time in Bahia, Brazil:

* **Rainy season:** January–May
* **Dry season:** June–December

Spatial coordinates are projected to UTM Zone 24S (EPSG:32724).

## Repository structure

```text
.
├── README.md
├── .gitignore
│
├── R/
│   ├── 0.Base_code.R
│   ├── 1.Movement_analysis.R
│   ├── 2.Graph_level_metrics.R
│   └── 4.Node_level_metrics.R
│
├── data/
│   └── README.md
│
├── results/
│
└── figures/
```

## Analysis workflow

### `0.Base_code.R`

Main data-processing and spatial-network construction script.

The script:

1. retrieves GPS tracking data from Movebank;
2. removes duplicated timestamps;
3. standardizes GPS records to 30-minute intervals;
4. converts timestamps to local Bahia time;
5. classifies observations into dry and rainy seasons;
6. projects coordinates to UTM Zone 24S;
7. constructs a common 2 × 2 km spatial grid covering the study area;
8. assigns each grid cell a permanent spatial identifier (`cell_id`);
9. divides trajectories into independent temporal bursts based on gaps between consecutive observations;
10. creates observed transitions between consecutive locations within the same individual and burst;
11. constructs directed seasonal adjacency matrices;
12. generates seasonal node lookup tables linking matrix positions (`node_id`) to permanent spatial cells (`cell_id`); and
13. checks spatial correspondence and consistency between the dry- and rainy-season networks.

This script must be run before the downstream analyses.

### `1.Movement_analysis.R`

Performs exploratory movement-data summaries and visualization.

The current workflow includes:

* number of GPS locations per month;
* distribution of sampling effort among individuals;
* seasonal and individual sampling summaries;
* step-length calculations;
* minimum, maximum, mean, and median step length per individual;
* graphical comparison of step-length distributions;
* spatial visualization of GPS locations; and
* summaries and plots of individual movement trajectories.

This script is currently under development. Some movement calculations and plotting sections require further revision before they are considered part of the final analytical workflow.

### `2.Graph_level_metrics.R`

Calculates descriptive properties of the observed dry- and rainy-season spatial networks.

The script:

1. performs diagnostic checks on the seasonal adjacency matrices;
2. removes self-loops for analyses of inter-cell connectivity;
3. converts transition frequencies into binary directed networks;
4. constructs undirected versions of the networks for clustering analyses;
5. calculates:

   * network size;
   * number of directed links;
   * connectance;
   * global clustering coefficient;
   * number of isolated nodes;
   * number of connected nodes;
   * size of the giant component; and
   * proportion of nodes belonging to the giant component;
6. calculates in-degree, out-degree, and total degree; and
7. visualizes seasonal degree distributions.

These metrics currently describe the observed networks. Seasonal differences should not be interpreted as formal statistical differences until differences in sampling effort between seasons are explicitly accounted for.

### `4.Node_level_metrics.R`

Calculates spatial metrics for individual network nodes and explores the seasonal distribution of key habitat patches.

The script calculates:

* self-loop count;
* contribution of each node to the total number of transitions;
* local self-loop rate;
* in-degree;
* out-degree;
* total degree;
* normalized degree;
* normalized unweighted betweenness centrality;
* betweenness centrality based on the inverse of transition frequency; and
* local clustering coefficient.

Results calculated independently with `igraph` are also compared with metrics returned by `moveNT`.

The script additionally:

1. associates node-level metrics with permanent spatial coordinates;
2. compares the use and centrality of corresponding habitat patches between seasons;
3. produces spatial maps of node-level metrics;
4. classifies highly central habitat patches using the upper 10% of each metric distribution within each season;
5. calculates descriptive Jaccard overlap between key habitat patches in the dry and rainy seasons;
6. classifies transitions according to local-time activity periods;
7. calculates activity-specific self-loops for roosting and feeding periods;
8. identifies the upper 10% of activity-specific self-loop habitat patches;
9. calculates overlap between key roosting and feeding patches within each season;
10. calculates seasonal overlap of key patches separately for roosting and feeding; and
11. maps patches classified as important for roosting, feeding, or both activities.

The Jaccard analyses are currently treated as descriptive summaries. Additional analyses are required before formal seasonal inference, particularly because sampling effort differs between the dry and rainy seasons.

## Activity-specific analyses

Transitions are currently classified according to local Bahia time into behavioral periods including:

* **Roosting:** early-morning and evening periods;
* **Feeding:** daytime period; and
* **Other:** observations outside the focal time windows.

A transition is assigned to an activity only when both its starting and ending GPS records fall within the same activity period.

Activity-specific self-loop analyses are then used to identify habitat patches with high local use during roosting and feeding periods.

## Current analytical development

The following components are already implemented:

* common spatial grid construction;
* permanent spatial node identification;
* seasonal network construction;
* graph-level descriptive metrics;
* node-level centrality metrics;
* activity-specific self-loop analyses;
* spatial mapping of central habitat patches; and
* descriptive overlap analyses.

The following components remain under development or require additional validation:

* final movement and step-length analyses;
* formal treatment of unequal sampling effort between seasons;
* statistical inference for seasonal differences in network metrics;
* validation of thresholds used to classify key habitat patches; and
* final integration of network results with the hypotheses of the associated manuscript.

The analytical workflow and this README will be updated as these components are completed.

## Main R dependencies

The current scripts use packages including:

`move`, `move2`, `moveNT`, `adehabitatLT`, `sf`, `terra`, `raster`, `sp`, `lubridate`, `dplyr`, `tidyr`, `igraph`, `ggplot2`, `ggpubr`, `geosphere`, and related dependencies.

Individual scripts specify the packages required for each analysis.

## Data availability

Raw GPS tracking data are not distributed through this repository.

Movement records used by the project originate from the Movebank study:

**"Cathartid ecology - Dénes - Brazil"**

Access to the original movement data should follow the permissions and data-access conditions established for the Movebank study.

Users wishing to reproduce the analyses should obtain appropriate access to the underlying tracking data.

## Reproducibility

Before running the analyses:

1. configure Movebank credentials securely;
2. replace machine-specific input and output paths with appropriate local or project-relative paths;
3. install the required R packages; and
4. run `0.Base_code.R` before the downstream scripts.

Passwords, API keys, access tokens, confidential tracking data, and machine-specific configuration files should never be committed to the repository.

## Suggested execution order

```text
0.Base_code.R
       │
       ├── 1.Movement_analysis.R
       ├── 2.Graph_level_metrics.R
       └── 4.Node_level_metrics.R
```

The downstream scripts depend on objects created by `0.Base_code.R`, including seasonal movement data, the common spatial grid, adjacency matrices, transition tables, and node lookup tables.

## Author

**Daniela Pinto-Coelho**

Postdoctoral research project on movement ecology and spatial networks in Raso da Catarina, Bahia, Brazil.

## Citation

Citation information for the associated manuscript will be added when available.

Because the analyses are currently under development, users are encouraged to contact the author before using preliminary analytical components.

## License

A software license has not yet been specified.
