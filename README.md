# Additive partitioning of multi-species distributional aggregation with empirical applications

## Partitioning framework

Here, we provide a self-explained R code (Alternative additive partitioning frameworks of distributional aggregation.R）for caculating the data dispersion (represented by a Poisson index of dispersion *I*) and spatial proximity (represented by an index as the average of pairwise Euclidean distances *P*) presented in our manuscript. 

###  **Folders**
#### *functions* - Custom R Script
1. Dirichlet-Multinomial distribution-S1.R
* **data_creation()**               - A function to generate a random dataset for spatial ecology simulation
* **species.distribution()**        - Simulation of species distribution using Poisson cluster process
* **I_sp_partitioning()**           - Partitioning abundance-based data dispersion into within and between-species components (Eq.1)
* **I_genus_sp_partitioning()**     - Partitioning of abundance-based data dispersion into within-genus within-species, within-genus between-species and between-genus components (Eq.3)
* **I_site_partitioning_simple()**  - Partitioning abundance dispersion into within and between-site components without using parallel computing (Eq.2)
* **I_site_partitioning_parall()**  - Partitioning abundance dispersion into within and between-site components using parallel computing (Eq.2)
* **I_site_genus_partitioning_1()** - Partitioning of abundance dispersion into within-site within-genus, within-site between-genus and between-site components (Eq.4)
* **I_site_genus_partitioning_2()** - Partitioning of abundance dispersion into within-site, between-site within-genus and between-site between-genus components (Eq.5)
* **D_sp_partitioning()**           - Partitioning pairwise distances into within and between-species components (Eq.6)
* **D_site_partitioning()**         - Partitioning pairwise distances into within-site and between-site components (Eq.6)
* **D_genus_sp_partitioning()**     - Partitioning of pairwise distances into within-genus within-species distance, within-genus between-species distance and between-genus distance components (Eq.7)
* **D_site_genus_partitioning_1()** - Partitioning of pairwise distances into within-site within-genus distance, within-site between-genus distance and between-site distance components (one form of Eq.7)
* **D_site_genus_partitioning_2()** - Partitioning of pairwise distances into within-site within-genus distance, within-site between-genus distance and between-site distance components (another form of Eq.7)

#### *input*
* Data_Figure1.xlsx   - An example 


#### *output*


###  Code Availability
#### Software
1. R version 4.0.2 (2020-06-22)
#### Required R packages
1. MASS version 7.3-51.6
2. spatstat version 1.64-1
3. stringr version 1.4.0
4. readxl version 1.0.15
5. doParallel version 1.0.15
## Citation


## Author(s)

## References: 

