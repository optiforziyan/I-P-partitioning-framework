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
Anderson, M., Gorley, R., Clarke, K. (2008). PERMANOVA+ for PRIMER: guide to software and statistical methods. Plymouth, United Kingdom: PRIMER-E Ltd.
Brown, J. H., Mehlman, D. W., Stevens, G. C. (1995). Spatial variation in abundance. Ecology, 76(7), 2028-2043. 
Cadotte, M. W., Jonathan Davies, T., Regetz, J., Kembel, S. W., Cleland, E., & Oakley, T. H. (2010). Phylogenetic diversity metrics for ecological communities: integrating species richness, abundance and evolutionary history. Ecology letters, 13(1), 96-105.
Cappelatti, L., Mauffrey, A. R., Griffin, J. N. (2020). Functional diversity of habitat formers declines scale-dependently across an environmental stress gradient. Oecologia, 194(1-2), 135-149.
Cavender-Bares, J., Kozak, K. H., Fine, P. V., & Kembel, S. W. (2009). The merging of community ecology and phylogenetic biology. Ecology letters, 12(7), 693-715.
Chen Y., Wu Y., Zhou J., Zhang W., Lin H., Liu X., Pan K., Shen T., Pan Z. (2021a). Effectively inferring overall spatial distribution pattern of species in a map when exact coordinate information is missing. Methods in Ecology and Evolution. 12, 971-984.
Chen, Y, Shen, T-J, Condit, R. (2021b). Disentangling multi-species aggregate versus overlapping distributions. Journal of Vegetation Science, 32, e12973.
Chen, Y. (2015). Biodiversity and biogeographic patterns in Asia-Pacific Region I: statistical methods and case studies. Bentham Science Publishers.
Chen, Y., Shen, T. J., Condit, R., Hubbell, S.P. (2018). Community-level species’ correlated distribution can be scale-independent and related to the evenness of abundance. Ecology, 99, 2787-2800.
Chen, Y., Shen, T. J., Van Chung, H., Shi, S., Jiang, J., Condit, R., Hubbell, S. P. (2019). Inferring multispecies distributional aggregation level from limited line transect‐derived biodiversity data. Methods in Ecology and Evolution, 10(7), 1015-1023.
Chen, Y., Shen, T.-J. (2024). Distributional ecology: Opening new research windows by addressing aggregation-related puzzles. Diversity and Distributions, 30, e13818. 
Condit, R., Perez, R., Aguilar, S., Lao, S., Foster, R. and Hubbell, S. (2019a). BCI 50-ha plot taxonomy. DOI: 10.15146/R3FH61.
Condit, R., Perez, R., Aguilar, S., Lao, S., Foster, R. and Hubbell, S. (2019b). Complete data from the Barro Colorado 50-ha plot: 423617 trees, 35 years, 2019 version. DOI: 10.15146/5xcp-0d46.
Conlisk, E., Conlisk, J. and Harte, J. (2007). The impossibility of estimating a negative binomial clustering parameter from presence-absence data: a comment on He and Gaston. American Naturalist, 170, 651-654.
Dale M. (1999). Spatial pattern analysis in plant ecology. Cambridge University Press, Cambridge.
Daleo, P., Alberti, J., Chaneton, E. J., Iribarne, O., Tognetti, P. M., Bakker, J. D., ... & Hautier, Y. (2023). Environmental heterogeneity modulates the effect of plant diversity on the spatial variability of grassland biomass. Nature communications, 14(1), 1809.
Darwin J. (1957). The power of the Poisson index of dispersion. Biometrika. 44:286-289.
Döring, T. F., Knapp, S., & Cohen, J. E. (2015). Taylor’s power law and the stability of crop yields. Field Crops Research, 183, 294-302. 
Dos Santos D., Fernandez H., Cuezzo M., Dominguez E. (2008). Sympatry inference and network analysis in biogeography. Systematic Biology. 57, 432-448.
Fine, P. V., & Kembel, S. W. (2011). Phylogenetic community structure and phylogenetic turnover across space and edaphic gradients in western Amazonian tree communities. Ecography, 34(4), 552-565.
Hanski I., Tiainen J. (1989). Bird ecology and Taylor’s variance-mean regression. Annales Zoologici Fennici. 26, 213–217.
Helmus, M. R., Keller, W., Paterson, M. J., Yan, N. D., Cannon, C. H., & Rusak, J. A. (2010). Communities contain closely related species during ecosystem disturbance. Ecology letters, 13(2), 162-174.
Jin, Y., & Qian, H. (2022). V. PhyloMaker2: An updated and enlarged R package that can generate very large phylogenies for vascular plants. Plant Diversity, 44(4), 335-339. 
Kathirgamatamby N. (1953). Note on the Poisson index of dispersion. Biometrika. 40, 225–228.
Keil, P., Wiegand, T., Tóth, A. B., McGlinn, D. J., & Chase, J. M. (2021). Measurement and analysis of interspecific spatial associations as a facet of biodiversity. Ecological Monographs, 91(3), e01452.
Kembel, S. W., Cowan, P. D., Helmus, M. R., Cornwell, W. K., Morlon, H., Ackerly, D. D., ... & Webb, C. O. (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics, 26(11), 1463-1464.
Krebs C. 1998. Ecological Methodology. Boston, MA, USA: Pearson Education. Kremen, C., Cameron, A., Moilanen, A., Phillips, S. J., Thomas, C. D., Beentje, H., ... & Zjhra, M. L. (2008). Aligning conservation priorities across taxa in Madagascar with high-resolution planning tools. Science, 320(5873), 222-226.
Laliberté, E., & Legendre, P. (2010). A distance‐based framework for measuring functional diversity from multiple traits. Ecology, 91(1), 299-305.
Lamy, T., Wisnoski, N. I., Andrade, R., Castorani, M. C., Compagnoni, A., Lany, N., ... & Sokol, E. R. (2021). The dual nature of metacommunity variability. Oikos, 130(12), 2078-2092.
Liao, Z., Zhou, J., Shen, T-J, Chen, Y. (2023). Inferring single- and multi-species distributional aggregation using quadrat sampling. Ecological Indicators, 156, 111085
Liaw, A., & Wiener, M. 2002. Classification and regression by randomForest. R News, 2, 18-22.
McCullagh P., Nelder J. 1989. Generalized linear models. Chapman and Hall/CRC.
Miguet, P., Jackson, H. B., Jackson, N. D., Martin, A. E., & Fahrig, L. (2016). What determines the spatial extent of landscape effects on species?. Landscape ecology, 31, 1177-1194.
Pielou E. 1977. Mathematical ecology. John Wiley & Sons.
Roberts D. 2017. Distance, dissimilarity, and mean-variance ratios in ordination. Methods in Ecology and Evolution. 8, 1398–1407.
Salazar, D., Jaramillo, M. A., & Marquis, R. J. (2016). Chemical similarity and local community assembly in the species rich tropical genus Piper. Ecology, 97(11), 3176-3183.
Sarkar, D., & Andrews, F. (2013). latticeExtra: extra graphical utilities based on lattice. R package version 0.6-26.
Singh, A., Yadav, A., & Rana, A. (2013). K-means with Three different Distance Metrics. International Journal of Computer Applications, 67(10).
Smets, T., Verbeeck, N., Claesen, M., Asperger, A., Griffioen, G., Tousseyn, T., ... & De Moor, B. (2019). Evaluation of distance metrics and spatial autocorrelation in uniform manifold approximation and projection applied to mass spectrometry imaging data. Analytical chemistry, 91(9), 5706-5714.
Smith, S. A., and Brown, J. W. (2018). Constructing a broadly inclusive seed plant phylogeny. American Journal of Botany, 105, 302-314.
Stiteler W., Patii G. (1969). Variance-to-mean ratio and Morisita’s index as measures of spatial patterns in ecological populations. Statistical Ecology. 1, 423–459.
Suárez‐Castro, A. F., Raymundo, M., Bimler, M., & Mayfield, M. M. (2022). Using multi‐scale spatially explicit frameworks to understand the relationship between functional diversity and species richness. Ecography, 2022(6), e05844.
Szumik C., Pereyra V., Casagranda M. (2019). Areas of endemism: to overlap or not to overlap, that is the question. Cladistics, 35, 198–229.
Taylor, L. R. (1961). Aggregation, variance and the mean. Nature, 189(4766), 732-735.
Tsianou, M. A., & Kallimanis, A. S. (2020). Geographical patterns and environmental drivers of functional diversity and trait space of amphibians of Europe. Ecological Research, 35(1), 123-138.
Wang, S., & Loreau, M. (2014). Ecosystem stability in space: α, β and γ variability. Ecology letters, 17(8), 891-901.
Wang, S., & Loreau, M. (2016). Biodiversity and ecosystem stability across scales in metacommunities. Ecology letters, 19(5), 510-518.
Warton D., Wright S., Wang Y. (2012). Distance-based multivariate analyses confound location and dispersion effects. Methods in Ecology and Evolution. 3:89-101.
Webb, C. O. (2000). Exploring the phylogenetic structure of ecological communities: an example for rain forest trees. The American Naturalist, 156(2), 145-155.
Webb, C. O., Ackerly, D. D. & Kembel, S. W. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 24, 98-100.
Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002). Phylogenies and community ecology. Annual review of ecology and systematics, 33(1), 475-505.
Wickham, H., (2016). Ggplot2: Elegant Graphics for Data 581 Analysis, Use R! Springer-Verlag New York.
Xing, D., & He, F. (2018). Environmental filtering explains a U‐shape latitudinal pattern in regional β‐deviation for eastern North American trees. Ecology letters, 22, 284-291.
Zanne, A. E., Tank, D. C., Cornwell, W. K., Eastman, J. M., Smith, S. A., FitzJohn, R. G., ... & Beaulieu, J. M. (2014). Three keys to the radiation of angiosperms into freezing environments. Nature, 506(7486), 89-92. 

