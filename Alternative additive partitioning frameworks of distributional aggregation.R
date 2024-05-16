# Partitioning workflow
# Encoding = CP936
# R version 4.0.2
################################################################################ .
# This is an R function to generate a random dataset for spatial ecology         ----
# ZL - 2023-02-15
# Arguments:
# i:   'S' = the number of species to simulate
# ii:  'R' = the number of sampling sites for each species
# iii: 'G' = S species have been found in G genus
# iv: 'lambda' = the mean abundance of each species at each site using Poisson

# The output called 'sxyz', which contains:
# 'genus' = a numeric identifier for the genus (ranging from 1 to G)
# 'sp' = a numeric identifier for the species (ranging from 1 to S).
# 'x'  = the x-coordinate of the sampling site.
# 'y'  = the y-coordinate of the sampling site.
# 'site' = a numeric identifier for the site (ranging from 1 to R).
# 'abun' = the observed abundance of the species at the site,                    

data_creation <- function(S = 8, G = 3, R = 6, lambda = 5, fix.rd = T)
{
  if(fix.rd){
    set.seed(1226)
  }
  x_raw = matrix(rpois(S * R, lambda = lambda), nrow = S)
  sxyz <- matrix(ncol = 6, nrow = S * R )
  colnames(sxyz) = c("genus","sp","x","y","site", "abun")
  if(fix.rd){
    set.seed(2018)
    # Generate a vector of S elements representing the genus to which each species belongs
    species_genus <- sample(1:G, S, replace = TRUE)
    if(length(unique(species_genus)) != G){
      set.seed(1226)
      species_genus <- sample(1:G, S, replace = TRUE)
    }
    if(length(unique(species_genus)) != G){
      set.seed(520)
      species_genus <- sample(1:G, S, replace = TRUE)
    }
  }else{
    species_genus <- sample(1:G, S, replace = TRUE)
    if(length(unique(species_genus)) != G){
      species_genus <- sample(1:G, S, replace = TRUE)
    }
    if(length(unique(species_genus)) != G){
      species_genus <- sample(1:G, S, replace = TRUE)
    }
  }
  sp = genus = x = y = NULL
  for (i in 1:S) {
    sp = c(sp, rep(i, R))
    genus = c(genus,rep(species_genus[i],R))
  }
  
  site.comb <- list() # Initialize an empty list to store all factor combinations
  # Loop to find all factor combinations
  for (i in 1:(R/2)) {
    if (R %% i == 0) {
      site.comb[[length(site.comb) + 1]] <- c(i, R/i)
    }
  }
  
  # All factor combinations
  x = expand.grid(1:site.comb[[round(length(site.comb)/2)]][1], 1:site.comb[[round(length(site.comb)/2)]][2])[,1]
  y = expand.grid(1:site.comb[[round(length(site.comb)/2)]][1], 1:site.comb[[round(length(site.comb)/2)]][2])[,2]
  
  
  sxyz[,"genus"] = genus
  sxyz[,"sp"] = sp
  sxyz[,"x"] = rep(x,S)
  sxyz[,"y"] = rep(y,S)
  sxyz[,"site"] = rep(1:R,S)
  sxyz[,"abun"] = as.array(x_raw)
  return(sxyz)
}

# EOF                                                                            
################################################################################ ----

################################################################################ .
# partitioning abundance dispersion into within and between-species components   ----
# ZL - 2023-02-15 (Eq.1)
# 'sxyz' is a six-column matrix
# the first column represents genus ID
# the second column represents species ID
# the third and fourth columns represented x, and y coordinates
# the five column represent the numeric identifier for each site
# the final column represent abundance in each location 
# (only used for quadrat sampling and dispersion partitioning)
# The adj.w parameter offers the possibility to adjust the two partitioning components with ecological meaning
# If adj.w is assigned as "not", the output of the function will return the original value, regardless of whether the two partition components are negative.
# If adj.w is assigned as "equal", the bisection method will be used to deal with the negative value.
# If adj.w is assigned as "weight", while ignoring negative values, the remaining components are reallocated according to their weights.

I_sp_partitioning <- function(sxyz, adj.w = c("not","weight"))
{ 
  # total number of genus
  G = length(unique(sxyz[,"genus"]))
  cat("Total number of genus:", G, "\n")
  
  # total number of species
  S = length(unique(sxyz[,"sp"]))
  cat("Total number of species:", S, "\n")
  
  # total number of sampling sites
  R = length(unique(sxyz[,"site"]))
  cat("Total number of sampling sites:", R, "\n")
  
  # total number of individuals of all species and all sites
  x_bullet_bullet <- sum(sxyz[,"abun"])
  cat("Total individuals:", x_bullet_bullet, "\n")
  
  # the total number of individuals of a species i in all the sites
  x_i_bullet <- c(tapply(sxyz[,"abun"], sxyz[,"sp"], sum))
  # cat("Total individuals by species:", x_i_bullet, "\n")
  
  # the total number of individuals of different species found in the sampling site j
  x_bullet_j <- c(tapply(sxyz[,"abun"], sxyz[,"site"], sum))
  # cat("Total individuals by site:", x_bullet_j, "\n")
  
  # Mean total number of individuals across all sampling sites
  x_bar_bullet_bullet <- x_bullet_bullet / R
  
  # Mean total number of individuals of each species across all sampling sites
  x_bar_i_bullet <- x_i_bullet / R
  
  # Mean total number of individuals in each sampling site
  x_bar_bullet_j <- x_bullet_j / S
  
  # Partitioning of multi-species Poisson index of dispersion
  I_total <- (sum(x_bullet_j^2) - R * x_bar_bullet_bullet^2) / (x_bar_bullet_bullet)
  cat("Multi-species Poisson Index of dispersion:", I_total , "\n")
  
  # Calculate the within-species intraspecific abundance variation component
  I_within_species <- sum((sapply(1:S, function(i) sum(sxyz[sxyz[,"sp"] == i, "abun"]^2) - R * x_bar_i_bullet[i]^2)) / x_bar_bullet_bullet)
  cat("Within-species intraspecific abundance variation component:", I_within_species, "\n")
  
  if(I_total - I_within_species > 0) {
    
  # Calculate the between-species intraspecific abundance synchrony component
  I_between_species <- sum((sapply(1:S, function(i) sum(sapply(1:S, 
                                                       function(k) { 
    if (k != i) { sum(sxyz[sxyz[,"sp"] == i, "abun"] * sxyz[sxyz[,"sp"] == k, "abun"]) - R * x_bar_i_bullet[i] * x_bar_i_bullet[k]} else {0}
                                                                   }))) / x_bar_bullet_bullet))
  
  cat("Between-species intraspecific abundance synchrony component:", I_between_species, "\n")
  
  } else {
    I_between_species = I_total - I_within_species
    cat("Between-species intraspecific abundance synchrony component:", I_between_species, "\n")
  }
  
  # Check equality
  if (abs(I_within_species + I_between_species - I_total) < 1e-10) {
    message("Total equals within component plus between component")
  } else {
    message("Total not equal to within component plus between component!!!")
  }
  
  # Adjustment
  if(!adj.w %in% c("not","weight")){
    stop(paste0("The adjust parameter adj.w - ", adj.w, " is not supported yet"))
  }
  
  if(adj.w == "not"){
  Res <- data.frame(I.total=c(I_total,I_total),I.within.species=c(I_within_species,I_within_species),
               I.between.species=c(I_between_species,I_between_species))
  rownames(Res) <- c("original", adj.w)
  return(Res) 
  } 
  if(adj.w == "weight"){
    if(sum(sign(c(I_within_species,I_between_species)))==2){
    Res <- data.frame(I.total=c(I_total,I_total),I.within.sp=c(I_within_species,I_within_species),
                        I.between.sp=c(I_between_species,I_between_species))
    rownames(Res) <- c("original", adj.w)
    return(Res)  
    }
    if(I_within_species < 0){
      I_within_species_aj = 0; I_between_species_aj = I_total
    }
    if(I_between_species < 0){
      I_between_species_aj = 0; I_within_species_aj = I_total
    }
    Res <- data.frame(I.total=c(I_total,I_total),I.within.sp=c(I_within_species,I_within_species_aj),
                      I.between.sp=c(I_between_species,I_between_species_aj))
    rownames(Res) <- c("original", adj.w)
  return(Res)  
  }
}

# EOF
################################################################################ ----


################################################################################ .
# partitioning of abundance dispersion into within-genus within-species,         ----
# within-genus between-species and between-genus components
# ZL - 2023-03-02 (Eq.3)
# 'sxyz' is a six-column matrix
# the first column represents genus ID
# the second column represents species ID
# the third and fourth columns represented x, and y coordinates
# the five column represent the numeric identifier for each site
# the final column represent abundance in each location 
# (only used for quadrat sampling and dispersion partitioning)
# The adj.w parameter offers the possibility to adjust the two partitioning components with ecological meaning
# If adj.w is assigned as "not", the output of the function will return the original value, regardless of whether the two partition components are negative.
# If adj.w is assigned as "equal", the bisection method will be used to deal with the negative value.
# If adj.w is assigned as "weight", while ignoring negative values, the remaining components are reallocated according to their weights.

I_genus_sp_partitioning <- function(sxyz, adj.w = c("not","weight"))
{ 
  options(digits = 8)
  
  # total number of genus
  G = length(unique(sxyz[,"genus"]))
  cat("Total number of genus:", G, "\n")
  
  # total number of species
  S = length(unique(sxyz[,"sp"]))
  cat("Total number of species:", S, "\n")
  
  # total number of sampling sites
  R = length(unique(sxyz[,"site"]))
  cat("Total number of sampling sites:", R, "\n")
  
  # Calculate total number of species per genus
  Sg <- tapply(sxyz[,"sp"], sxyz[,"genus"], FUN = length)/R
  Sgn <- tapply(sxyz[,"sp"], sxyz[,"genus"], FUN = unique)
  cat("Total number of species per genus:", Sg, "\n")
  
  # total number of individuals of all species and all sites
  x_bullet_bullet_bullet <- sum(sxyz[,"abun"])
  cat("Total individuals:", x_bullet_bullet_bullet, "\n")
  
  # the total number of individuals of a species j of genus i in all the sites
  x_ij_bullet <- tapply(sxyz[,"abun"], list(sxyz[,"genus"],sxyz[,"sp"]), sum)
  x_ij_bullet[is.na(x_ij_bullet)] = 0 
  # cat("Total individuals by species:", x_ij_bullet, "\n")
  
  # the total number of individuals of all species of different genus found in the sampling site k
  x_bullet_bullet_k <- c(tapply(sxyz[,"abun"], sxyz[,"site"], sum))
  # cat("Total individuals by site:", x_bullet_bullet_k, "\n")
  
  
  # Mean total number of individuals across all sampling sites
  x_bar_bullet_bullet_bullet <- x_bullet_bullet_bullet / R
  
  # Mean total number of individuals of each species across all sampling sites
  x_bar_ij_bullet <- x_ij_bullet / R
  
  # Mean total number of individuals in each sampling site
  x_bar_bullet_bullet_k <- x_bullet_bullet_k / S
  
  # Partitioning of multi-species Poisson index of dispersion
  I_total <- (sum(x_bullet_bullet_k^2) - R * x_bar_bullet_bullet_bullet^2) / (x_bar_bullet_bullet_bullet)
  cat("Multi-species Poisson Index of dispersion:", I_total , "\n")
  
  # Calculate the within-genus within-species component
  I_within_genus_within_species <-  (sum(sxyz[,"abun"]^2) - R * sum(x_bar_ij_bullet^2))/x_bar_bullet_bullet_bullet
  cat("within-genus within-species component:", I_within_genus_within_species, "\n")
  
  # within-genus between-species synchrony component
  p1 = 0
  for (k in 1:R) {
    for (i in 1:G) {
      for (j in Sgn[[i]]) {
        for (l in Sgn[[i]]) {
          if (l != j) {
            x_ij_k = sxyz[sxyz[,"genus"]==i&sxyz[,"sp"]==j&sxyz[,"site"]==k,"abun"]
            if(length(x_ij_k)  == 0){
              x_ij_k = 0
            } else {
              x_ij_k = x_ij_k
            }
            x_il_k = sxyz[sxyz[,"genus"]==i&sxyz[,"sp"]==l&sxyz[,"site"]==k,"abun"]
            if(length(x_il_k)  == 0){
              x_il_k = 0
            } else {
              x_il_k = x_il_k
            }
            p1 = p1 + x_ij_k * x_il_k
            # cat("site:",k," genus:",i," species(j):",j," species(l):",l, " x_ij_k:",round(x_ij_k,1)," x_il_k:",round(x_il_k,1)," p1:", p1, "\n")
          } # if
        } # l species loop
      } # j species loop
    }  # genus loop
  } # site loop
  
  p2 = 0
  for (i in 1:G) {
    for (j in Sgn[[i]]) {
      for (l in Sgn[[i]]) {
        if (l != j) {
          p2 = p2 + x_bar_ij_bullet[i,j] * x_bar_ij_bullet[i,l]
          # cat("genus:",i," species(j):",j," species(l):",l, " x_bar_ij_bullet:",round(x_bar_ij_bullet[i,j],1)," x_bar_il_bullet:",round(x_bar_ij_bullet[i,l],1)," p2:", p2, "\n")
        } # if
      } # l species loop
    } # j species loop
  }  # genus loop
  
  I_within_genus_between_species = as.numeric((p1 - R * p2)/x_bar_bullet_bullet_bullet)
  cat("within-genus between-species synchrony component:", I_within_genus_between_species, "\n")
  
  # between-genus synchrony component
  p1 = 0
  for (k in 1:R) {
    for (i in 1:G) {
      for (u in 1:G) {
        if (u != i) {
          for (j in Sgn[[i]]) {
            for (v in Sgn[[u]]) {
              x_ij_k = sxyz[sxyz[,"genus"]==i & sxyz[,"sp"]==j & sxyz[,"site"]==k,"abun"]     
              if(length(x_ij_k)  == 0){
                x_ij_k = 0
              } else {
                x_ij_k = x_ij_k
              }
              x_uv_k = sxyz[sxyz[,"genus"]==u & sxyz[,"sp"]==v & sxyz[,"site"]==k,"abun"]
              if(length(x_uv_k)  == 0){
                x_uv_k = 0
              } else {
                x_uv_k = x_uv_k
              }
              p1 = p1 + x_ij_k * x_uv_k
              # cat("site:",k," genus(i):",i," genus(u):",u," species(j):",j," species(v):",v, " x_ij_k:",round(x_ij_k,1)," x_uv_k:",round(x_uv_k,1)," p1:", p1, "\n")
            } # v species loop
          } # j species loop
        } # if loop
      } # u genus loop
    }  # i genus loop
  } # k site loop
  
  p2 = 0
  for (i in 1:G) {
    for (u in 1:G) {
      if (u != i) {
        for (j in Sgn[[i]]) {
          for (v in Sgn[[u]]) {
            p2 = p2 + x_bar_ij_bullet[i,j] * x_bar_ij_bullet[u,v]
            # cat(" genus(i):",i," genus(u):",u," species(j):",j," species(v):",v, " x_bar_ij_bullet:",round(x_bar_ij_bullet[i,j],1)," x_bar_uv_bullet:",round(x_bar_ij_bullet[u,v],1)," p2:", p2, "\n")
          } 
        } 
      } 
    }  
  } 
  
  I_between_genus = as.numeric((p1 - R * p2)/x_bar_bullet_bullet_bullet)
  cat("between-genus synchrony component:", I_between_genus, "\n")
  
  
  # Check equality
  if (abs(I_within_genus_within_species +  I_within_genus_between_species + I_between_genus - I_total) < 1e-10) {
    message("Total equals within genus within species plus within genus between species plus between genus")
  } else {
    message("Total not equal to within_genus_within_species plus within_genus_between_species plus between_genus!!!")
  }
  
  #
  if(!adj.w %in% c("not","weight")){
    stop(paste0("The adjust parameter adj.w - ", adj.w, " is not supported yet"))
  }
  
  if(adj.w == "not"){
    Res <- data.frame(I.total=c(I_total,I_total),
                      I.within.genus.within.sp=c(I_within_genus_within_species,I_within_genus_within_species),
                      I.within.genus.between.sp =c(I_within_genus_between_species,I_within_genus_between_species),
                      I.between.genus=c(I_between_genus,I_between_genus))
    rownames(Res) <- c("original", adj.w)
    return(Res) 
  }  # if not adjust
  
  #
  # if(adj.w == "equal"){
  #   I_comp = c(I_within_genus_within_species,I_within_genus_between_species,I_between_genus)
  #   
  #   if(sum(sign(I_comp))== 3 | sum(sign(I_comp))== 2){
  #     Res <- data.frame(I_total=c(I_total,I_total),
  #                       I_within.genus.within.species=c(I_within_genus_within_species,I_within_genus_within_species),
  #                       I_within.genus.between.species =c(I_within_genus_between_species,I_within_genus_between_species),
  #                       I_between.genus=c(I_between_genus,I_between_genus))
  #     rownames(Res) <- c("original", adj.w)
  #     return(Res) 
  #     } # if no negative or two positive with zero
  #   
  #   if(sum(sign(I_comp)) == 1){
  #       neg_pos <- which(I_comp < 0)
  #       I_comp[-neg_pos] <- I_comp[-neg_pos] + (I_comp[neg_pos]/2)
  #       I_comp[neg_pos] <- 0
  #       if(sum(sign(I_comp)) == 0){
  #         neg_pos <- which(I_comp < 0)
  #         I_comp[-neg_pos] <- I_comp[-neg_pos] + (I_comp[neg_pos]/2)
  #         I_comp[neg_pos] <- 0
  #         Res <- data.frame(I_total=c(I_total,I_total),
  #                           I_within.genus.within.species=c(I_within_genus_within_species,I_comp[1]),
  #                           I_within.genus.between.species =c(I_within_genus_between_species,I_comp[2]),
  #                           I_between.genus=c(I_between_genus,I_comp[3]))
  #         rownames(Res) <- c("original", adj.w)
  #         return(Res) 
  #       }
  #     } # if one negative 
  #   
  #   if(sum(sign(I_comp))==-1){
  #     neg_pos <- which(I_comp < 0)
  #     I_comp[neg_pos] = 0 
  #     I_comp[-neg_pos] = I_total
  #     Res <- data.frame(I_total=c(I_total,I_total),
  #                   I_within.genus.within.species=c(I_within_genus_within_species,I_comp[1]),
  #                   I_within.genus.between.species =c(I_within_genus_between_species,I_comp[2]),
  #                   I_between.genus=c(I_between_genus,I_comp[3]))
  #     rownames(Res) <- c("original", adj.w)
  #     return(Res)
  #     } # if two negative 
  #  } # if equal
  
  #  
  if(adj.w == "weight"){
    I_comp = c(I_within_genus_within_species,I_within_genus_between_species,I_between_genus)
    
    if(sum(sign(I_comp))== 3 | sum(sign(I_comp))== 2){
      Res <- data.frame(I.total=c(I_total,I_total),
                        I.within.genus.within.sp=c(I_within_genus_within_species,I_within_genus_within_species),
                        I.within.genus.between.sp =c(I_within_genus_between_species,I_within_genus_between_species),
                        I.between.genus=c(I_between_genus,I_between_genus))
      rownames(Res) <- c("original", adj.w)
      return(Res) 
      } # if no negative or two positive with zero
    
    if(sum(sign(I_comp))==1){
      neg_pos <- which(I_comp < 0)
      opt_pos <- which(I_comp >= 0)
      I_comp[neg_pos] = 0
      I_comp <- I_comp * I_total / sum(I_comp)
        
      Res <- data.frame(I.total=c(I_total,I_total),
                      I.within.genus.within.sp=c(I_within_genus_within_species,I_comp[1]),
                      I.within.genus.between.sp =c(I_within_genus_between_species,I_comp[2]),
                      I.between.genus=c(I_between_genus,I_comp[3]))
      rownames(Res) <- c("original", adj.w)
      return(Res) 
    } # if one negative 
    
    if(sum(sign(I_comp))==-1){
      neg_pos <- which(I_comp < 0)
      I_comp[neg_pos] = 0 
      I_comp[-neg_pos] = I_total
     
      Res <- data.frame(I.total=c(I_total,I_total),
                      I.within.genus.within.sp = c(I_within_genus_within_species,I_comp[1]),
                      I.within.genus.between.sp =c(I_within_genus_between_species,I_comp[2]),
                      I.between.genus=c(I_between_genus,I_comp[3]))
      rownames(Res) <- c("original", adj.w)
      return(Res) 
      } # if two negative 
   } # if weight
    
} 

# EOF
################################################################################ ----


################################################################################ .
# partitioning abundance dispersion into within and between-site components      ----
# ZL - 2023-02-15 (Eq.2)
# 'sxyz' is a six-column matrix
# the first column represents genus ID
# the second column represents species ID
# the third and fourth columns represented x, and y coordinates
# the five column represent the numeric identifier for each site
# the final column represent abundance in each location 
# (only used for quadrat sampling and dispersion partitioning)
# The adj.w parameter offers the possibility to adjust the two partitioning components with ecological meaning
# If adj.w is assigned as "not", the output of the function will return the original value, regardless of whether the two partition components are negative.
# If adj.w is assigned as "equal", the bisection method will be used to deal with the negative value.
# If adj.w is assigned as "weight", while ignoring negative values, the remaining components are reallocated according to their weights.

# simple for loop version
I_site_partitioning_simple <- function(sxyz, adj.w = c("not","weight"))
{ 
 
  # total number of genus
  G = length(unique(sxyz[,"genus"]))
  cat("Total number of genus:", G, "\n")
  
  # total number of species
  S = length(unique(sxyz[,"sp"]))
  cat("Total number of species:", S, "\n")
  
  # total number of sampling sites
  R = length(unique(sxyz[,"site"]))
  cat("Total number of sampling sites:", R, "\n")
  
  # total number of individuals of all species and all sites
  x_bullet_bullet <- sum(sxyz[,"abun"])
  cat("Total individuals:", x_bullet_bullet, "\n")
  
  # the total number of individuals of a species i in all the sites
  x_i_bullet <- tapply(sxyz[,"abun"], sxyz[,"sp"], sum)
  # cat("Total individuals by species:", x_i_bullet, "\n")
  
  # the total number of individuals of different species found in the sampling site j
  x_bullet_j <- c(tapply(sxyz[,"abun"], sxyz[,"site"], sum))
  # cat("Total individuals by site:", x_bullet_j, "\n")
  
  # Mean total number of individuals across all sampling sites
  x_bar_bullet_bullet <- x_bullet_bullet / R
  
  # Mean total number of individuals of each species across all sampling sites
  x_bar_i_bullet <- x_i_bullet / R
  
  # Mean total number of individuals in each sampling site
  x_bar_bullet_j <- x_bullet_j / S
  
  # Partitioning of multi-species Poisson index of dispersion
  I_total <- (sum(x_bullet_j^2) - R * x_bar_bullet_bullet^2) / (x_bar_bullet_bullet)
  cat("Multi-species Poisson Index of dispersion:", I_total , "\n")
  
  # Calculate the within-site intraspecific abundance variation component
  I_within_site <- sum(sapply(1:R, 
                              function(j){(R/S)*sum((sxyz[sxyz[,"site"] == j, "abun"])^2)}) - S*(x_bar_bullet_j^2)) / sum(x_bar_bullet_j)
  cat("Within-site interspecific abundance variation component:", I_within_site, "\n")
  
  if(I_total - I_within_site > 0){
  # Calculate the between-site intraspecific abundance synchrony component
  # component 1
  p1 = 0
  for (j in 1:R) {
    for (i in 1:S) {
      for (k in 1:S) {
        if (k != i) {
          x_ij = sxyz[sxyz[,"sp"]==i&sxyz[,"site"]==j,"abun"]
          x_kj = sxyz[sxyz[,"sp"]==k&sxyz[,"site"]==j,"abun"]
          p1 = p1 + x_ij * x_kj
          # cat("site(j):",j," sp1(i):",i," sp2(k):",k, " x_ij:",x_ij,"x_kj:",x_kj," p1:", p1, "\n")
        } # if
      } # k
    } # i
  } # j
  p1 = p1 * (R/S)
  # component 2
  p2 = 0
  for (j in 1:R) {
    for (l in 1:R) {
      if (l != j) {
        x_bar_bullet_j_j = x_bar_bullet_j[j]
        x_bar_bullet_j_l = x_bar_bullet_j[l]
        p2 = p2 +  x_bar_bullet_j_j * x_bar_bullet_j_l
        # cat("site_1(j):",j," site_2(l):",l, " x_bar_bullet_j_j:",round(x_bar_bullet_j_j,1)," x_bar_bullet_j_l:",round(x_bar_bullet_j_l,1)," p2:", p2, "\n")
      }
    }
  }
  p2 = p2 * S
  # component 3
  p3 = sum(x_bar_bullet_j)
  
  
  I_between_site = (p1-p2)/p3
  cat("between-site interspecific abundance synchrony component:", I_between_site, "\n")
  } else{
  I_between_site = I_total - I_within_site
  cat("between-site interspecific abundance synchrony component:", I_between_site, "\n")
    
  }
  # Check equality
  if (abs(I_within_site + I_between_site - I_total) < 1e-10) {
    message("Total equals within component plus between component")
  } else {
    message("Total not equal to within component plus between component!!!")
  }
  
  # Adjustment
  if(!adj.w %in% c("not","weight")){
    stop(paste0("The adjust parameter adj.w - ", adj.w, " is not supported yet"))
  }
  
  if(adj.w == "not"){
    Res <- data.frame(I.total=c(I_total,I_total),I.within.site=c(I_within_site,I_within_site),
                      I.between.site=c(I_between_site,I_between_site))
    rownames(Res) <- c("original", adj.w)
    return(Res) 
  } 
  if(adj.w == "weight"){
    if(sum(sign(c(I_within_site,I_between_site)))==2){
      Res <- data.frame(I.total=c(I_total,I_total),I.within.site=c(I_within_site,I_within_site),
                        I.between.site=c(I_between_site,I_between_site))
      rownames(Res) <- c("original", adj.w)
      return(Res)  
    }
    if(I_within_site < 0){
      I_within_site_aj = 0; I_between_site_aj = I_total
    }
    if(I_between_site < 0){
      I_between_site_aj = 0; I_within_site_aj = I_total
    }
    Res <- data.frame(I.total=c(I_total,I_total),I.within.site=c(I_within_site,I_within_site_aj),
                      I.between.site=c(I_between_site,I_between_site_aj))
    rownames(Res) <- c("original", adj.w)
    return(Res)  
  }
}

# parallel computing version
I_site_partitioning_parall <- function(sxyz, adj.w = c("not","weight"),core.use = 0.9)
{ 
  
  # total number of genus
  G = length(unique(sxyz[,"genus"]))
  cat("Total number of genus:", G, "\n")
  
  # total number of species
  S = length(unique(sxyz[,"sp"]))
  cat("Total number of species:", S, "\n")
  
  # total number of sampling sites
  R = length(unique(sxyz[,"site"]))
  cat("Total number of sampling sites:", R, "\n")
  
  # total number of individuals of all species and all sites
  x_bullet_bullet <- sum(sxyz[,"abun"])
  cat("Total individuals:", x_bullet_bullet, "\n")
  
  # the total number of individuals of a species i in all the sites
  x_i_bullet <- c(tapply(sxyz[,"abun"], sxyz[,"sp"], sum))
  # cat("Total individuals by species:", x_i_bullet, "\n")
  
  # the total number of individuals of different species found in the sampling site j
  x_bullet_j <- c(tapply(sxyz[,"abun"], sxyz[,"site"], sum))
  # cat("Total individuals by site:", x_bullet_j, "\n")
  
  # Mean total number of individuals across all sampling sites
  x_bar_bullet_bullet <- x_bullet_bullet / R
  
  # Mean total number of individuals of each species across all sampling sites
  x_bar_i_bullet <- x_i_bullet / R
  
  # Mean total number of individuals in each sampling site
  x_bar_bullet_j <- x_bullet_j / S
  
  # Partitioning of multi-species Poisson index of dispersion
  I_total <- (sum(x_bullet_j^2) - R * x_bar_bullet_bullet^2) / (x_bar_bullet_bullet)
  cat("Multi-species Poisson Index of dispersion:", I_total , "\n")
  
  # Calculate the within-site intraspecific abundance variation component
  I_within_site <- sum(sapply(1:R, 
                              function(j){(R/S)*sum((sxyz[sxyz[,"site"] == j, "abun"])^2)}) - S*(x_bar_bullet_j^2)) / sum(x_bar_bullet_j)
  cat("Within-site interspecific abundance variation component:", I_within_site, "\n")
  
  if(I_total - I_within_site > 0){
  
  # Calculate the between-site intraspecific abundance synchrony component
  library(doParallel)
  
  # set up parallel computing with number of cores
  no_cores <- floor(detectCores() * core.use)
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  # Calculate the between-site intraspecific abundance synchrony component
  # component 1
  p1 = foreach(j = 1:R, .combine = "+") %dopar% {
    p1_j = 0
    for (i in 1:S) {
      for (k in 1:S) {
        if (k != i) {
          x_ij = sxyz[sxyz[,"sp"]==i&sxyz[,"site"]==j,"abun"]
          x_kj = sxyz[sxyz[,"sp"]==k&sxyz[,"site"]==j,"abun"]
          p1_j = p1_j + x_ij * x_kj
        } # if
      } # k
    } # i
    return(p1_j)
  }
  p1 = sum(p1) * (R/S)
  
  # component 2
  p2 = foreach(j = 1:R, .combine = "+") %dopar% {
    p2_j = 0
    for (l in 1:R) {
      if (l != j) {
        x_bar_bullet_j_j = x_bar_bullet_j[j]
        x_bar_bullet_j_l = x_bar_bullet_j[l]
        p2_j = p2_j +  x_bar_bullet_j_j * x_bar_bullet_j_l
      }
    }
    return(p2_j)
  }
  p2 = sum(p2) * S
  
  # component 3
  p3 = sum(x_bar_bullet_j)
  
  I_between_site = (p1 - p2) / p3
  cat("between-site interspecific abundance synchrony component:", I_between_site, "\n")
  stopCluster(cl)
  
  } else {
  
  I_between_site = I_total - I_within_site
  cat("between-site interspecific abundance synchrony component:", I_between_site, "\n")  
      
  }
  # Check equality
  if (abs(I_within_site + I_between_site - I_total) < 1e-10) {
    message("Total equals within component plus between component")
  } else {
    message("Total not equal to within component plus between component!!!")
  }
  
  # Adjustment
  if(!adj.w %in% c("not","weight")){
    stop(paste0("The adjust parameter adj.w - ", adj.w, " is not supported yet"))
  }
  
  if(adj.w == "not"){
    Res <- data.frame(I.total=c(I_total,I_total),I_within.site=c(I_within_site,I_within_site),
                      I.between.site=c(I_between_site,I_between_site))
    rownames(Res) <- c("original", adj.w)
    return(Res) 
  } 
  if(adj.w == "weight"){
    if(sum(sign(c(I_within_site,I_between_site)))==2){
      Res <- data.frame(I.total=c(I_total,I_total),I_within.site=c(I_within_site,I_within_site),
                        I.between.site=c(I_between_site,I_between_site))
      rownames(Res) <- c("original", adj.w)
      return(Res)  
    }
    if(I_within_site < 0){
      I_within_site_aj = 0; I_between_site_aj = I_total
    }
    if(I_between_site < 0){
      I_between_site_aj = 0; I_within_site_aj = I_total
    }
    Res <- data.frame(I.total=c(I_total,I_total),I_within.site=c(I_within_site,I_within_site_aj),
                      I.between.site=c(I_between_site,I_between_site_aj))
    rownames(Res) <- c("original", adj.w)
    return(Res)  
  }
}

# EOF
################################################################################ ----


################################################################################ .
# partitioning of abundance dispersion into within-site within-genus,            ----
# within-site between-genus and between-site components
# ZL - 2023-03-09 (Eq.4)
# 'sxyz' is a six-column matrix
# the first column represents genus ID
# the second column represents species ID
# the third and fourth columns represented x, and y coordinates
# the five column represent the numeric identifier for each site
# the final column represent abundance in each location 
# (only used for quadrat sampling and dispersion partitioning)
# The adj.w parameter offers the possibility to adjust the two partitioning components with ecological meaning
# If adj.w is assigned as "not", the output of the function will return the original value, regardless of whether the two partition components are negative.
# If adj.w is assigned as "equal", the bisection method will be used to deal with the negative value.
# If adj.w is assigned as "weight", while ignoring negative values, the remaining components are reallocated according to their weights.

I_site_genus_partitioning_1 <- function(sxyz,adj.w = c("not","weight"))
{ 
  # total number of genus
  G = length(unique(sxyz[,"genus"]))
  cat("Total number of genus:", G, "\n")
  
  # total number of species
  S = length(unique(sxyz[,"sp"]))
  cat("Total number of species:", S, "\n")
  
  # total number of sampling sites
  R = length(unique(sxyz[,"site"]))
  cat("Total number of sampling sites:", R, "\n")
  
  # Calculate total number of species per genus
  Sg <- tapply(sxyz[,"sp"], sxyz[,"genus"], length)/R
  Sgn <- tapply(sxyz[,"sp"], sxyz[,"genus"], FUN = unique)
  cat("Total number of species per genus:", Sg, "\n")
  
  # total number of individuals of all species and all sites
  x_bullet_bullet_bullet <- sum(sxyz[,"abun"])
  cat("Total individuals:", x_bullet_bullet_bullet, "\n")
  
  # the total number of individuals of a species j of genus i in all the sites
  x_ij_bullet <- tapply(sxyz[,"abun"], list(sxyz[,"genus"],sxyz[,"sp"]), sum)
  x_ij_bullet[is.na(x_ij_bullet)] = 0 
  cat("Total individuals by species:", x_ij_bullet, "\n")
  
  # the total number of individuals of different species of different genus found in the sampling site k
  x_bullet_bullet_k <- c(tapply(sxyz[,"abun"], sxyz[,"site"], sum))
  cat("Total individuals by site:", x_bullet_bullet_k, "\n")
  
  # the total number of individuals of all species of genus i in site k
  x_i_bullet_k <- tapply(sxyz[,"abun"], list(sxyz[,"genus"],sxyz[,"site"]), sum)
  x_i_bullet_k[is.na(x_i_bullet_k)] = 0 
  cat("Total individuals by genus and site:", x_i_bullet_k, "\n")
  
  # Mean total number of individuals of all species of genus i in site k
  x_bar_i_bullet_k = as.matrix(x_i_bullet_k) / matrix(Sg,nrow=length(Sg),ncol=R)
  
  # Mean total number of individuals across all sampling sites
  x_bar_bullet_bullet_bullet <- x_bullet_bullet_bullet / R
  
  # Mean total number of individuals of each species across all sampling sites
  x_bar_ij_bullet <- x_ij_bullet / R
  
  # Mean total number of individuals in each sampling site
  x_bar_bullet_bullet_k <- x_bullet_bullet_k / S
  
  # Partitioning of multi-species Poisson index of dispersion
  I_total <- (sum(x_bullet_bullet_k^2) - R * x_bar_bullet_bullet_bullet^2) / (x_bar_bullet_bullet_bullet)
  cat("Multi-species Poisson Index of dispersion:", I_total , "\n")
  
  # Calculate the within-site within-genus component
  p1 = 0
  for (k in 1:R) {
    for (i in 1:G) {
      for(j in Sgn[[i]]){
        x_ij_k = sxyz[sxyz[,"genus"]==i & sxyz[,"sp"]==j & sxyz[,"site"]==k,"abun"]
        if(length(x_ij_k)  == 0){
          x_ij_k = 0
        } else {
          x_ij_k = x_ij_k
        }
        p1 = p1 + x_ij_k^2
      } # j loop
    } # i loop
  } # k loop
  p1 = p1 * R / S 
  
  p2 = 0
  for (k in 1:R) {
    for (i in 1:G) {
      p2 = p2 + (Sg[i]/S)^2 * x_bar_i_bullet_k[i,k]^2
    }
  }
  p2 = S*p2
  
  p3 = 0
  for (k in 1:R) {
    p3 = p3 + x_bar_bullet_bullet_k[k]
  }
  
  I_within_site_within_genus <-  as.numeric((p1 - p2) / p3)
  cat("within-site within-genus component:", I_within_site_within_genus, "\n")
  
  # within-site between-genus synchrony component
  p1 = 0
  for (k in 1:R) {
    for (i in 1:G) {
      for (j in Sgn[[i]]) {
        for (l in Sgn[[i]]) {
          if (l != j) {
            x_ij_k = sxyz[sxyz[,"genus"]==i&sxyz[,"sp"]==j&sxyz[,"site"]==k,"abun"]
            if(length(x_ij_k)  == 0){
              x_ij_k = 0
            } else {
              x_ij_k = x_ij_k
            }
            x_il_k = sxyz[sxyz[,"genus"]==i&sxyz[,"sp"]==l&sxyz[,"site"]==k,"abun"]
            if(length(x_il_k)  == 0){
              x_il_k = 0
            } else {
              x_il_k = x_il_k
            }
            p1 = p1 + x_ij_k * x_il_k
            # cat("site:",k," genus:",i," species(j):",j," species(l):",l, " x_ij_k:",round(x_ij_k,1)," x_il_k:",round(x_il_k,1)," p1:", p1, "\n")
          } # if
        } # l species loop
      } # j species loop
    }  # genus loop
  } # site loop
  p1 = p1 * R / S 
  
  
  p2 = 0
  for (k in 1:R) {
    for (i in 1:G) {
      for (j in 1:G) {
        if (j != i) {
          p2 = p2 + ((Sg[i]*Sg[j])/(S*S))*(x_bar_i_bullet_k[i,k] * x_bar_i_bullet_k[j,k])
        } # if
      } # l species loop
    }  # genus loop
  }
  p2 = S*p2
  
  p3 = 0
  for (k in 1:R) {
    p3 = p3 + x_bar_bullet_bullet_k[k]
  }
  
  I_within_site_between_genus = as.numeric((p1 - p2) / p3)
  cat("within-site between-genus synchrony component:", I_within_site_between_genus, "\n")
  
  
  # between-site synchrony component
  p1 = 0
  for (k in 1:R) {
    for (i in 1:G) {
      for (u in 1:G) {
        if (u != i) {
          for (j in Sgn[[i]]) {
            for (v in Sgn[[u]]) {
              x_ij_k = sxyz[sxyz[,"genus"]==i & sxyz[,"sp"]==j & sxyz[,"site"]==k,"abun"]     
              if(length(x_ij_k)  == 0){
                x_ij_k = 0
              } else {
                x_ij_k = x_ij_k
              }
              x_uv_k = sxyz[sxyz[,"genus"]==u & sxyz[,"sp"]==v & sxyz[,"site"]==k,"abun"]
              if(length(x_uv_k)  == 0){
                x_uv_k = 0
              } else {
                x_uv_k = x_uv_k
              }
              p1 = p1 + x_ij_k * x_uv_k
              # cat("site:",k," genus(i):",i," genus(u):",u," species(j):",j," species(v):",v, " x_ij_k:",round(x_ij_k,1)," x_uv_k:",round(x_uv_k,1)," p1:", p1, "\n")
            } # v species loop
          } # j species loop
        } # if loop
      } # u genus loop
    }  # i genus loop
  } # k site loop
  
  p1 = p1 * R / S 
  
  p2 = 0
  for (k in 1:R) {
    for (l in 1:R) {
      if (l != k) {
        for (i in 1:G) {
          for (j in 1:G) {
            p2 = p2 + ((Sg[i]*Sg[j])/(S*S))*(x_bar_i_bullet_k[i,k] * x_bar_i_bullet_k[j,l])
          }
        }
      }
    }
  }
  p2 = S * p2
  
  p3 = 0
  for (k in 1:R) {
    p3 = p3 + x_bar_bullet_bullet_k[k]
  }
  
  
  I_between_site = as.numeric((p1 - p2) / p3)
  
  cat("between-site synchrony component:", I_between_site, "\n")
  
  
  # Check equality
  if (abs(I_within_site_within_genus +  I_within_site_between_genus + I_between_site - I_total) < 1e-10) {
    message("Total equals within site within genus plus within site between genus plus between site")
  } else {
    message("Total not equal to within site within genus plus within site between genus plus between site!!!")
  }
  
  
  #
  if(!adj.w %in% c("not","weight")){
    stop(paste0("The adjust parameter adj.w - ", adj.w, " is not supported yet"))
  }
  
  if(adj.w == "not"){
    Res <- data.frame(I.total=c(I_total,I_total),
                      I.within.site.within.genus=c(I_within_site_within_genus,I_within_site_within_genus),
                      I.within.site.between.genus =c(I_within_site_between_genus,I_within_site_between_genus),
                      I.between.site=c(I_between_site,I_between_site))
    rownames(Res) <- c("original", adj.w)
    return(Res) 
  }  # if not adjust
  
  #  
  if(adj.w == "weight"){
    I_comp = c(I_within_site_within_genus,I_within_site_between_genus,I_between_site)
    
    if(sum(sign(I_comp))== 3 | sum(sign(I_comp))== 2){
      Res <- data.frame(I.total=c(I_total,I_total),
                        I.within.site.within.genus=c(I_within_site_within_genus,I_within_site_within_genus),
                        I.within.site.between.genus =c(I_within_site_between_genus,I_within_site_between_genus),
                        I.between.site=c(I_between_site,I_between_site))
      rownames(Res) <- c("original", adj.w)
      return(Res) 
    } # if no negative or two positive with zero
    
    if(sum(sign(I_comp))==1){
      neg_pos <- which(I_comp < 0)
      opt_pos <- which(I_comp >= 0)
      I_comp[neg_pos] = 0
      I_comp <- I_comp * I_total / sum(I_comp)
      
      Res <- data.frame(I.total=c(I_total,I_total),
                        I.within.site.within.genus=c(I_within_site_within_genus,I_comp[1]),
                        I.within.site.between.genus =c(I_within_site_between_genus,I_comp[2]),
                        I.between.site=c(I_between_site,I_comp[3]))
      rownames(Res) <- c("original", adj.w)
      return(Res) 
    } # if one negative 
    
    if(sum(sign(I_comp))==-1){
      neg_pos <- which(I_comp < 0)
      I_comp[neg_pos] = 0 
      I_comp[-neg_pos] = I_total
      
      Res <- data.frame(I.total=c(I_total,I_total),
                        I.within.site.within.genus=c(I_within_site_within_genus,I_comp[1]),
                        I.within.site.between.genus =c(I_within_site_between_genus,I_comp[2]),
                        I.between.site=c(I_between_site,I_comp[3]))
      rownames(Res) <- c("original", adj.w)
      return(Res) 
    } # if two negative 
  } # if weight
  
} 

# EOF
################################################################################ ----


################################################################################ .
# partitioning of abundance dispersion into within-site,                          ----
# between-site within-genus and between-site between-genus components
# ZL - 2023-03-12 (Eq.5)
# 'sxyz' is a six-column matrix
# the first column represents genus ID
# the second column represents species ID
# the third and fourth columns represented x, and y coordinates
# the five column represent the numeric identifier for each site
# the final column represent abundance in each location 
# (only used for quadrat sampling and dispersion partitioning)
# The adj.w parameter offers the possibility to adjust the two partitioning components with ecological meaning
# If adj.w is assigned as "not", the output of the function will return the original value, regardless of whether the two partition components are negative.
# If adj.w is assigned as "equal", the bisection method will be used to deal with the negative value.
# If adj.w is assigned as "weight", while ignoring negative values, the remaining components are reallocated according to their weights.

I_site_genus_partitioning_2 <- function(sxyz, adj.w = c("not","weight"))
{ 
  # total number of genus
  G = length(unique(sxyz[,"genus"]))
  cat("Total number of genus:", G, "\n")
  
  # total number of species
  S = length(unique(sxyz[,"sp"]))
  cat("Total number of species:", S, "\n")
  
  # total number of sampling sites
  R = length(unique(sxyz[,"site"]))
  cat("Total number of sampling sites:", R, "\n")
  
  # Calculate total number of species per genus
  Sg <- tapply(sxyz[,"sp"], sxyz[,"genus"], length)/R
  Sgn <- tapply(sxyz[,"sp"], sxyz[,"genus"], FUN = unique)
  cat("Total number of species per genus:", Sg, "\n")
  
  # total number of individuals of all species and all sites
  x_bullet_bullet_bullet <- sum(sxyz[,"abun"])
  cat("Total individuals:", x_bullet_bullet_bullet, "\n")
  
  # the total number of individuals of a species j of genus i in all the sites
  x_ij_bullet <- tapply(sxyz[,"abun"], list(sxyz[,"genus"],sxyz[,"sp"]), sum)
  x_ij_bullet[is.na(x_ij_bullet)] = 0 
  cat("Total individuals by species:", x_ij_bullet, "\n")
  
  # the total number of individuals of different species of different genus found in the sampling site k
  x_bullet_bullet_k <- c(tapply(sxyz[,"abun"], sxyz[,"site"], sum))
  cat("Total individuals by site:", x_bullet_bullet_k, "\n")
  
  # the total number of individuals of all species of genus i in site k
  x_i_bullet_k <- tapply(sxyz[,"abun"], list(sxyz[,"genus"],sxyz[,"site"]), sum)
  x_i_bullet_k[is.na(x_i_bullet_k)] = 0 
  cat("Total individuals by species:", x_i_bullet_k, "\n")
  
  # Mean total number of individuals of all species of genus i in site k
  x_bar_i_bullet_k = as.matrix(x_i_bullet_k) / matrix(Sg,nrow=length(Sg),ncol=R)
  
  # Mean total number of individuals across all sampling sites
  x_bar_bullet_bullet_bullet <- x_bullet_bullet_bullet / R
  
  # Mean total number of individuals of each species across all sampling sites
  x_bar_ij_bullet <- x_ij_bullet / R
  
  # Mean total number of individuals in each sampling site
  x_bar_bullet_bullet_k <- x_bullet_bullet_k / S
  
  # Partitioning of multi-species Poisson index of dispersion
  I_total <- (sum(x_bullet_bullet_k^2) - R * x_bar_bullet_bullet_bullet^2) / (x_bar_bullet_bullet_bullet)
  cat("Multi-species Poisson Index of dispersion:", I_total , "\n")
  
  # Calculate the within-site within-genus component
  p1 = 0
  for (k in 1:R) {
    for (i in 1:G) {
      for(j in Sgn[[i]]){
        x_ij_k = sxyz[sxyz[,"genus"]==i & sxyz[,"sp"]==j & sxyz[,"site"]==k,"abun"]
        if(length(x_ij_k)  == 0){
          x_ij_k = 0
        } else {
          x_ij_k = x_ij_k
        }
        p1 = p1 + x_ij_k^2
      } # j loop
    } # i loop
  } # k loop
  p1 = p1 * R / S 
  
  p2 = 0
  for (k in 1:R) {
    p2 = p2 + x_bar_bullet_bullet_k[k]^2
  }
  p2 = S*p2
  
  p3 = 0
  for (k in 1:R) {
    p3 = p3 + x_bar_bullet_bullet_k[k]
  }
  
  I_within_site <-  as.numeric((p1 - p2) / p3)
  cat("within-site component:", I_within_site, "\n")
  
  # within-site between-genus synchrony component
  p1 = 0
  for (k in 1:R) {
    for (i in 1:G) {
      for (j in Sgn[[i]]) {
        for (l in Sgn[[i]]) {
          if (l != j) {
            x_ij_k = sxyz[sxyz[,"genus"]==i&sxyz[,"sp"]==j&sxyz[,"site"]==k,"abun"]
            if(length(x_ij_k)  == 0){
              x_ij_k = 0
            } else {
              x_ij_k = x_ij_k
            }
            x_il_k = sxyz[sxyz[,"genus"]==i&sxyz[,"sp"]==l&sxyz[,"site"]==k,"abun"]
            if(length(x_il_k)  == 0){
              x_il_k = 0
            } else {
              x_il_k = x_il_k
            }
            p1 = p1 + x_ij_k * x_il_k
            # cat("site:",k," genus:",i," species(j):",j," species(l):",l, " x_ij_k:",round(x_ij_k,1)," x_il_k:",round(x_il_k,1)," p1:", p1, "\n")
          } # if
        } # l species loop
      } # j species loop
    }  # genus loop
  } # site loop
  p1 = p1 * R / S 
  
  
  p2 = 0
  for (k in 1:R) {
    for (l in 1:R) {
      if(l != k){
        for (i in 1:G) {
          p2 = p2 + ((Sg[i]/S)^2)*(x_bar_i_bullet_k[i,k] * x_bar_i_bullet_k[i,l])
        } 
      }  
    }
  } 
  p2 = S*p2
  
  p3 = 0
  for (k in 1:R) {
    p3 = p3 + x_bar_bullet_bullet_k[k]
  }
  
  I_between_site_within_genus = as.numeric((p1 - p2) / p3)
  cat("between-site within-genus synchrony component:",I_between_site_within_genus, "\n")
  
  
  # between-site synchrony component
  p1 = 0
  for (k in 1:R) {
    for (i in 1:G) {
      for (u in 1:G) {
        if (u != i) {
          for (j in Sgn[[i]]) {
            for (v in Sgn[[u]]) {
              x_ij_k = sxyz[sxyz[,"genus"]==i & sxyz[,"sp"]==j & sxyz[,"site"]==k,"abun"]     
              if(length(x_ij_k)  == 0){
                x_ij_k = 0
              } else {
                x_ij_k = x_ij_k
              }
              x_uv_k = sxyz[sxyz[,"genus"]==u & sxyz[,"sp"]==v & sxyz[,"site"]==k,"abun"]
              if(length(x_uv_k)  == 0){
                x_uv_k = 0
              } else {
                x_uv_k = x_uv_k
              }
              p1 = p1 + x_ij_k * x_uv_k
              # cat("site:",k," genus(i):",i," genus(u):",u," species(j):",j," species(v):",v, " x_ij_k:",round(x_ij_k,1)," x_uv_k:",round(x_uv_k,1)," p1:", p1, "\n")
            } # v species loop
          } # j species loop
        } # if loop
      } # u genus loop
    }  # i genus loop
  } # k site loop
  
  p1 = p1 * R / S 
  
  p2 = 0
  for (k in 1:R) {
    for (l in 1:R) {
      if (l != k) {
        for (i in 1:G) {
          for (j in 1:G) {
            if(i != j){
              p2 = p2 + ((Sg[i]*Sg[j])/(S*S))*(x_bar_i_bullet_k[i,k] * x_bar_i_bullet_k[j,l])
            }
          }
        }
      }
    }
  }
  p2 = S * p2
  
  p3 = 0
  for (k in 1:R) {
    p3 = p3 + x_bar_bullet_bullet_k[k]
  }
  
  
  I_between_site_between_genus = as.numeric((p1 - p2) / p3)
  cat("between-site between-genus synchrony component:", I_between_site_between_genus, "\n")
  
  
  # Check equality
  if (abs(I_within_site +  I_between_site_within_genus + I_between_site_between_genus - I_total) < 1e-10) {
    message("Total equals to the three partioning components")
  } else {
    message("Total not equal to the three partioning components!!!")
  }
  
  
  #
  if(!adj.w %in% c("not","weight")){
    stop(paste0("The adjust parameter adj.w - ", adj.w, " is not supported yet"))
  }
  
  if(adj.w == "not"){
    Res <- data.frame(I.total=c(I_total,I_total),
                      I.within.site = c(I_within_site,I_within_site),
                      I.between.site.within.genus = c(I_between_site_within_genus,I_between_site_within_genus),
                      I.between.site.between.genus =c(I_between_site_between_genus,I_between_site_between_genus))
    rownames(Res) <- c("original", adj.w)
    return(Res) 
  }  # if not adjust
  
  #  
  if(adj.w == "weight"){
    I_comp = c(I_within_site,I_between_site_within_genus,I_between_site_between_genus)
    
    if(sum(sign(I_comp))== 3 | sum(sign(I_comp))== 2){
      Res <- data.frame(I.total=c(I_total,I_total),
                        I.within.site = c(I_within_site,I_within_site),
                        I.between.site.within.genus = c(I_between_site_within_genus,I_between_site_within_genus),
                        I.between.site.between.genus =c(I_between_site_between_genus,I_between_site_between_genus))
      rownames(Res) <- c("original", adj.w)
      return(Res) 
    } # if no negative or two positive with zero
    
    if(sum(sign(I_comp))==1){
      neg_pos <- which(I_comp < 0)
      opt_pos <- which(I_comp >= 0)
      I_comp[neg_pos] = 0
      I_comp <- I_comp * I_total / sum(I_comp)
      
      Res <- data.frame(I.total=c(I_total,I_total),
                        I.within.site = c(I_within_site,I_comp[1]),
                        I.between.site.within.genus = c(I_between_site_within_genus,I_comp[2]),
                        I.between.site.between.genus =c(I_between_site_between_genus,I_comp[3]))
      rownames(Res) <- c("original", adj.w)
      return(Res) 
    } # if one negative 
    
    if(sum(sign(I_comp))==-1){
      neg_pos <- which(I_comp < 0)
      I_comp[neg_pos] = 0 
      I_comp[-neg_pos] = I_total
      
      Res <- data.frame(I.total=c(I_total,I_total),
                        I.within.site = c(I_within_site,I_comp[1]),
                        I.between.site.within.genus = c(I_between_site_within_genus,I_comp[2]),
                        I.between.site.between.genus =c(I_between_site_between_genus,I_comp[3]))
      
      rownames(Res) <- c("original", adj.w)
      return(Res) 
    } # if two negative 
  } # if weight
  
} 

# EOF
################################################################################ ----

################################################################################ .
# partitioning pairwise distances into within and between-species components     ----
# ZL - 2023-03-16 (Eq.6)
# 'sxyz' is a six-column matrix
# the first column represents genus ID
# the second column represents species ID
# the third and fourth columns represented x, and y coordinates
# the five column represent the numeric identifier for each site
# the final column represent abundance in each location 
# (only used for quadrat sampling and dispersion partitioning)

D_sp_partitioning <- function(sxyz,rescaled=FALSE,expit=TRUE)
{
  # select only quadrat with individuals
  sxyz = sxyz[which(sxyz[,"abun"]!=0),]
  
  # the species id of input sxyz matrix
  id = as.character(sxyz[,"sp"])
  
  # pairwise distance for sxyz matrix
  dxy =  parallelDist::parDist(sxyz[,c("x","y")], method = "euclidean")
  
  # how many combinations of species i in site j, n = S*R
  n = dim(sxyz)[1] # n is x_bullet_bullet
  
  # Choose any combination of 2 rows from S*R rows
  cb=t(combn(1:n,2)) 
  
  # Which ids are the rows that indicate species equivalence
  ids=which(id[cb[,1]]==id[cb[,2]]) 
  
  # How many rows in total of species equivalence are there?
  ni=length(ids) 
  
  # How many rows of unequal species are there?
  no=dim(cb)[1]-ni 
  
  # If ni = 0, i.e., the case where each species is distributed at a different site
  # total pairwise distance equals between-species species component
  
  # If no = 0, i.e., the case where all species are distributed at the same sites
  # in such a case, the total pairwise distance equals to within-species component
  
  if(ni==0) # 
  {
    din=0
    dout=sum(dxy)
  }else if(no==0)
  {
    din=sum(dxy)
    dout=0
  }else
  {
    din=sum(dxy[ids])
    dout=sum(dxy[-ids])
  }
  # rescaled
  if(rescaled==FALSE)
  {
    dtot=(din+dout)/(n*(n-1)/2)
    din=din/(n*(n-1)/2)
    dout=dout/(n*(n-1)/2)
  }else
  {
    dtot=(din+dout)/(n*(n-1)/2)
    din=din/(ni*(ni-1)/2+1e-300)
    dout=dout/(no*(no-1)/2+1e-300)
  }
  # exponential transformation
  if(expit==TRUE)
  {
    dtot=exp(-dtot)
    din=exp(-din)
    dout=exp(-dout)
  }
  # results
  return(list(D.total=dtot,D.within.sp=din,D.between.sp=dout))
}


# EOF
################################################################################ ----


################################################################################ .
# partitioning pairwise distances into within-site and between-site components   ----
# ZL - 2023-03-23 (Eq.6)
# 'sxyz' is a six-column matrix
# the first column represents genus ID
# the second column represents species ID
# the third and fourth columns represented x, and y coordinates
# the five column represent the numeric identifier for each site
# the final column represent abundance in each location 
# (only used for quadrat sampling and dispersion partitioning)

D_site_partitioning <- function(sxyz,rescaled=FALSE,expit=TRUE)
{
  
  # select only quadrat with individuals
  sxyz = sxyz[which(sxyz[,"abun"]!=0),]
  
  # the site id of input sxyz matrix
  id = as.character(sxyz[,"site"])
  
  # pairwise distance for sxyz matrix
  dxy = parallelDist::parDist(sxyz[,c("x","y")], method = "euclidean")
  
  # how many combinations of species i in site j, n = S*R
  n = dim(sxyz)[1] # n is x_bullet_bullet
  
  # Choose any combination of 2 rows from S*R rows
  cb=t(combn(1:n,2)) 
  
  # Which ids are the rows that indicate site equivalence
  ids=which(id[cb[,1]]==id[cb[,2]]) 
  
  # How many rows in total of site equivalence are there?
  ni=length(ids) 
  
  # How many rows of unequal site are there?
  no=dim(cb)[1]-ni 
  
  # If ni = 0, i.e., the case where each species is distributed at a different site
  # total pairwise distance equals between-site component
  
  # If no = 0, i.e., the case where all species are distributed at the same sites
  # in such a case, the total pairwise distance equals to within-site component
  
  if(ni==0) # 
  {
    din=0
    dout=sum(dxy)
  }else if(no==0)
  {
    din=sum(dxy)
    dout=0
  }else
  {
    din=sum(dxy[ids])
    dout=sum(dxy[-ids])
  }
  # rescaled
  if(rescaled==FALSE)
  {
    dtot=(din+dout)/(n*(n-1)/2)
    din=din/(n*(n-1)/2)
    dout=dout/(n*(n-1)/2)
  }else
  {
    dtot=(din+dout)/(n*(n-1)/2)
    din=din/(ni*(ni-1)/2+1e-300)
    dout=dout/(no*(no-1)/2+1e-300)
  }
  # exponential transformation
  if(expit==TRUE)
  {
    dtot=exp(-dtot)
    din=exp(-din)
    dout=exp(-dout)
  }
  # results
  return(list(D.total=dtot,D.within.site=din,D.between.site=dout))
}

# EOF
################################################################################ ----

################################################################################ .
# partitioning of pairwise distances into within-genus within-species distance,  ----
# within-genus between-species distance and between-genus distance components
# ZL - 2023-03-23 (Eq.7)
# 'sxyz' is a six-column matrix
# the first column represents genus ID
# the second column represents species ID
# the third and fourth columns represented x, and y coordinates
# the five column represent the numeric identifier for each site
# the final column represent abundance in each location 
# (only used for quadrat sampling and dispersion partitioning)

D_genus_sp_partitioning <- function(sxyz,rescaled=FALSE,expit=TRUE)
{
 
  # select only quadrat with individuals
  sxyz = sxyz[which(sxyz[,"abun"]!=0),]
  
  # the site id of input sxyz matrix
  id_site = as.character(sxyz[,"site"])
  
  # the species id of input sxyz matrix
  id_sp = as.character(sxyz[,"sp"])
  
  # the genus id of input sxyz matrix
  id_genus = as.character(sxyz[,"genus"])
  
  # pairwise distance for sxyz matrix
  dxy = parallelDist::parDist(sxyz[,c("x","y")], method = "euclidean")
  
  # how many combinations of species i in site j, n = S*R
  n = dim(sxyz)[1] # n is x_bullet_bullet
  
  # Choose any combination of 2 rows from S*R rows
  cb=t(combn(1:n,2)) 
  
  # within-genus within-species distance component
  id1=which(id_genus[cb[,1]]==id_genus[cb[,2]] & id_sp[cb[,1]]==id_sp[cb[,2]]) 
  
  # within-genus between-species distance component
  id2 = which(id_genus[cb[,1]]==id_genus[cb[,2]] & id_sp[cb[,1]]!=id_sp[cb[,2]]) 
  
  # between genus distance component
  id3 = which(id_genus[cb[,1]]!=id_genus[cb[,2]]) 
 
  dtot=sum(dxy)
  din.g.in.s=sum(dxy[id1])
  din.g.out.s=sum(dxy[id2])
  dout.g=sum(dxy[id3])
  
  # rescaled
  if(rescaled==FALSE)
  {
    dtot=(dtot)/(n*(n-1)/2)
    din.g.in.s=din.g.in.s/(n*(n-1)/2)
    din.g.out.s=din.g.out.s/(n*(n-1)/2)
    dout.g =dout.g/(n*(n-1)/2)
  }else
  {
    dtot=(dtot)/(n*(n-1)/2)
    din.g.in.s=din.g.in.s/(length(id1)*(length(id1)-1)/2+1e-300)
    din.g.out.s=din.g.out.s/(length(id2)*(length(id2)-1)/2+1e-300)
    dout.g =dout.g/(length(id3)*(length(id3)-1)/2)
  }
  # exponential transformation
  if(expit==TRUE)
  {
    dtot=exp(-dtot)
    din.g.in.s=exp(-din.g.in.s)
    din.g.out.s=exp(-din.g.out.s)
    dout.g =exp(-dout.g)
  }
  # results
  return(list(D.total=dtot,D.within.genus.within.species=din.g.in.s,
              D.within.genus.between.species=din.g.out.s,
              D.between.genus=dout.g))
}

# EOF
################################################################################ ----

################################################################################ .
# partitioning of pairwise distances into within-site within-genus distance,     ----
# within-site between-genus distance and between-site distance components   
# ZL - 2023-03-27 (Eq.7)
# 'sxyz' is a six-column matrix
# the first column represents genus ID
# the second column represents species ID
# the third and fourth columns represented x, and y coordinates
# the five column represent the numeric identifier for each site
# the final column represent abundance in each location 
# (only used for quadrat sampling and dispersion partitioning)

D_site_genus_partitioning_1 <- function(sxyz,rescaled=FALSE,expit=TRUE)
{
 
  # select only quadrat with individuals
  sxyz = sxyz[which(sxyz[,"abun"]!=0),]
  
  # the site id of input sxyz matrix
  id_site = as.character(sxyz[,"site"])
  
  # the genus id of input sxyz matrix
  id_genus = as.character(sxyz[,"genus"])
  
  # pairwise distance for sxyz matrix
  dxy = parallelDist::parDist(sxyz[,c("x","y")], method = "euclidean")
  
  # how many combinations of species i in site j, n = S*R
  n = dim(sxyz)[1] # n is x_bullet_bullet
  
  # Choose any combination of 2 rows from S*R rows
  cb=t(combn(1:n,2)) 
  
  # within-site within-genus distance component
  id1=which(id_site[cb[,1]] == id_site[cb[,2]] & id_genus[cb[,1]]==id_genus[cb[,2]] ) 
  
  # within-site between-genus distance component
  id2 = which(id_site[cb[,1]] == id_site[cb[,2]] & id_genus[cb[,1]]!=id_genus[cb[,2]]) 
  
  # between site distance component
  id3 = which(id_site[cb[,1]]!=id_site[cb[,2]]) 
  
  dtot=sum(dxy)
  din.si.in.g=sum(dxy[id1])
  din.si.out.g=sum(dxy[id2])
  dout.si=sum(dxy[id3])
  
  # rescaled
  if(rescaled==FALSE)
  {
    dtot=(dtot)/(n*(n-1)/2)
    din.si.in.g= din.si.in.g/(n*(n-1)/2)
    din.si.out.g=din.si.out.g/(n*(n-1)/2)
    dout.si =dout.si/(n*(n-1)/2)
  }else
  {
    dtot=(dtot)/(n*(n-1)/2)
    din.si.in.g=din.si.in.g/(length(id1)*(length(id1)-1)/2+1e-300)
    din.si.out.g=din.si.out.g/(length(id2)*(length(id2)-1)/2+1e-300)
    dout.si =dout.si/(length(id3)*(length(id3)-1)/2)
  }
  # exponential transformation
  if(expit==TRUE)
  {
    dtot=exp(-dtot)
    din.si.in.g=exp(-din.si.in.g)
    din.si.out.g=exp(-din.si.out.g)
    dout.si =exp(-dout.si)
  }
  # results
  return(list(D.total=dtot,D.within.site.within.genus=din.si.in.g,
              D.within.site.between.genus=din.si.out.g,
              D.between.site=dout.si))
}

# EOF
################################################################################ ----

################################################################################ .
# partitioning of pairwise distances into within-site within-genus distance,     ----
# within-site between-genus distance and between-site distance components   
# ZL - 2023-03-27 (Eq.7)
# 'sxyz' is a six-column matrix
# the first column represents genus ID
# the second column represents species ID
# the third and fourth columns represented x, and y coordinates
# the five column represent the numeric identifier for each site
# the final column represent abundance in each location 
# (only used for quadrat sampling and dispersion partitioning)

D_site_genus_partitioning_2 <- function(sxyz,rescaled=FALSE,expit=TRUE)
{
 
  # select only quadrat with individuals
  sxyz = sxyz[which(sxyz[,"abun"]!=0),]
  
  # the site id of input sxyz matrix
  id_site = as.character(sxyz[,"site"])
  
  # the genus id of input sxyz matrix
  id_genus = as.character(sxyz[,"genus"])
  
  # pairwise distance for sxyz matrix
  dxy = parallelDist::parDist(sxyz[,c("x","y")], method = "euclidean")
  
  # how many combinations of species i in site j, n = S*R
  n = dim(sxyz)[1] # n is x_bullet_bullet
  
  # Choose any combination of 2 rows from S*R rows
  cb=t(combn(1:n,2)) 
  
  # within-site distance component
  id1=which(id_site[cb[,1]] == id_site[cb[,2]] ) 
  
  # between-site within-genus distance component
  id2 = which(id_site[cb[,1]]!=id_site[cb[,2]] & id_genus[cb[,1]]==id_genus[cb[,2]]) 
  
  # between site between genus distance component
  id3 = which(id_site[cb[,1]]!=id_site[cb[,2]] & id_genus[cb[,1]]!=id_genus[cb[,2]]) 
  
  dtot=sum(dxy)
  din.si=sum(dxy[id1])
  dout.si.in.g=sum(dxy[id2])
  dout.si.out.g=sum(dxy[id3])
  
  # rescaled
  if(rescaled==FALSE)
  {
    dtot=(dtot)/(n*(n-1)/2)
    din.si= din.si/(n*(n-1)/2)
    dout.si.in.g=dout.si.in.g/(n*(n-1)/2)
    dout.si.out.g =dout.si.out.g/(n*(n-1)/2)
  }else
  {
    dtot=(dtot)/(n*(n-1)/2)
    din.si=din.si/(length(id1)*(length(id1)-1)/2+1e-300)
    dout.si.in.g=dout.si.in.g/(length(id2)*(length(id2)-1)/2+1e-300)
    dout.si.out.g =dout.si.out.g/(length(id3)*(length(id3)-1)/2)
  }
  # exponential transformation
  if(expit==TRUE)
  {
    dtot=exp(-dtot)
    din.si=exp(-din.si)
    dout.si.in.g=exp(-dout.si.in.g)
    dout.si.out.g =exp(-dout.si.out.g)
  }
  # results
  return(list(D.total=dtot,
              D.within.site=din.si,
              D.between.site.within.genus=dout.si.in.g,
              D.between.site.between.genus=dout.si.out.g))
}

# EOF

################################################################################ ----



# numeric test                                                                   ----
sxyz = data_creation(S=5, G=3, R=12, lambda=3,fix.rd = T)

# data dispersion
system.time(I_sp_partitioning(sxyz,adj.w = "weight")) # 0.05 sec.
system.time(I_site_partitioning_simple(sxyz,adj.w = "weight")) # 0.14 sec.
system.time(I_site_partitioning_parall(sxyz,adj.w = "weight",core.use = 0.8)) # 0.20 sec.
system.time(I_genus_sp_partitioning(sxyz,adj.w = "weight")) # 307 sec. 
system.time(I_site_genus_partitioning_1(sxyz,adj.w = "weight")) # 398 sec.
system.time(I_site_genus_partitioning_2(sxyz,adj.w = "weight")) # 392 sec.

# spatial proximity
system.time(D_sp_partitioning(sxyz,rescaled=FALSE,expit=FALSE)) # 151.7 sec.
system.time(D_site_partitioning(sxyz,rescaled=FALSE,expit=FALSE)) # 199.7 sec.
system.time(D_genus_sp_partitioning(sxyz,rescaled = FALSE,expit=FALSE)) # 192.58
system.time(D_site_genus_partitioning_1(sxyz,rescaled = FALSE,expit=FALSE)) # 344.36
system.time(D_site_genus_partitioning_2(sxyz,rescaled = FALSE,expit=FALSE)) # 345.76

################################################################################ ----
# Fig.1                                                                          ----

library(readxl)
sxyz <- read_excel("./input/Data_Figure1.xlsx")
sxyz <- as.matrix(sxyz)

# data dispersion
system.time(I_sp_partitioning(sxyz,adj.w = "weight")) # 0.05 sec.
system.time(I_site_partitioning_simple(sxyz,adj.w = "weight")) # 0.14 sec.
system.time(I_site_partitioning_parall(sxyz,adj.w = "weight",core.use = 0.8)) # 0.20 sec.
system.time(I_genus_sp_partitioning(sxyz,adj.w = "weight")) # 307 sec. 
system.time(I_site_genus_partitioning_1(sxyz,adj.w = "weight")) # 398 sec.
system.time(I_site_genus_partitioning_2(sxyz,adj.w = "weight")) # 392 sec.

# spatial proximity
system.time(D_sp_partitioning(sxyz,rescaled=FALSE,expit=T)) # 151.7 sec.
system.time(D_site_partitioning(sxyz,rescaled=FALSE,expit=T)) # 199.7 sec.
system.time(D_genus_sp_partitioning(sxyz,rescaled = FALSE,expit=T)) # 192.58
system.time(D_site_genus_partitioning_1(sxyz,rescaled = FALSE,expit=T)) # 344.36
system.time(D_site_genus_partitioning_2(sxyz,rescaled = FALSE,expit=T)) # 345.76


################################################################################ ----
# Fig.2                                                                          ----
# Fig.2 A                                                                        
sxyz = data_creation(S=1, G=1, R=8, lambda=2, fix.rd = F)
sxyz[,"abun"] = c(5,0,0,5,5,0,0,5)
I_sp_partitioning(sxyz,adj.w = "not")
I_site_partitioning_simple(sxyz,adj.w = "weight")
D_sp_partitioning(sxyz,rescaled=FALSE, expit=TRUE)
D_site_partitioning(sxyz,rescaled=FALSE, expit=TRUE)



sxyz[,"abun"] = c(4,0,0,6,3,0,0,7)
sxyz
I_sp_partitioning(sxyz,adj.w = "not")
I_site_partitioning_simple(sxyz,adj.w = "not")
D_sp_partitioning(sxyz,rescaled=FALSE,expit=TRUE)
D_site_partitioning(sxyz,rescaled=FALSE,expit=TRUE)

sxyz[,"abun"] = c(9,0,0,1,2,0,0,8)
sxyz
I_sp_partitioning(sxyz,adj.w = "not")
I_site_partitioning_simple(sxyz,adj.w = "not")
D_sp_partitioning(sxyz,rescaled=FALSE,expit=TRUE)
D_site_partitioning(sxyz,rescaled=FALSE,expit=TRUE)


sxyz[,"abun"] = c(5,0,5,5,5,0,0,0)
sxyz
I_sp_partitioning(sxyz,adj.w = "not")
I_site_partitioning_simple(sxyz,adj.w = "not")
D_sp_partitioning(sxyz,rescaled=FALSE,expit=TRUE)
D_site_partitioning(sxyz,rescaled=FALSE,expit=TRUE)

sxyz[,"abun"] = c(4,0,6,7,3,0,0,0)
sxyz
I_sp_partitioning(sxyz,adj.w = "not")
I_site_partitioning_simple(sxyz,adj.w = "not")
D_sp_partitioning(sxyz,rescaled=FALSE,expit=TRUE)
D_site_partitioning(sxyz,rescaled=FALSE,expit=TRUE)

sxyz[,"abun"] = c(9,0,2,8,1,0,0,0)
sxyz
I_sp_partitioning(sxyz,adj.w = "not")
I_site_partitioning_simple(sxyz,adj.w = "not")
D_sp_partitioning(sxyz,rescaled=FALSE,expit=TRUE)
D_site_partitioning(sxyz,rescaled=FALSE,expit=TRUE)

sxyz[,"abun"] = c(5,5,5,5,0,0,0,0)
sxyz
I_sp_partitioning(sxyz,adj.w = "not")
I_site_partitioning_simple(sxyz,adj.w = "not")
D_sp_partitioning(sxyz,rescaled=FALSE,expit=TRUE)
D_site_partitioning(sxyz,rescaled=FALSE,expit=TRUE)

sxyz[,"abun"] = c(4,3,6,7,0,0,0,0)
sxyz
I_sp_partitioning(sxyz,adj.w = "not")
I_site_partitioning_simple(sxyz,adj.w = "not")
D_sp_partitioning(sxyz,rescaled=FALSE,expit=TRUE)
D_site_partitioning(sxyz,rescaled=FALSE,expit=TRUE)

sxyz[,"abun"] = c(9,8,2,1,0,0,0,0)
sxyz
I_sp_partitioning(sxyz,adj.w = "not")
I_site_partitioning_simple(sxyz,adj.w = "not")
D_sp_partitioning(sxyz,rescaled=FALSE,expit=TRUE)
D_site_partitioning(sxyz,rescaled=FALSE,expit=TRUE)

# Fig.2B    
# (1,1)
sxyz = data_creation(S=2, G=1, R=16, lambda=2, fix.rd = F)
sxyz[,"x"] = rep(1:4,4)
sxyz[,"y"] = rep(1:4, times = 2, each = 4)
sxyz[,"abun"] = c(c(5,5,1,0,5,2,0,1,0,0,1,0,1,0,0,0),
                  c(5,5,0,1,5,2,1,0,0,0,0,1,0,0,0,1))

res <- I_sp_partitioning(sxyz,adj.w = "weight")
scales::percent(res$I.within.sp/res$I.total)
scales::percent(res$I.between.sp/res$I.total)

res <- D_sp_partitioning(sxyz,rescaled=FALSE, expit=T)
scales::percent(log(res$D.within.sp)/log(res$D.total))
scales::percent(log(res$D.between.sp)/log(res$D.total))

# (1,2)
sxyz[,"abun"] = c(c(9,1,1,0,5,2,0,1,0,0,1,0,1,0,0,0),
                  c(5,5,0,1,5,2,1,0,0,0,0,1,0,0,0,1))

res <- I_sp_partitioning(sxyz,adj.w = "weight")
scales::percent(res$I.within.sp/res$I.total)
scales::percent(res$I.between.sp/res$I.total)

res <- D_sp_partitioning(sxyz,rescaled=FALSE, expit=T)
scales::percent(log(res$D.within.sp)/log(res$D.total))
scales::percent(log(res$D.between.sp)/log(res$D.total))

# (1,3)
sxyz[,"abun"] = c(c(9,1,1,0,5,2,0,1,0,0,1,0,1,0,0,0),
                  c(9,1,0,1,5,2,1,0,0,0,0,1,0,0,0,1))

res <- I_sp_partitioning(sxyz,adj.w = "weight")
scales::percent(res$I.within.sp/res$I.total)
scales::percent(res$I.between.sp/res$I.total)

res <- D_sp_partitioning(sxyz,rescaled=FALSE, expit=T)
scales::percent(log(res$D.within.sp)/log(res$D.total))
scales::percent(log(res$D.between.sp)/log(res$D.total))

# (2,1)
sxyz[,"abun"] = c(c(5,5,1,1,5,2,1,0,0,0,0,1,0,0,0,0),
                  c(5,5,0,0,5,2,0,1,0,0,1,0,0,0,1,1))

res <- I_sp_partitioning(sxyz,adj.w = "weight")
scales::percent(res$I.within.sp/res$I.total)
scales::percent(res$I.between.sp/res$I.total)

res <- D_sp_partitioning(sxyz,rescaled=FALSE, expit=T)
scales::percent(log(res$D.within.sp)/log(res$D.total))
scales::percent(log(res$D.between.sp)/log(res$D.total))

# (2,2)
sxyz[,"abun"] = c(c(9,1,1,1,5,2,1,0,0,0,0,1,0,0,0,0),
                  c(5,5,0,0,5,2,0,1,0,0,1,0,0,0,1,1))

res <- I_sp_partitioning(sxyz,adj.w = "weight")
scales::percent(res$I.within.sp/res$I.total)
scales::percent(res$I.between.sp/res$I.total)

res <- D_sp_partitioning(sxyz,rescaled=FALSE, expit=T)
scales::percent(log(res$D.within.sp)/log(res$D.total))
scales::percent(log(res$D.between.sp)/log(res$D.total)) 

# (2,3)
sxyz[,"abun"] = c(c(9,1,1,1,5,2,1,0,0,0,0,1,0,0,0,0),
                  c(9,1,0,0,5,2,0,1,0,0,1,0,0,0,1,1))

res <- I_sp_partitioning(sxyz,adj.w = "weight")
scales::percent(res$I.within.sp/res$I.total)
scales::percent(res$I.between.sp/res$I.total)

res <- D_sp_partitioning(sxyz,rescaled=FALSE, expit=T)
scales::percent(log(res$D.within.sp)/log(res$D.total))
scales::percent(log(res$D.between.sp)/log(res$D.total))

# (3,1)
sxyz = data_creation(S=2, G=1, R=16, lambda=2, fix.rd = F)
sxyz[,"x"] = rep(1:4,4)
sxyz[,"y"] = rep(1:4, times = 2, each = 4)
sxyz[,"abun"] = c(c(5,5,1,1,5,2,1,1,0,0,0,0,0,0,0,0),
                  c(5,5,1,1,5,0,1,1,0,0,2,0,0,0,0,0))

res <- I_sp_partitioning(sxyz,adj.w = "weight")
scales::percent(res$I.within.sp/res$I.total)
scales::percent(res$I.between.sp/res$I.total)

res <- D_sp_partitioning(sxyz,rescaled=FALSE, expit=TRUE)
scales::percent(log(res$D.within.sp)/log(res$D.total))
scales::percent(log(res$D.between.sp)/log(res$D.total))

# (3,2)
sxyz = data_creation(S=2, G=1, R=16, lambda=2, fix.rd = F)
sxyz[,"x"] = rep(1:4,4)
sxyz[,"y"] = rep(1:4, times = 2, each = 4)
sxyz[,"abun"] = c(c(9,1,1,1,5,2,1,1,0,0,0,0,0,0,0,0),
                  c(5,5,1,1,5,0,1,1,0,0,2,0,0,0,0,0))

res <- I_sp_partitioning(sxyz,adj.w = "weight")
scales::percent(res$I.within.sp/res$I.total)
scales::percent(res$I.between.sp/res$I.total)

res <- D_sp_partitioning(sxyz,rescaled=FALSE, expit=TRUE)
scales::percent(log(res$D.within.sp)/log(res$D.total))
scales::percent(log(res$D.between.sp)/log(res$D.total))

# (3,3)
sxyz = data_creation(S=2, G=1, R=16, lambda=2, fix.rd = F)
sxyz[,"x"] = rep(1:4,4)
sxyz[,"y"] = rep(1:4, times = 2, each = 4)
sxyz[,"abun"] = c(c(9,1,1,1,5,2,1,1,0,0,0,0,0,0,0,0),
                  c(9,1,1,1,5,0,1,1,0,0,2,0,0,0,0,0))

res <- I_sp_partitioning(sxyz,adj.w = "weight")
scales::percent(res$I.within.sp/res$I.total)
scales::percent(res$I.between.sp/res$I.total)

res <- D_sp_partitioning(sxyz,rescaled=FALSE, expit=TRUE)
scales::percent(log(res$D.within.sp)/log(res$D.total))
scales::percent(log(res$D.between.sp)/log(res$D.total))



################################################################################ ----
# Numeric simulation                                                             ----

library(MASS)
#simulation of species distribution using Poisson cluster process
#generate a matrix with cnum and rnum, species presence in each cell was recorded as 1
#using bivariate probability function to generate offsprings
#Poisson distribution to generate number of parental points
#geometric distribution to generate number of offsprings for each parental point
#uniform distribution to obtain locations of parental points
#cnum determines number of columns, while rnum determines number of rows
#lambda is the parameter for Poisson distribution, prob is the parameter for geometric distribution
species.distribution<-function(delta,lambda,prob,cnum,rnum)
{
  pnum<-rpois(1,lambda)+1 # avoiding zeros, Poission distribution,lambda =  parameter intensity ??
  offnum<-rgeom(pnum,prob)+1 # avoiding zeros,The Geometric Distribution, prob = 0.01
  # generate the locations of parental pointsdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
  x<-sample(cnum,pnum,replace=FALSE)
  y<-sample(rnum,pnum,replace=FALSE)
  #
  #now generating offspring locations
  sig<-matrix(c(delta,0,0,delta),2,2)
  pp<-vector()
  for(i in 1:pnum)
  {
    this<-mvrnorm(offnum[i],c(x[i],y[i]),sig)#bivariate normal distribution,delta link to ??2 
    this<-rbind(this,c(x[i],y[i])) #including parental point is tricky, otherwise sometimes this is NULL
    #exluding out-of-bound points
    idd<-which(this[,1]>cnum | this[,2]>rnum | this[,1]<1 | this[,2]<1)
    if(length(idd)!=0)
    {
      this<-this[-idd,]
    }
    pp<-rbind(pp,this)
  }
  #
  return(pp)
}#end

### example
par(mfrow=c(1,1))
test1 <- species.distribution(delta = 5, lambda = 3, prob = 0.01, cnum = 200, rnum = 200)
plot(test1, xlim = c(0, 200), ylim = c(0, 200), col=2, pch = 1, cex = 0.5, xlab = "Coordinate X",ylab = "Coordinate Y")
abline(v= seq( 0, 200, by=20),h=seq(0, 200,by=20),col = grey(3/8),lty = "dotted")
test2 <- species.distribution(delta = 50, lambda = 3, prob = 0.01, cnum = 200, rnum = 200)
test3 <- species.distribution(delta = 500, lambda = 3, prob = 0.01, cnum = 200, rnum = 200)
test4 <- species.distribution(delta = 5000, lambda = 3, prob = 0.01, cnum = 200, rnum = 200)


points(test2, col =1, pch = 22,cex = 0.5)
points(test3, col =3, pch =4,cex = 0.5)
points(test4, col =5, pch =3,cex = 0.5)


par(mfrow=c(2,2))
library(spatstat) 

xy <- as.ppp(test1, W = c(0,200,0,200))
q1 <- quadratcount(xy, nx=10, ny=10)
plot(q1, main = "")
plot(xy, col=2, pch = 1, cex = 0.5, alpha = 0.05, add = TRUE)

xy <- as.ppp(test2, W = c(0,200,0,200))
q1 <- quadratcount(xy, nx=10, ny=10)
plot(q1, main = "")
plot(xy, col=1, pch = 1, cex = 0.5, add = TRUE)

xy <- as.ppp(test3, W = c(0,200,0,200))
q1 <- quadratcount(xy, nx=10, ny=10)
plot(q1, main = "")
plot(xy, col=3, pch = 1, cex = 0.5, add = TRUE)

xy <- as.ppp(test4, W = c(0,200,0,200))
q1 <- quadratcount(xy, nx=10, ny=10)
plot(q1, main = "")
plot(xy, col=5, pch = 1, cex = 0.5, add = TRUE)

Mul_SDM_site <- function(x_length = 20, # Indicates the length of each sample plot (in meters)
                         y_length = 20, # Indicates the width of each sample plot (in meters)
                         delta = 5,     # Parameter delta is determined by the variance-covariance matrix of the distribution and controls the dispersion parameter of the data points.
                         lambda = 3,    # lambda is the parameter for Poisson distribution
                         prob = 0.01,   # prob is the parameter for geometric distribution
                         cnum = 200,    # cnum determines number of columns  (Length of entire study area in meters)
                         rnum = 200,    # while rnum determines number of rows  (Width of entire study area in meters)
                         Spnum = c(5,4,4,3,2), # Indicates the number of species per genus
                         Gnum = 5       # How many genera are there?
){
  library(spatstat) 
  options(warn = -1)
  # 判断样地数量
  R = (cnum/x_length)*(rnum/y_length)
  
  # 初始化属和物种序号矩阵
  SP_result <- matrix(nrow = sum(Spnum), ncol = 2)
  
  species_counter <- 1 # 用于连续计数的物种序号
  
  for (u in 1:Gnum) {
    for (spp in 1:Spnum[u]) {
      SP_result[species_counter, 1] <- u # 属序号
      SP_result[species_counter, 2] <- species_counter # 物种序号
      species_counter <- species_counter + 1
    }
  }
  
  # 新建sxyz数据框
  sxyz <- matrix(ncol = 6)
  colnames(sxyz) = c("genus","sp","x","y","site", "abun")
  sxyz = sxyz[-1,]
  
  library(stringr)
  
  # 定义函数以提取范围的下限和上限
  extract_limits <- function(range_str) {
    str_extract_all(range_str, "\\d+")[[1]]
  }
  
  
  plot(c(0,0),xlim = c(0, cnum), ylim = c(0, rnum), col=0, pch = 1, cex = 0.5, xlab = "Coordinate X",ylab = "Coordinate Y")
  
  # 对于每一个属循环
  for (u in 1:nrow(SP_result)) {
    dat = species.distribution(delta = delta, lambda = lambda, prob = prob, cnum = cnum, rnum = rnum)
    points(dat, col = u, pch = 22,cex = 0.5)
    xy_ppp <- as.ppp(dat, W = c(0,cnum,0,rnum))
    q1 <- quadratcount(xy_ppp, nx=cnum/x_length, ny=rnum/y_length)
    q1 <- as.data.frame(q1)
    names(q1) = c("or_y","or_x","abun")
    q1$x  <- (sapply(q1$or_x, function(x) as.numeric(extract_limits(x)[1])) +  sapply(q1$or_x, function(x) as.numeric(extract_limits(x)[2])))/2000
    q1$y  <- (sapply(q1$or_y, function(y) as.numeric(extract_limits(y)[1])) +  sapply(q1$or_y, function(y) as.numeric(extract_limits(y)[2])))/2000
    q1$genus = rep(SP_result[u,1],nrow(q1))
    q1$sp = rep(SP_result[u,2],nrow(q1))
    q1$site = 1:nrow(q1)
    sxyz = rbind(sxyz, q1[,c("genus","sp","x","y","site", "abun")])
  }
  
  return(as.matrix(sxyz))
}

# Test 
Gnum <- 1
min_val <- 2 
max_val <- 5 
Spnum <- sample(min_val:max_val, Gnum , replace = TRUE)


sxyz = 
  Mul_SDM_site(x_length = 20, 
               y_length = 20,
               delta = 50,     
               lambda = 3,  
               prob = 0.01, 
               cnum = 200,   
               rnum = 200,    
               Spnum = Spnum, 
               Gnum = Gnum)       


# data dispersion
system.time(I_sp_partitioning(sxyz,adj.w = "weight")) # 0.05 sec.
system.time(I_site_partitioning_simple(sxyz,adj.w = "weight")) # 0.14 sec.
system.time(I_site_partitioning_parall(sxyz,adj.w = "weight",core.use = 0.8)) # 0.20 sec.
system.time(I_genus_sp_partitioning(sxyz,adj.w = "weight")) # 307 sec. 
system.time(I_site_genus_partitioning_1(sxyz,adj.w = "weight")) # 398 sec.
system.time(I_site_genus_partitioning_2(sxyz,adj.w = "weight")) # 392 sec.

# spatial proximity
system.time(D_sp_partitioning(sxyz,rescaled=FALSE,expit=FALSE)) # 151.7 sec.
system.time(D_site_partitioning(sxyz,rescaled=FALSE,expit=FALSE)) # 199.7 sec.
system.time(D_genus_sp_partitioning(sxyz,rescaled = FALSE,expit=FALSE)) # 192.58
system.time(D_site_genus_partitioning_1(sxyz,rescaled = FALSE,expit=FALSE)) # 344.36
system.time(D_site_genus_partitioning_2(sxyz,rescaled = FALSE,expit=FALSE)) # 345.76


