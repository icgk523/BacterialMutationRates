require(phytools)
require(readxl)
require(nlme)
require(tidyverse)
require(evobiR)
require(phyr)

final_tree <- read.tree("20161115_pruned_tree_cutoff_0_03.nwk")
enz_oi <- c("uvrA", "nfi", "mutT", "mutY", "mutL", "mutS", "mutH", "uvrD", "dnaQ", "ung", "mutM")

enz_transitions <- sapply(enz_oi,function(x) NULL)

for(i in 1:length(enz_oi)) {
  enz_transitions[[enz_oi[i]]] <- read.csv(paste0("enz_transitions/20180206_",enz_oi[i],"_transitions.csv"),header = T,row.names = 1)
}

orig_pa_matrix <- read.csv("phylo_full_table.csv", header = T, row.names = NULL)

#### changing the tree tip labels to species names ####
final_tree$tip.label <- orig_pa_matrix$ncbi_name[match(final_tree$tip.label, orig_pa_matrix$phylotree_tip_label)]

#### filter based on species list ####
complete_data <- read.csv("mutation_rate_data.csv")
species_list <- complete_data$Species


for (species_name in final_tree$tip.label){  # Loop through each species name in the tree
  partial_species_name <- paste(strsplit(species_name, " ")[[1]][1:2], collapse = " ")  # Get the first two words of the species name
  if (!any(grepl(partial_species_name, species_list, ignore.case = TRUE))) {  # Check if the partial species name is not in the species list (case insensitive)
    final_tree <- drop.tip(final_tree, species_name, trim.internal = TRUE)  # Remove the species from the tree if not found in the list
  }
}

#### Collapse strains into one and then randomly sample any species which are not directly next to each other ####
collapse_identical_tips <- function(phy,tip_label){
  matching_tips <- which(sapply(strsplit(phy$tip.label, " "), function(x) paste(x[1:2], collapse = " ")) == tip_label)
  nt <- length(phy$tip.label) # Number of tips in the tree
  nm <- length(matching_tips) # Number of tips matching the label
  keep <- numeric(nm)# Initialize a vector to keep track of which matching tips to keep
  
  cur_tip <- 1 # Start from the first matching tip
  while(cur_tip<=nm){ # Iterate over all matching tips
    if(cur_tip == nm){ # If it's the last tip, mark it to keep and break the loop
      keep[cur_tip] <- 1
      break
    }
    next_tip <- cur_tip + 1
    mrca_ <- getMRCA(phy,c(matching_tips[cur_tip],matching_tips[next_tip])) # Find the most recent common ancestor (MRCA) of the current tip and the next tip
    descendants <- getDescendants(phy, mrca_)  # Find all descendants of the MRCA
    descendant_tips <- descendants[descendants<=nt]  # Identify which of these descendants are tips (leaf nodes)
    if(all(descendant_tips %in% matching_tips)){ # If all descendant tips are in the list of matching tips
      keep[cur_tip] <- 1  # Mark the current tip to be kept
      cur_tip <- cur_tip + length(descendant_tips) # Move the current tip index forward by the number of descendant tips
    }else{
      keep[cur_tip] <- 1  # Mark the current tip to be kept
      cur_tip <- cur_tip + 1 # Move to the next tip
    }
  }
  to_drop <- matching_tips[!keep] # Identify the tips that are not marked to be kept
  new_phy <- drop.tip(phy,to_drop) # Drop these tips from the phylogeny
  
  return(new_phy) 
}
collapse_all_species <- function(phy) {
  species <- unique(sapply(strsplit(phy$tip.label, " "), function(x) paste(x[1:2], collapse = " "))) # Extract unique first two words as species names
  for(species_name in species) {
    if(sum(grepl(species_name, sapply(strsplit(phy$tip.label, " "), function(x) paste(x[1:2], collapse = " ")))) > 1) { # Check if the species name occurs more than once
      phy <- collapse_identical_tips(phy, species_name) # Collapse tips for the current species
    }
  }
  return(phy)
}

final_tree_collapsed <- collapse_all_species(final_tree)

# for (i in 1:20){ # Plot the tree to see differences across 20 iterations
#   final_tree_collapsed = collapse_all_species(final_tree_collapsed)
#   plot(final_tree_collapsed)
# }

remove_final_duplicates <- function(phy) {
  species_counts <- table(sapply(strsplit(phy$tip.label, " "), function(x) paste(x[1], x[2], sep = " "))) # Count occurrences of each species based on the first two words of the tip labels
  species_to_remove <- names(species_counts)[species_counts > 1]  # Species with more than one presence
  for(species_name in species_to_remove) { # Iterate over each species that has duplicates
    removed_species <- sample(phy$tip.label[grepl(species_name, phy$tip.label)], 1) # Randomly select one tip label from the duplicates to remove
    phy <- drop.tip(phy, removed_species) # Remove the selected tip from the phylogeny
  }
  return(phy)
}

#### Merge the tree with the data set ####
final_tree_collapsed$tip.label <- gsub("^([^[:space:]]+\\s+[^[:space:]]+).*", "\\1", final_tree_collapsed$tip.label)
final_tree_collapsed$tip.label <- gsub("[^[:alpha:][:space:]]", "", final_tree_collapsed$tip.label)
final_tree_collapsed <- keep.tip(final_tree_collapsed, species_list[species_list %in% final_tree_collapsed$tip.label]) 

# Subset trees by the number of strains
two_strain_tree <- keep.tip(final_tree_collapsed, complete_data$Species[complete_data$Strains >= 2][complete_data$Species[complete_data$Strains >= 2] %in% final_tree_collapsed$tip.label])
three_strain_tree <- keep.tip(final_tree_collapsed, complete_data$Species[complete_data$Strains > 2][complete_data$Species[complete_data$Strains == 2] %in% final_tree_collapsed$tip.label])

# Subset tree by those with polymorphisms
two_strain_tree <- keep.tip(two_strain_tree, complete_data$Species[complete_data$Overall_Rate > 0][complete_data$Species[complete_data$Overall_Rate > 0] %in% two_strain_tree$tip.label])
three_strain_tree <- keep.tip(three_strain_tree, complete_data$Species[complete_data$Overall_Rate > 0][complete_data$Species[complete_data$Overall_Rate >0] %in% three_strain_tree$tip.label])

# Plot the trees
plot(two_strain_tree, show.tip.label = T, edge.width = 1, cex = 0.5)
plot(three_strain_tree, show.tip.label = T, edge.width = 1, cex = 0.5)

# Save the trees
write.tree(two_strain_tree, file = "two_strain_tree.nwk")
write.tree(three_strain_tree, file = "three_strain_tree.nwk")


#### Watterson's Analysis ####
subset_data <- complete_data %>% filter(Wattersons_Corrected_Overall_Rate != "#DIV/0!") %>% mutate(Log_Wattersons_Corrected_Overall_Rate = log10(as.numeric(Wattersons_Corrected_Overall_Rate))) # Keep only Watterson's data which has values

two_strain_tree <- keep.tip(two_strain_tree, subset_data$Species[subset_data$Species %in% two_strain_tree$tip.label]) # Subset the tree to the species present in the data
subset_data = subset_data[subset_data$Species %in% two_strain_tree$tip.label, ] # Subset the data to the species in the tree for redundancy
subset_data$Log_Wattersons_Corrected_Overall_Rate <- as.numeric(subset_data$Log_Wattersons_Corrected_Overall_Rate) # Make sure Watterson's is numeric

rownames(subset_data) = subset_data$Species # Set rownames for reordering to phylogeny
reordered_subset_data <- ReorderData(two_strain_tree, subset_data, taxa.names="row names") # Reorder the data to be the same as the order in the phylogeny; 'row names' reflects the where the species names are
reordered_subset_data$Log_Wattersons_Corrected_Overall_Rate <- as.numeric(reordered_subset_data$Log_Wattersons_Corrected_Overall_Rate) # Keep it numeric

genes = names(subset_data[4:14]) # names of genes we're looking at
results_watterson <- data.frame(gene = character(), lambda = numeric(), linear_AIC = numeric(), pglmm_AIC = numeric(), linear_p_value = numeric(), pglmm_p_value = numeric(), stringsAsFactors = FALSE) # An empty data frame to store outputs


for (gene in genes) {
  gene_data <- reordered_subset_data[[gene]]  # Extract the column corresponding to the current gene
  linear_model <- lm(Log_Wattersons_Corrected_Overall_Rate ~ gene_data, data = reordered_subset_data)
  
  residuals_test <- residuals(linear_model)  # Get residuals and compute lambda
  estimated_lambda <- phylosig(two_strain_tree, residuals_test, method = "lambda")[[1]]
  
  # Linear model statistics
  linear_AIC <- AIC(linear_model)
  linear_p_value <- coef(summary(linear_model))[8]
  
  # Add linear model results to the table
  results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "linear", AIC = linear_AIC, p_value = linear_p_value, stringsAsFactors = FALSE))
  
  # Linear Null
  linear_null_model <- lm(Log_Wattersons_Corrected_Overall_Rate ~ 1, data = reordered_subset_data)
  linear_null_AIC <- AIC(linear_null_model)
  linear_null_p_value <- coef(summary(linear_null_model))[8]
  
  # Add linear null model results to the table
  results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "linear_null", AIC = linear_null_AIC, p_value = NA, stringsAsFactors = FALSE))
  
  
  # Fit PGLMM model
  pglmm_model <- pglmm(Log_Wattersons_Corrected_Overall_Rate ~ gene_data + (1 | Species), data = reordered_subset_data, family = "gaussian", cov_ranef = list(Species = two_strain_tree), REML = FALSE)
  pglmm_AIC <- pglmm_model$AIC
  pglmm_p_value <- pglmm_model$B.pvalue[[2]]
  
  # Add PGLMM model results to the table
  results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "pglmm", AIC = pglmm_AIC, p_value = pglmm_p_value, stringsAsFactors = FALSE))
  
  # Adaptive lambda starting value search
  if (estimated_lambda <= 0.5) {
    search_lambdas <- seq(0, 1, by = 0.1)
  } else {
    search_lambdas <- rev(seq(0, 1, by = 0.1))
  }
  
  # Fit models with different correlation structures and their null models
  for (lambda_starting in search_lambdas) {
    pagel_model <- tryCatch({    # Try corPagel
      gls(Log_Wattersons_Corrected_Overall_Rate ~ gene_data, data = reordered_subset_data, correlation = corPagel(lambda_starting, two_strain_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    pagel_null_model <- tryCatch({    # Try null model with corPagel
      gls(Log_Wattersons_Corrected_Overall_Rate ~ 1, data = reordered_subset_data, correlation = corPagel(lambda_starting, two_strain_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(pagel_model)) {
      pagel_AIC <- AIC(pagel_model)
      pagel_p_value <- coef(summary(pagel_model))[8]
      results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "pagel", AIC = pagel_AIC, p_value = pagel_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(pagel_null_model)) {
      pagel_null_AIC <- AIC(pagel_null_model)
      results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "pagel_null", AIC = pagel_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    blomberg_model <- tryCatch({    # Try corBlomberg
      gls(Log_Wattersons_Corrected_Overall_Rate ~ gene_data, data = reordered_subset_data, correlation = corBlomberg(lambda_starting, two_strain_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    blomberg_null_model <- tryCatch({    # Try null model with corBlomberg
      gls(Log_Wattersons_Corrected_Overall_Rate ~ 1, data = reordered_subset_data, correlation = corBlomberg(lambda_starting, two_strain_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(blomberg_model)) {
      blomberg_AIC <- AIC(blomberg_model)
      blomberg_p_value <- coef(summary(blomberg_model))[8]
      results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "blomberg", AIC = blomberg_AIC, p_value = blomberg_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(blomberg_null_model)) {
      blomberg_null_AIC <- AIC(blomberg_null_model)
      results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "blomberg_null", AIC = blomberg_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    brownian_model <- tryCatch({    # Try corBrownian
      gls(Log_Wattersons_Corrected_Overall_Rate ~ gene_data, data = reordered_subset_data, correlation = corBrownian(lambda_starting, two_strain_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    brownian_null_model <- tryCatch({    # Try null model with corBrownian
      gls(Log_Wattersons_Corrected_Overall_Rate ~ 1, data = reordered_subset_data, correlation = corBrownian(lambda_starting, two_strain_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(brownian_model)) {
      brownian_AIC <- AIC(brownian_model)
      brownian_p_value <- coef(summary(brownian_model))[8]
      results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "brownian", AIC = brownian_AIC, p_value = brownian_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(brownian_null_model)) {
      brownian_null_AIC <- AIC(brownian_null_model)
      results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "brownian_null", AIC = brownian_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    martins_model <- tryCatch({    # Try corMartins
      gls(Log_Wattersons_Corrected_Overall_Rate ~ gene_data, data = reordered_subset_data, correlation = corMartins(lambda_starting, two_strain_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    martins_null_model <- tryCatch({    # Try null model with corMartins
      gls(Log_Wattersons_Corrected_Overall_Rate ~ 1, data = reordered_subset_data, correlation = corMartins(lambda_starting, two_strain_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(martins_model)) {
      martins_AIC <- AIC(martins_model)
      martins_p_value <- coef(summary(martins_model))[8]
      results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "martins", AIC = martins_AIC, p_value = martins_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(martins_null_model)) {
      martins_null_AIC <- AIC(martins_null_model)
      results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "martins_null", AIC = martins_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    grafen_model <- tryCatch({    # Try corGrafen
      gls(Log_Wattersons_Corrected_Overall_Rate ~ gene_data, data = reordered_subset_data, correlation = corGrafen(lambda_starting, two_strain_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    grafen_null_model <- tryCatch({    # Try null model with corGrafen
      gls(Log_Wattersons_Corrected_Overall_Rate ~ 1, data = reordered_subset_data, correlation = corGrafen(lambda_starting, two_strain_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(grafen_model)) {
      grafen_AIC <- AIC(grafen_model)
      grafen_p_value <- coef(summary(grafen_model))[8]
      results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "grafen", AIC = grafen_AIC, p_value = grafen_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(grafen_null_model)) {
      grafen_null_AIC <- AIC(grafen_null_model)
      results_watterson <- rbind(results_watterson, data.frame(gene = gene, lambda = estimated_lambda, model = "grafen_null", AIC = grafen_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
}


#### GC & Ts/Tv Analysis ####
complete_data = complete_data %>% 
  mutate(GC_rate = log10((TC + TG + AG + AC)/totalAT/((CT + CA + GA + GT)/totalGC + (TC + TG + AG + AC)/totalAT)),
         tvts_ratio = log10((CT + TC + GA + AG)/(CA + CG + TG + TA + GT + GC + AC + AT)))
subset_data <- complete_data %>% 
  filter(is.finite(GC_rate) & is.finite(tvts_ratio))

gc_subset = subset_data %>% filter(Strains > 2)
tvts_subset = subset_data %>% filter(Strains >= 2)

gc_tree <- keep.tip(two_strain_tree, gc_subset$Species[gc_subset$Species %in% two_strain_tree$tip.label]) # Subset the tree to the species present in the data
tvts_tree <- keep.tip(two_strain_tree, tvts_subset$Species[tvts_subset$Species %in% two_strain_tree$tip.label]) # Subset the tree to the species present in the data

gc_subset = gc_subset[gc_subset$Species %in% gc_tree$tip.label, ] # Subset the data to the species in the tree for redundancy
tvts_subset = tvts_subset[tvts_subset$Species %in% tvts_tree$tip.label, ] # Subset the data to the species in the tree for redundancy

rownames(gc_subset) = gc_subset$Species # Set rownames for reordering to phylogeny
gc_subset <- ReorderData(gc_tree, gc_subset, taxa.names="row names") # Reorder the data to be the same as the order in the phylogeny; 'row names' reflects the where the species names are

rownames(tvts_subset) = tvts_subset$Species # Set rownames for reordering to phylogeny
tvts_subset <- ReorderData(tvts_tree, tvts_subset, taxa.names="row names") # Reorder the data to be the same as the order in the phylogeny; 'row names' reflects the where the species names are

genes_gc = names(gc_subset[4:14] %>%   select_if(~ length(unique(.)) >= 2 & !all(. == 0) & !all(. == 1)))
genes_tvts = names(tvts_subset[4:14] %>%   select_if(~ length(unique(.)) >= 2 & !all(. == 0) & !all(. == 1)))

results_gc_rate <- data.frame(gene = character(), lambda = numeric(), model = character(), AIC = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
results_tvts_ratio <- data.frame(gene = character(), lambda = numeric(), model = character(), AIC = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each gene
for (gene in genes_gc) {
  gene_data <- gc_subset[[gene]]  # Extract the column corresponding to the current gene
  linear_model <- lm(GC_rate ~ gene_data, data = gc_subset)
  
  residuals_test <- residuals(linear_model)  # Get residuals and compute lambda
  estimated_lambda <- phylosig(gc_tree, residuals_test, method = "lambda")[[1]]
  
  # Linear model statistics
  linear_AIC <- AIC(linear_model)
  linear_p_value <- coef(summary(linear_model))[8]
  
  # Add linear model results to the table
  results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "linear", AIC = linear_AIC, p_value = linear_p_value, stringsAsFactors = FALSE))
  
  # Linear Null
  linear_null_model <- lm(GC_rate ~ 1, data = gc_subset)
  linear_null_AIC <- AIC(linear_null_model)
  linear_null_p_value <- coef(summary(linear_null_model))[8]
  
  # Add linear null model results to the table
  results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "linear_null", AIC = linear_null_AIC, p_value = NA, stringsAsFactors = FALSE))
  
  # Fit PGLMM model
  pglmm_model <- pglmm(GC_rate ~ gene_data + (1 | Species), data = gc_subset, family = "gaussian", cov_ranef = list(Species = gc_tree), REML = FALSE)
  pglmm_AIC <- pglmm_model$AIC
  pglmm_p_value <- pglmm_model$B.pvalue[[2]]
  
  # Add PGLMM model results to the table
  results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "pglmm", AIC = pglmm_AIC, p_value = pglmm_p_value, stringsAsFactors = FALSE))
  
  # Adaptive lambda starting value search
  if (estimated_lambda <= 0.5) {
    search_lambdas <- seq(0, 1, by = 0.1)
  } else {
    search_lambdas <- rev(seq(0, 1, by = 0.1))
  }
  
  # Fit models with different correlation structures and their null models
  for (lambda_starting in search_lambdas) {
    pagel_model <- tryCatch({    # Try corPagel
      gls(GC_rate ~ gene_data, data = gc_subset, correlation = corPagel(lambda_starting, gc_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    pagel_null_model <- tryCatch({    # Try null model with corPagel
      gls(GC_rate ~ 1, data = gc_subset, correlation = corPagel(lambda_starting, gc_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(pagel_model)) {
      pagel_AIC <- AIC(pagel_model)
      pagel_p_value <- coef(summary(pagel_model))[8]
      results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "pagel", AIC = pagel_AIC, p_value = pagel_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(pagel_null_model)) {
      pagel_null_AIC <- AIC(pagel_null_model)
      results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "pagel_null", AIC = pagel_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    blomberg_model <- tryCatch({    # Try corBlomberg
      gls(GC_rate ~ gene_data, data = gc_subset, correlation = corBlomberg(lambda_starting, gc_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    blomberg_null_model <- tryCatch({    # Try null model with corBlomberg
      gls(GC_rate ~ 1, data = gc_subset, correlation = corBlomberg(lambda_starting, gc_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(blomberg_model)) {
      blomberg_AIC <- AIC(blomberg_model)
      blomberg_p_value <- coef(summary(blomberg_model))[8]
      results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "blomberg", AIC = blomberg_AIC, p_value = blomberg_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(blomberg_null_model)) {
      blomberg_null_AIC <- AIC(blomberg_null_model)
      results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "blomberg_null", AIC = blomberg_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
    
  for (lambda_starting in search_lambdas) {
    brownian_model <- tryCatch({    # Try corBrownian
      gls(GC_rate ~ gene_data, data = gc_subset, correlation = corBrownian(lambda_starting, gc_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    brownian_null_model <- tryCatch({    # Try null model with corBrownian
      gls(GC_rate ~ 1, data = gc_subset, correlation = corBrownian(lambda_starting, gc_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(brownian_model)) {
      brownian_AIC <- AIC(brownian_model)
      brownian_p_value <- coef(summary(brownian_model))[8]
      results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "brownian", AIC = brownian_AIC, p_value = brownian_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(brownian_null_model)) {
      brownian_null_AIC <- AIC(brownian_null_model)
      results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "brownian_null", AIC = brownian_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    martins_model <- tryCatch({    # Try corMartins
      gls(GC_rate ~ gene_data, data = gc_subset, correlation = corMartins(lambda_starting, gc_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    martins_null_model <- tryCatch({    # Try null model with corMartins
      gls(GC_rate ~ 1, data = gc_subset, correlation = corMartins(lambda_starting, gc_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(martins_model)) {
      martins_AIC <- AIC(martins_model)
      martins_p_value <- coef(summary(martins_model))[8]
      results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "martins", AIC = martins_AIC, p_value = martins_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(martins_null_model)) {
      martins_null_AIC <- AIC(martins_null_model)
      results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "martins_null", AIC = martins_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    grafen_model <- tryCatch({    # Try corGrafen
      gls(GC_rate ~ gene_data, data = gc_subset, correlation = corGrafen(lambda_starting, gc_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    grafen_null_model <- tryCatch({    # Try null model with corGrafen
      gls(GC_rate ~ 1, data = gc_subset, correlation = corGrafen(lambda_starting, gc_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(grafen_model)) {
      grafen_AIC <- AIC(grafen_model)
      grafen_p_value <- coef(summary(grafen_model))[8]
      results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "grafen", AIC = grafen_AIC, p_value = grafen_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(grafen_null_model)) {
      grafen_null_AIC <- AIC(grafen_null_model)
      results_gc_rate <- rbind(results_gc_rate, data.frame(gene = gene, lambda = estimated_lambda, model = "grafen_null", AIC = grafen_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
}

# Loop through each gene
for (gene in genes_tvts) {
  gene_data <- tvts_subset[[gene]]  # Extract the column corresponding to the current gene
  linear_model <- lm(tvts_ratio ~ gene_data, data = tvts_subset)
  
  residuals_test <- residuals(linear_model)  # Get residuals and compute lambda
  estimated_lambda <- phylosig(tvts_tree, residuals_test, method = "lambda")[[1]]
  
  # Linear model statistics
  linear_AIC <- AIC(linear_model)
  linear_p_value <- coef(summary(linear_model))[8]
  
  # Add linear model results to the table
  results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "linear", AIC = linear_AIC, p_value = linear_p_value, stringsAsFactors = FALSE))
  
  # Linear Null
  linear_null_model <- lm(tvts_ratio ~ 1, data = tvts_subset)
  linear_null_AIC <- AIC(linear_null_model)
  linear_null_p_value <- coef(summary(linear_null_model))[8]
  
  # Add linear null model results to the table
  results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "linear_null", AIC = linear_null_AIC, p_value = NA, stringsAsFactors = FALSE))
  
  # Fit PGLMM model
  pglmm_model <- pglmm(tvts_ratio ~ gene_data + (1 | Species), data = tvts_subset, family = "gaussian", cov_ranef = list(Species = tvts_tree), REML = FALSE)
  pglmm_AIC <- pglmm_model$AIC
  pglmm_p_value <- pglmm_model$B.pvalue[[2]]
  
  # Add PGLMM model results to the table
  results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "pglmm", AIC = pglmm_AIC, p_value = pglmm_p_value, stringsAsFactors = FALSE))
  
  # Adaptive lambda starting value search
  if (estimated_lambda <= 0.5) {
    search_lambdas <- seq(0, 1, by = 0.1)
  } else {
    search_lambdas <- rev(seq(0, 1, by = 0.1))
  }
  
  # Fit models with different correlation structures and their null models
  for (lambda_starting in search_lambdas) {
    pagel_model <- tryCatch({    # Try corPagel
      gls(tvts_ratio ~ gene_data, data = tvts_subset, correlation = corPagel(lambda_starting, tvts_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    pagel_null_model <- tryCatch({    # Try null model with corPagel
      gls(tvts_ratio ~ 1, data = tvts_subset, correlation = corPagel(lambda_starting, tvts_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(pagel_model)) {
      pagel_AIC <- AIC(pagel_model)
      pagel_p_value <- coef(summary(pagel_model))[8]
      results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "pagel", AIC = pagel_AIC, p_value = pagel_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(pagel_null_model)) {
      pagel_null_AIC <- AIC(pagel_null_model)
      results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "pagel_null", AIC = pagel_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    blomberg_model <- tryCatch({    # Try corBlomberg
      gls(tvts_ratio ~ gene_data, data = tvts_subset, correlation = corBlomberg(lambda_starting, tvts_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    blomberg_null_model <- tryCatch({    # Try null model with corBlomberg
      gls(tvts_ratio ~ 1, data = tvts_subset, correlation = corBlomberg(lambda_starting, tvts_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(blomberg_model)) {
      blomberg_AIC <- AIC(blomberg_model)
      blomberg_p_value <- coef(summary(blomberg_model))[8]
      results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "blomberg", AIC = blomberg_AIC, p_value = blomberg_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(blomberg_null_model)) {
      blomberg_null_AIC <- AIC(blomberg_null_model)
      results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "blomberg_null", AIC = blomberg_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    brownian_model <- tryCatch({    # Try corBrownian
      gls(tvts_ratio ~ gene_data, data = tvts_subset, correlation = corBrownian(lambda_starting, tvts_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    brownian_null_model <- tryCatch({    # Try null model with corBrownian
      gls(tvts_ratio ~ 1, data = tvts_subset, correlation = corBrownian(lambda_starting, tvts_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(brownian_model)) {
      brownian_AIC <- AIC(brownian_model)
      brownian_p_value <- coef(summary(brownian_model))[8]
      results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "brownian", AIC = brownian_AIC, p_value = brownian_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(brownian_null_model)) {
      brownian_null_AIC <- AIC(brownian_null_model)
      results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "brownian_null", AIC = brownian_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    martins_model <- tryCatch({    # Try corMartins
      gls(tvts_ratio ~ gene_data, data = tvts_subset, correlation = corMartins(lambda_starting, tvts_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    martins_null_model <- tryCatch({    # Try null model with corMartins
      gls(tvts_ratio ~ 1, data = tvts_subset, correlation = corMartins(lambda_starting, tvts_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(martins_model)) {
      martins_AIC <- AIC(martins_model)
      martins_p_value <- coef(summary(martins_model))[8]
      results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "martins", AIC = martins_AIC, p_value = martins_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(martins_null_model)) {
      martins_null_AIC <- AIC(martins_null_model)
      results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "martins_null", AIC = martins_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
  
  for (lambda_starting in search_lambdas) {
    grafen_model <- tryCatch({    # Try corGrafen
      gls(tvts_ratio ~ gene_data, data = tvts_subset, correlation = corGrafen(lambda_starting, tvts_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    grafen_null_model <- tryCatch({    # Try null model with corGrafen
      gls(tvts_ratio ~ 1, data = tvts_subset, correlation = corGrafen(lambda_starting, tvts_tree, form = ~Species), method = "ML")
    }, error = function(e) NULL)
    
    if (!is.null(grafen_model)) {
      grafen_AIC <- AIC(grafen_model)
      grafen_p_value <- coef(summary(grafen_model))[8]
      results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "grafen", AIC = grafen_AIC, p_value = grafen_p_value, stringsAsFactors = FALSE))
    }
    
    if (!is.null(grafen_null_model)) {
      grafen_null_AIC <- AIC(grafen_null_model)
      results_tvts_ratio <- rbind(results_tvts_ratio, data.frame(gene = gene, lambda = estimated_lambda, model = "grafen_null", AIC = grafen_null_AIC, p_value = NA, stringsAsFactors = FALSE))
      break
    }
  }
}

# Sort the results by gene alphabetically and then by AIC in ascending order
results_tvts_ratio <- results_tvts_ratio %>%  arrange(gene, AIC)
results_gc_rate <- results_gc_rate %>%  arrange(gene, AIC)

summary_table_wattersons = data.frame(gene = character(), lambda = numeric(), Mean_with_Gene = numeric(), SE_with_Gene = numeric(), Mean_without_Gene = numeric(), SE_without_Gene = numeric(), Ratio_w_wo = numeric(), stringsAsFactors = FALSE)
for(gene in genes_tvts){
  lambda <- phylosig(two_strain_tree, residuals(lm(Log_Wattersons_Corrected_Overall_Rate ~ reordered_subset_data[[gene]], data = reordered_subset_data)), method = "lambda")[[1]]
  Mean_with_Gene = mean(reordered_subset_data$Log_Wattersons_Corrected_Overall_Rate[reordered_subset_data[[gene]]==1])
  Mean_without_Gene = mean(reordered_subset_data$Log_Wattersons_Corrected_Overall_Rate[reordered_subset_data[[gene]]==0])
  SE_with_Gene = sd(reordered_subset_data$Log_Wattersons_Corrected_Overall_Rate[reordered_subset_data[[gene]] == 1], na.rm = TRUE) / sqrt(sum(reordered_subset_data[[gene]] == 1, na.rm = TRUE))
  SE_without_Gene = sd(reordered_subset_data$Log_Wattersons_Corrected_Overall_Rate[reordered_subset_data[[gene]] == 0], na.rm = TRUE) / sqrt(sum(reordered_subset_data[[gene]] == 0, na.rm = TRUE))
  Ratio_w_wo = sum(reordered_subset_data[[gene]] == 0, na.rm = TRUE) / sum(reordered_subset_data[[gene]] == 1, na.rm = TRUE)
  summary_table_wattersons = rbind(summary_table_wattersons, data.frame(gene = gene, lambda = lambda, Mean_with_Gene = Mean_with_Gene, Mean_without_Gene = Mean_without_Gene, Ratio_w_wo = Ratio_w_wo))
}

summary_table_tvts = data.frame(gene = character(), lambda = numeric(), Mean_with_Gene = numeric(), SE_with_Gene = numeric(), Mean_without_Gene = numeric(), SE_without_Gene = numeric(), Ratio_w_wo = numeric(), stringsAsFactors = FALSE)
for(gene in genes_tvts){
  lambda <- phylosig(tvts_tree, residuals(lm(tvts_ratio ~ tvts_subset[[gene]], data = tvts_subset)), method = "lambda")[[1]]
  Mean_with_Gene = mean(tvts_subset$tvts_ratio[tvts_subset[[gene]]==1])
  Mean_without_Gene = mean(tvts_subset$tvts_ratio[tvts_subset[[gene]]==0])
  SE_with_Gene = sd(tvts_subset$tvts_ratio[tvts_subset[[gene]] == 1], na.rm = TRUE) / sqrt(sum(tvts_subset[[gene]] == 1, na.rm = TRUE))
  SE_without_Gene = sd(tvts_subset$tvts_ratio[tvts_subset[[gene]] == 0], na.rm = TRUE) / sqrt(sum(tvts_subset[[gene]] == 0, na.rm = TRUE))
  Ratio_w_wo = sum(tvts_subset[[gene]] == 0, na.rm = TRUE) / sum(tvts_subset[[gene]] == 1, na.rm = TRUE)
  summary_table_tvts = rbind(summary_table_tvts, data.frame(gene = gene, lambda = lambda, Mean_with_Gene = Mean_with_Gene, Mean_without_Gene = Mean_without_Gene, Ratio_w_wo = Ratio_w_wo))
}
summary_table_gc = data.frame(gene = character(), lambda = numeric(), Mean_with_Gene = numeric(), SE_with_Gene = numeric(), Mean_without_Gene = numeric(), SE_without_Gene = numeric(), Ratio_w_wo = numeric(), stringsAsFactors = FALSE)
for(gene in genes_gc){
  lambda <- phylosig(gc_tree, residuals(lm(GC_rate ~ gc_subset[[gene]], data = gc_subset)), method = "lambda")[[1]]
  Mean_with_Gene = mean(gc_subset$GC_rate[gc_subset[[gene]] == 1], na.rm = TRUE)
  Mean_without_Gene = mean(gc_subset$GC_rate[gc_subset[[gene]] == 0], na.rm = TRUE)
  SE_with_Gene = sd(gc_subset$GC_rate[gc_subset[[gene]] == 1], na.rm = TRUE) / sqrt(sum(gc_subset[[gene]] == 1, na.rm = TRUE))
  SE_without_Gene = sd(gc_subset$GC_rate[gc_subset[[gene]] == 0], na.rm = TRUE) / sqrt(sum(gc_subset[[gene]] == 0, na.rm = TRUE))
  Ratio_w_wo = sum(gc_subset[[gene]] == 0, na.rm = TRUE) / sum(gc_subset[[gene]] == 1, na.rm = TRUE)
  summary_table_gc = rbind(summary_table_gc, data.frame(gene = gene, lambda = lambda, Mean_with_Gene = Mean_with_Gene, SE_with_Gene = SE_with_Gene, Mean_without_Gene = Mean_without_Gene, SE_without_Gene = SE_without_Gene, Ratio_w_wo = Ratio_w_wo))
}

cat("The models with the best fits for log(Watterson's Diversity) are:\n"); data.frame(results_watterson %>% group_by(gene) %>% slice(which.min(AIC)) %>% select(gene, lambda, model, AIC, p_value))
cat("The summary table for log(Watterson's Diversity) is:\n"); print(summary_table_wattersons)

cat("The models with the best fits for log(Ts:Tv) are:\n"); data.frame(results_tvts_ratio %>% group_by(gene) %>% slice(which.min(AIC)) %>% select(gene, lambda, model, AIC, p_value))
cat("The summary table for log(Ts:Tv) is:\n"); print(summary_table_tvts)

cat("The models with the best fitsfor log(GC > AT) are:\n"); data.frame(results_gc_rate %>% group_by(gene) %>% slice(which.min(AIC)) %>% select(gene, lambda, model, AIC, p_value))
cat("The summary table for log(GC > AT) is:\n"); print(summary_table_gc)

cat("\n For all AICs that have not been included in the above summaries, look at the files: results_watterson, results_gc_rate, results_tvts_ratio")



