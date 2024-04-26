# First, load necessary libraries
library(ape)
library(phangorn)
library('TreeDist')
library('tidyverse')

# Generate all the simulated trees we care about

command = paste0("iqtree2 --alisim alignment_bd -t -rlen 0.04 0.08 0.16",
            "RANDOM{bd{0.1/0.05}/10} -m JC --length 10 -seed 1")

args <- commandArgs(trailingOnly = TRUE)
print(args)  # This will print all the arguments

bound <- as.numeric(args[1])

seeds = (bound-1):(bound)

print(seeds)

greaters = c(16, 32, 64, 128, 256, 512, 1024, 2048)
models = c("JC") # , "F81", "K80", "HKY") # Jukes Cantor, Felsenstion, Kimura, HKY,

lengths = c(4:10, greaters)

# Function for loading a simulation
load_simulation <- function(model, length, seed) {
    prefix <- paste0("simulations/m" , model, "-l", length, "-s", seed)
    tree_filename <- paste0(prefix, ".treefile")
    phy_filename <- paste0(prefix, ".phy")
  
    tree <- read.tree(tree_filename)

    sequences <- read.phyDat(phy_filename, format = "interleaved", type = "DNA")

    return(list(tree = tree, sequences = sequences))
}

get_trees <- function(sequences) {
    # Compute distance matrix
    dna_dist <- dist.hamming(sequences)
    UPGMA_tree <- upgma(dna_dist)
    
    # Get most parsimonious
    parsimonious_tree <- optim.parsimony(UPGMA_tree, sequences, method="sankoff", rearrangements="NNI")
    
    # Get best likelihood
    fit <- pml(UPGMA_tree, sequences)
    fitML <- optim.pml(fit, model = "JC", rearrangements = "NNI", optNni = TRUE)
    likelihood_tree <- fitML$tree
    
    return(list(parsimonious_tree = parsimonious_tree, likelihood_tree = likelihood_tree))
}


evaluate_both <- function(optimal, sequences) {
    # Get trees
    trees <- get_trees(sequences)
    
    rf_ptree <- RF.dist(optimal, trees$parsimonious_tree, check.labels = TRUE, normalize = TRUE)
    rf_ltree <- RF.dist(optimal, trees$likelihood_tree, check.labels = TRUE, normalize = TRUE)
    
    nni_ptree <- NNIDist(optimal, trees$parsimonious_tree)[['best_upper']]
    nni_ltree <- NNIDist(optimal, trees$likelihood_tree)[['best_upper']]
    
    return(list(rf_ptree = rf_ptree, rf_ltree = rf_ltree, nni_ptree = nni_ptree, nni_ltree = nni_ltree))
    
}


# Evaluation function 1.0: naive distance measure
evaluate_RF <- function(optimal, sequences) {
    # Get trees
    trees <- get_trees(sequences)
    
    dist_ptree <- RF.dist(optimal, trees$parsimonious_tree, check.labels = TRUE, normalize = TRUE)
    dist_ltree <- RF.dist(optimal, trees$likelihood_tree, check.labels = TRUE, normalize = TRUE)
    
    return(list(dist_ptree = dist_ptree, dist_ltree = dist_ltree))
}

# Evaluation function 2 using NNID
evaluate_NNID <- function(optimal, sequences) {
    # Get trees
    trees <- get_trees(sequences)
    
    dist_ptree <- NNIDist(optimal, trees$parsimonious_tree)[['best_upper']]
    dist_ltree <- NNIDist(optimal, trees$likelihood_tree)[['best_upper']]
        
    return(list(dist_ptree = dist_ptree, dist_ltree = dist_ltree))
}

results <- tibble(length = integer(),
                  seed = integer(),
                  rf_ptree = double(),
                  rf_ltree = double(),
                  nni_ptree = double(),
                  nni_ltree = double())

# Loop over lengths and seeds
for (length in lengths) {
    for (seed in seeds) {
        # Load the simulation and evaluate the RF
        result <- load_simulation("JC", length = length, seed = seed)
        dist_metrics <- evaluate_both(result$tree, result$sequences)
        
        # Add the results to the tibble
        results <- results %>%
          add_row(length = length,
                  seed = seed,
                  rf_ptree = dist_metrics$rf_ptree,
                  rf_ltree = dist_metrics$rf_ltree,
                  nni_ptree = dist_metrics$nni_ptree,
                  nni_ltree = dist_metrics$nni_ltree)
    }
}

savename = paste0("results", args[1], ".csv")

write.csv(results, savename, row.names = FALSE)
