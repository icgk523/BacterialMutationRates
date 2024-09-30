# Project Summary
* **Authors:** George Kalogiannis (g.kalogiannis23@imperial.ac.uk) and Adam Eyre-Walker (a.c.eyre-walker@sussex.ac.uk)
* **Date:** August 2024

## Table of Contents
* [Code](#code)
* [Data](#data)
* [Results](#results)


## Code
Code directory contains the following R scripts. Please inspect each file before sourcing, as each may require various packages to be installed:
- **mass_gene_analysis.R:** A script for conducting data and phylogeny preparation and statistical analysis of bacterial mutation patterns and rates.
- **rsquared_gls_function.R:** Functions for calculating R<sup>2</sup>. These are fixed functions from the 'piecewiseSEM' package from Lefcheck (2016; Methods Ecol. Evol.)

## Data
The data directory contains the following subdirectories and files:
- **20161115_pruned_tree_cutoff_0_03.nwk:** Phylogenetic tree of bacterial strains obtained from Sane et al. (2023).
- ```enz_transitions/```:  Files of bacterial repair gene states (present or absent) at phylogeny nodes.
- **phylo_full_table.csv:** File used for mapping NCBI accession number of phylogeny tip label to NCBI strain name.
- **mutation_rate_data.csv:** File containing mutation information estimated for bacterial species.

## Results
The results directory contains the following:
- **results_[gc/tvts/watterson].csv:** Results from phylogenetic analyses of traits against individual repair gene presence and absence. 
- **results_[gc/tvts/watterson]_all_genes.csv:** Results from phylogenetic analyses of traits against the total number of repair genes.
- **summary_table_[gc/tvts/watterson].csv:** Summaries of mean and standard errors of results for each repair gene and measured trait.
- **two_strain_tree.nwk:** Newick phylogenetic tree of bacterial species with a number of strains greater or equal to 2 per species. 
- **three_strain_tree.nwk:** Newick phylogenetic tree of bacterial species with at least 2 strains per species.
- **gene_presence_phylogeny.pdf:** Phylogenetic tree of repair gene presence and absence for 240 bacterial species depicted in the two strain phylogenetic tree.