# Project Summary
* **Author:** George Kalogiannis (g.kalogiannis23@imperial.ac.uk), Adam Eyre-Walker (a.c.eyre-walker@sussex.ac.uk)
* **Date:** August 2024

## Table of Contents
* [Code](#code)
* [Data](#data)
* [Results](#results)


## Code
Code directory contains the following R scripts. Please inspect each file before sourcing, as each may require various packages to be installed:
- **mass_gene_analysis.R:** A script for analysing the relationship between mass & genome size in insects. The script uses phylogenetic data, mass and genome information from the directories: ```data/phylogeny```, ```data/mass```, and ```data/genome```. 
- **toga_output_analysis.R:** A script for performing model comparisons on patterns of gene loss relative to body mass in insects. The script inputs the outputs from the Tool to Infer Orthologs from Genome Alignments (TOGA) pipeline from Kirilenko et al. (2023; Science), and uses phylogenetic and mass information to estimate trends. Inputs data from ```data/phylogeny```, ```data/mass```, ```data/toga```.
- **hox_gene_analysis.R:** A script for performing model comparisons on patterns of Hox gene loss relative to body mass in insects. The script inputs the outputs from the HbxFinder pipeline from Mulhair et al. (2023; Genome Res.), and uses phylogenetic and mass information to estimate trends. Inputs data from ```data/phylogeny```, ```data/mass```, ```data/hox```.
- **transposable_element_analysis.R:**  A script for performing model comparisons on patterns of proportions of transposable elements in an insect species genome relative to its body mass. Inputs data from ```data/phylogeny```, ```data/mass```, ```data/repetitive_elements```.
- **rsquared_gls_function.R:** Functions for calculating R<sup>2</sup>. These are fixed functions from the 'piecewiseSEM' package from Lefcheck (2016; Methods Ecol. Evol.)

## Data
The data directory contains the following subdirectories and files:
- ```genome/```:  NCBI information on whole-genome sequence accession lengths.
- ```mass/```: Contains two files outputted from the [InsectMasses](https://github.com/icgk523/InsectMasses.git) GitHub page with information on >6000 insect species masses, and a subset with genomic information.
- ```phylogeny/```: Contains a Newick tree of 645 insect species. Time-uncalibrated tree was downloaded for the species using the 'rotl' R package (Michonneau et al., 2016; Methods Ecol. Evol.) and calibrated using a time-calibrated phylogeny from 'TimeTree' (Kumar et al., 2022; Mol. Biol. Evol.) and the 'congruification' function from Eastman et al. (2013; Methods Ecol. Evol.).
- ```toga/```: Contains four summary tables of orthology metrics for different insect orders, outputted by TOGA.
- ```hox/```: Contains numbers of individual Hox genes and their identity reports outputted from HbxFinder pipeline for four orders of insects. Also contains the N50 scores for the genomes analysed.

## Results
An empty directory used as a place to output files from the scripts run in code. While the scripts do not necessarily produce outputs here, users including code that does this should specify here as an output location.

## HPC
- **download_genomes.sh**: A script for downloading genomes from the NCBI genomes ftp site from a list of genome annotations.
- **make_lastz_chains_job_submission.sh**: An example job script for running the 'make_lastz_chains' genome alignment pipeline from Kirilenko et al. (2023; Science) on the Imperial College London high-performance computer with Hymenoptera genomes.
- **toga_job_submission.sh**: An example job script for running the 'TOGA' orthology annotation pipeline from Kirilenko et al. (2023; Science) on the Imperial College London high-performance computer with Lepidoptera genomes.
- **output_analysis.py**: A python script for extracting the orthology information from the TOGA '.log' files.

## Proposal
- **proposal.tex:** A LaTex document containing the proposal. It creates a Bibliography using the ```bibliography.bib``` file.
- **bibliography.bib:** A Bibliography used when compiling the ```proposal.tex``` file.
- **run_proposal.sh:** A bash shell script for compiling ```proposal.tex``` through the command line. Running the script produces the output pdf file ```proposal.pdf```.
- **proposal.pdf:** The output from compiling the ```proposal.tex``` file using the ```proposal.sh``` bash script.
- **gantt.sty:** A LaTex package for creating Gantt charts.

## Project
- **project.tex:** A LaTex document containing the proposal. It creates a Bibliography using the ```bibliography.bib``` file.
- **bibliography.bib:** A Bibliography used when compiling the ```project.tex``` file.
- **run_project.sh:** A bash shell script for compiling ```project.tex``` through the command line. Running the script produces the output pdf file ```project.pdf```.
- **project.pdf:** The output from compiling the ```project.tex``` file using the ```project.sh``` bash script.
- ```figures/```: Logo and figures used to compile the ```project.pdf``` report.