# Phylogenetics

## Requirements

To run the pipelines, ETE Toolkit and entrez-direct must be installed.

This can be done using conda:

    conda create -n phylogeny python=3.6

    conda install -c bioconda entrez-direct

    conda install -c etetoolkit ete3 ete_toolchain

## Phylogenetic trees

To produce a Maximum Likelihood tree of free-living and symbiotic Chlorella species, move into the 'scripts' directory and run:

    ./download_seqs.sh && ./alignment.sh && ./ml_tree.sh

To produce a more in-depth Maximum Likelihood tree of symbiotic Chlorella, instead run:

    ./download_seqs.sh sym && ./alignment.sh && ml_tree.sh

- Results will be saved to phylogeny/trees

- The alignments in phylogeny/temp can be used to produce Bayesian trees using BEAST

## KEGG pathway analysis

### Heatmap of genes-per-orthogroup in a pathway 

- In pathways.R, the 'pathway' variable can be changed to choose a KEGG pathway

- pathways.R will produce a heatmap of the number of genes-per-species in each orthogroup associated with the pathway

- EC numbers for each orthogroup in the pathway can be found in the 'pathway_ogs' variable

- The gene names of single-copy orthologues in the pathway will be saved for use in the dN/dS calculation

### dN/dS ratio calculation for single-copy orthologues in a pathway

To compare the free-branch model to the M0 model with *C. sp.* A99 set as the foreground branch, run:

    ./dnds.sh

If the free-branch model has a significantly higher likelihood, a tree of the calculated dN/dS ratios will be saved to genomes/evolution_rates/results.

The foreground branch can be changed; run:

    ./dnds.sh mcon    
to set the foreground branch to *M. conductrix*, or: 

    ./dnds.sh cvar
to set the foreground branch to *C. variabilis*.

## General genome statistics

- genome_statistics.R can be used to produce graphs showing general statistics about the 5 *Chlorella* genomes
    


