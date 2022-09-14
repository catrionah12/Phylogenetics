install.packages("rstudioapi")
install.packages("stringr")
install.packages("tidyr")
install.packages("dplyr")
install.packages("ggplot2")
library(rstudioapi)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)

setwd(dirname(getActiveDocumentContext()$path))
setwd('../genomes/data/proteins/function_annotations')

## Pathway heatmap
pathway <- 'map00910'  ## CHANGE THIS TO CHOOSE PATHWAY ##

#Read in annotation files
a99_annotations <-  read.csv("a99_annotations.tsv", sep = '\t')
vulgaris_annotations <- read.csv("c_vulgaris_annotations.tsv", sep = '\t')
conductrix_annotations <- read.csv("m_conductrix_annotations.tsv", sep = '\t')
variabilis_annotations <- read.csv("c_variabilis_annotations.tsv", sep = '\t')
sorokiniana_annotations <- read.csv("c_sorokiniana_annotations.tsv", sep = '\t')
annotations <- rbind(a99_annotations, vulgaris_annotations, conductrix_annotations, variabilis_annotations,sorokiniana_annotations)
annotations <- data.frame(annotations$X.query, annotations$EC, annotations$KEGG_Pathway, annotations$BRITE)

#Read in file for orthogroups in each species 
orthogroups <- read.csv('../../../OrthoFinder/Results_Aug28/Orthogroups/Orthogroups.tsv', sep = '\t')

# Label OGs from each species with KEGG pathways and EC numbers
annotated_ogs = data.frame()
for (i in 2:6) {
  og_genes <- data.frame(orthogroups[ , 1], orthogroups[ , i])
  separated_genes <- separate_rows(og_genes, orthogroups...i., sep = ',')
  trimmed_genes <- separated_genes %>% mutate(across(where(is.character), str_trim))
  complete_og_genes <- trimmed_genes[!trimmed_genes$orthogroups...i.=="", ]
  colnames(complete_og_genes) <- c('orthogroup', 'protein_name')
  annotated_og_species <- merge(complete_og_genes, annotations,  by.x = 'protein_name', by.y = 'annotations.X.query')
  annotated_ogs <- rbind(annotated_ogs, annotated_og_species[-c(1)])
}
annotated_ogs = unique(annotated_ogs)

colnames(annotated_ogs) = c('Orthogroup', 'EC', 'KEGG_pathway', 'KEGG_BRITE')

# Get orthogroups annotated with chosen pathway
pathway_ogs <- data.frame(annotated_ogs$Orthogroup[str_detect(annotated_ogs$KEGG_pathway, pathway)],
                         annotated_ogs$EC[str_detect(annotated_ogs$KEGG_pathway, pathway)])
colnames(pathway_ogs) <- c("Orthogroup", "EC")
pathway_ogs <- pathway_ogs %>% arrange(EC)

#Read in number of genes per species for each orthogroup
gene_count <- read.csv('../../../OrthoFinder/Results_Aug28/Orthogroups/Orthogroups.GeneCount.tsv', sep = '\t')

#Get number of genes for each orthogroup in the pathway
pathway_counts <- merge(pathway_ogs, gene_count, by = 'Orthogroup')
pathway_counts <- pathway_counts[-c(8,2)]
colnames(pathway_counts) <- c('Orthogroup', 'a99', 'c_sorokiniana', 'c_variabilis', 'c_vulgaris', 'm_conductrix')

#Make heatmap of counts
heatmap <- cbind(pathway_counts[1], stack(pathway_counts[2:6]))
colnames(heatmap) <- c('Orthogroup', 'Count', 'Species')
positions <- c("m_conductrix","a99","c_variabilis", "c_vulgaris", "c_sorokiniana")
p = ggplot(heatmap, aes(y = Orthogroup, x = Species , fill = Count)) + geom_tile() + scale_x_discrete(limits = positions) + scale_fill_viridis_c(option = 'D', limits = c(0,4)) 
p + scale_y_discrete(limits = pathway_ogs$Orthogroup)


## Get ids of one-to-one orthologues in pathway for dn/ds calculation

setwd('../../../evolution_rates/ids')

one_to_one <- data.frame(pathway_counts$Orthogroup[pathway_counts$a99 ==1 & pathway_counts$c_sorokiniana
                                                  ==1 & pathway_counts$c_variabilis == 1
                                                  & pathway_counts$c_vulgaris == 1
                                                  & pathway_counts$m_conductrix ==1])
one_to_one <- cbind(one_to_one, 'm_conductrix', 'A99', 'c_variabilis', 'c_vulgaris', 'c_sorokiniana')
colnames(one_to_one) <- c('Orthogroup', 'm_conductrix', 'A99', 'c_variabilis', 'c_vulgaris', 'c_sorokiniana')

for (i in one_to_one$Orthogroup)
{
  pathway_ogs = as.character(orthogroups[orthogroups$Orthogroup==i,2:6])
  writeLines(pathway_ogs, paste(i, '.txt', sep = ''))
}
