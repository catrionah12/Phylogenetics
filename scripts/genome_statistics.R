install.packages("rstudioapi")
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

## Genome statistics
setwd('../genomes/data/stats')

order <- c('C_sorokiniana', 'M_conductrix', 'A99', 'C_variabilis', 'C_vulgaris')
stats <- read.csv('genome_stats.txt', sep = '\t')
q <- ggplot(stats, aes(x = Species, y = Genome_size)) + geom_bar(stat = 'identity', width = 0.8, , fill = '#F5D440') + coord_flip() + scale_x_discrete(limits = rev(order))
q + theme_light() + labs(y = "Genome size (Mbp)")

protein_coding <- data.frame(stats$Species, stats$Protein_coding, 'zzz')
pseudo <- data.frame(stats$Species, stats$Pseudogenes, 'Pseudogenes')
colnames(protein_coding) <- c('Species', 'Number', 'Type')
colnames(pseudo) <- c('Species', 'Number', 'Type')
genes_chart <- rbind(pseudo, protein_coding)

p = ggplot(genes_chart, aes(x = Species, y = Number, group = Type, fill = Type)) + geom_bar(position = 'dodge', stat = 'identity', width = 0.8) + coord_flip() + scale_x_discrete(limits = rev(order)) 
p +theme_light() + labs(y = 'Number of genes (thousand)') + scale_fill_manual(values = c('#21918c', '#3b528b'), labels = c('Pseudogenes', 'Protein-coding genes'))


## Genes-per-orthogroup-per-species


duplication_stats <- read.csv("duplication_statistics.txt", sep = '\t')
genes_per_species = cbind(duplication_stats[1], stack(duplication_stats[2:6]))
colnames(genes_per_species) = c('genes_per_species', 'Percent', "Species")
p =ggplot(genes_per_species, aes(x= genes_per_species, y = Percent, group = Species, color = Species)) + geom_line() + labs(y = '% of orthogroups', x = "Number of genes in orthogroup") + xlim(0,8)
p
