library(tidyverse)

# Read in the orthogroup table
ortho <- read_tsv("results/N0.tsv")

# pivot to long format.
# if value is NA (i.e. not present in that genome), remove it
ortho_long <- ortho %>%
  pivot_longer(cols = -c("HOG", "OG", "Gene Tree Parent Clade"), names_to = "genome", values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  separate_rows(gene, sep = ",")

# calculate the number of genomes that each OG is present in and add to the table
ortho_long <- ortho_long %>%
  group_by(OG) %>%
  mutate(n_genomes = n_distinct(genome)) %>%
  ungroup()

# get NLRtracker results
genomes <- list.files("NLRtracker")

nlrtracker <- list()

for (i in genomes) {
  nlrtracker[[i]] <- read_tsv(paste0("NLRtracker/", i, "/", i, "_NLRtracker.tsv"))
  nlrtracker[[i]]$genome <- genomes[i]
  nlrtracker[[i]] <- nlrtracker[[i]] %>%
    rename(gene = "seqname", subclass = "Subclass (putative)")
}

nlrtracker <- bind_rows(nlrtracker, .id = "genome")

# join by gene and genome
ortho_long <- ortho_long %>%
  left_join(nlrtracker, by = c("gene", "genome"))

# if subclass is NA then remove it
ortho_long <- ortho_long %>%
  filter(!is.na(subclass))

# calculate number of each subclass for each og
ortho_long <- ortho_long %>%
  group_by(OG, subclass) %>%
  summarise(n = n_distinct(gene)) %>%
  ungroup()
