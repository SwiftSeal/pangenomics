library(ggtree)
library(gggenes)
library(ggplot2)
library(cowplot)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)

# Configuration

genomes <- c(
    "annuum" = "C. annuum",
    "candolleanum" = "S. candolleanum",
    "etuberosum" = "S. etuberosum",
    "lycopersicum" = "S. lycopersicum",
    "tuberosum_phureja_E4_63" = "S. tuberosum\n(phureja E4-63)",
    "tuberosum_tuberosum_RH10_15" = "S. tuberosum\n(tuberosum RH10-15)",
    "chilense" = "S. chilense",
    "galapagense" = "S. galapagense",
    "neorickii" = "S. neorickii",
    "tuberosum_phureja_E86_69" = "S. tuberosum\n(phureja E86-69)",
    "tuberosum_tuberosum_RH" = "S. tuberosum\n(tuberosum RH)",
    "chmielewskii" = "S. chmielewskii",
    "habrochaites" = "S. habrochaites",
    "peruvianum" = "S. peruvianum",
    "tuberosum_stenotomum_A6_26" = "S. tuberosum\n(stenotomum A6-26)",
    "verrucosum" = "S. verrucosum",
    "corneliomulleri" = "S. corneliomulleri",
    "lycopersicoides" = "S. lycopersicoides",
    "pimpinellifolium" = "S. pimpinellifolium",
    "tuberosum_stenotomum_PG6359" = "S. tuberosum\n(stenotomum PG6359)"
)

tuberising <- c(
    "annuum" = FALSE,
    "candolleanum" = TRUE,
    "etuberosum" = FALSE,
    "lycopersicum" = FALSE,
    "tuberosum_phureja_E4_63" = TRUE,
    "tuberosum_tuberosum_RH10_15" = TRUE,
    "chilense" = FALSE,
    "galapagense" = FALSE,
    "neorickii" = FALSE,
    "tuberosum_phureja_E86_69" = TRUE,
    "tuberosum_tuberosum_RH" = TRUE,
    "chmielewskii" = FALSE,
    "habrochaites" = FALSE,
    "peruvianum" = FALSE,
    "tuberosum_stenotomum_A6_26" = TRUE,
    "verrucosum" = TRUE,
    "corneliomulleri" = FALSE,
    "lycopersicoides" = FALSE,
    "pimpinellifolium" = FALSE,
    "tuberosum_stenotomum_PG6359" = TRUE
)

paired_colours = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")
nice_colours = c("#2596be", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

nice_theme <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank(),

  strip.text = element_blank())

# Import data

## Merge all resistify results into a single dataframe
resistify <- map(names(genomes), ~ read.table(paste0("results/resistify/", .x, "/results.tsv"), header = TRUE) %>% mutate(genome = .x)) %>%
  bind_rows()

## Orthofinder tree
genome_tree <- read.tree("results/orthofinder/Species_Tree/SpeciesTree_rooted.txt")

## EDTA calculated % of TE content
te_content <- map(names(genomes), ~ read.table(paste0("results/edta/", .x, ".mod.EDTA.TEanno.sum"), skip = 6, fill = TRUE, nrows = 12, col.names = c("te_classification", "te_count", "te_masked_bp", "percentage")) %>% filter(te_count != "--") %>%  mutate(genome = .x, percentage = as.numeric(str_remove(percentage, "%")))) %>%
  bind_rows()

## Orthofinder orthogroups
orthogroups <- read.table("results/orthofinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", sep = "\t", header = TRUE) %>%
  pivot_longer(!1:3, names_to = "genome", values_to = "gene") %>%
  separate_rows(gene, sep = ", ")

## Helixer gff
helixer_gff <- map(names(genomes), ~ read.table(paste0("results/helixer/", .x, ".gff")) %>% select(c(V3,V9)) %>% filter(V3 == "exon") %>% mutate(genome = .x)) %>%
  bind_rows()

helixer_summary <- helixer_gff %>%
  mutate(gene_id = str_extract(V9, "Parent=(([^;]+))")) %>%
  mutate(gene_id = str_replace(gene_id, "Parent=", "")) %>%
  group_by(genome, gene_id) %>%
  summarise(count = n())

ggplot(helixer_summary, aes(y = count, fill = genome)) +
  geom_boxplot()

## Overlaps

te_overlap <- map(names(genomes), ~ read.table(paste0("results/intersect/", .x, ".bed")) %>% mutate(genome = .x)) %>%
  bind_rows()

# Analysis

summarise(group_by(resistify, genome), count = n()) %>%
  mutate(tuberising = tuberising[match(genome, names(tuberising))]) %>%
  ggplot(aes(x = tuberising, y = count)) +
  geom_boxplot(fill = nice_colours[1]) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black"))

association_plot <- inner_join(
           te_content,
           summarise(group_by(resistify, genome), count = n()),
           by = "genome") %>%
  filter(te_classification != "unknown") %>%
  mutate(tuberising = tuberising[match(genome, names(tuberising))]) %>%
  ggplot(aes(x = percentage, y = count, colour = tuberising)) +
  geom_point() +
  facet_wrap(~ te_classification, scales = "free_x") +
  theme(
    legend.position = c(0.84, 0.15),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_blank()) +
  labs(x = "TE content (%)", y = "Number of resistance genes") +
  scale_colour_manual(values = nice_colours)

ggsave("association_plot.png", association_plot, width = 4, height = 4, units = "in")

## NLR orthogroups

orthogroup_histogram_plot <- orthogroups %>%
  filter(gene %in% resistify$Sequence) %>%
  group_by(genome, OG) %>%
  summarise(count = n()) %>%
  group_by(OG) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = count)) +
  geom_bar(fill = nice_colours[1]) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Number of genomes", y = "Number of orthogroups")

ggsave("orthogroup_histogram.png", orthogroup_histogram_plot, width = 2, height = 2, units = "in")

orthogroup_classification <- resistify %>%
  left_join(orthogroups, by = c("Sequence" = "gene")) %>%
  group_by(Classification) %>%
  summarise(count = n_distinct(OG)) %>%
  ggplot(aes(x = Classification, y = count, fill = Classification)) +
  geom_bar(stat = "identity") +
    theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "none") +
  labs(x = "Classification", y = "Number of orthogroups") +
  scale_fill_manual(values = paired_colours)

ggsave("orthogroup_classification.png", orthogroup_classification, width = 2, height = 2, units = "in")

# of the CNLs, how many are True for MADA?
#resistify %>%
#    filter(Classification == "CNL") %>%
#    summarise(count = sum(MADA == "True"), total = n())

# Phylogenetic tree

## Summarise resistify results

resistify_summary <- resistify %>%
    group_by(genome, Classification) %>%
    summarise(count = n())

## create column of nice names
resistify_summary$GenomeFormatted <- genomes[match(resistify_summary$genome, names(genomes))]

## Plot tree

tree_plot <- ggtree(genome_tree) +
    geom_tiplab(align = TRUE, color = "white") +
    theme(plot.margin = unit(c(0,0,0.35,0), "in"))
order <- get_taxa_name(tree_plot)
resistify_summary$genome <- factor(resistify_summary$genome, levels = order)

## plot facet of resistify summaries
resistify_plot <- ggplot(resistify_summary %>% mutate(dummy = "dummy"), aes(y = dummy, x = count, fill = Classification)) +
    geom_bar(stat = "identity") +
    facet_grid(genome ~ .) +
    theme(legend.position = c(0.75, 0.8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()
    ) +
    scale_fill_manual(values = paired_colours)

# plot a facet of just the genome names
names_plot <- ggplot(resistify_summary) +
    geom_text(aes(x = 0, y = 0, label = GenomeFormatted, fontface = "italic"), size = 2, hjust = 0) +
    xlim(0,1) +
    facet_grid(genome ~ .) +
    theme(plot.margin = unit(c(0.06,0,0.4,0), "in"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()
    )

merged_plot <- plot_grid(tree_plot,
                         names_plot,
                         resistify_plot,
                         ncol = 3,
                         rel_widths = c(1, 0.7, 3))

ggsave("results/plot.png", merged_plot, width = 6, height = 6, units = "in")

## Overlap analysis

te_overlap %>%
  filter(V4 %in% resistify$Sequence) %>%
  ggplot(aes(y = genome, fill = V15)) +
    geom_bar()

