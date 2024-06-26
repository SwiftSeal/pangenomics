library(ggtree)
library(ggplot2)
#library(ggseqlogo)
library(cowplot)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(GenomicRanges)

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

refplantnlr_genomes <- c(
  "Ty-2" = "habrochaites",
  "Ptr1" = "lycopersicoides",
  "Bs4" = "lycopersicum",
  "CNL-11990" = "lycopersicum",
  "Hero-A" = "lycopersicum",
  "I2" = "lycopersicum",
  "Mi-1" = "lycopersicum",
  "NRC0" = "lycopersicum",
  "NRC1" = "lycopersicum",
  "NRC2" = "lycopersicum",
  "NRC3" = "lycopersicum",
  "NRC4a" = "lycopersicum",
  "NRC4b" = "lycopersicum",
  "NRC6" = "lycopersicum",
  "Rx4" = "lycopersicum",
  "Tm2" = "lycopersicum",
  "tm2_sus" = "lycopersicum",
  "Tm2-nv" = "lycopersicum",
  "Tm2^2" = "lycopersicum",
  "Prf" = "lycopersicum",
  "Sw5-b" = "lycopersicum",
  "Ph-3" = "pimpinellifolium",
  "SpNBS-LRR" = "pimpinellifolium",
  "Gpa2" = "tuberosum_tuberosum_RH10_15",
  "Gro1-4" = "tuberosum_tuberosum_RH10_15",
  "R3a" = "tuberosum_tuberosum_RH10_15",
  "Rx" = "tuberosum_tuberosum_RH10_15",
  "StrPtr1" = "tuberosum_tuberosum_RH10_15"
)

transposable_elements <- c(
"CACTA_TIR_transposon" = "TIR",
"Copia_LTR_retrotransposon" = "LTR",
"Gypsy_LTR_retrotransposon" = "LTR",
"hAT_TIR_transposon" = "TIR",
"helitron" = "Helitron",
"LTR_retrotransposon" = "LTR",
"Mutator_TIR_transposon" = "TIR",
"PIF_Harbinger_TIR_transposon" = "TIR",
"Tc1_Mariner_TIR_transposon" = "TIR"
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

## EDTA calculated % of TE contento
te_content <- map(names(genomes), ~ read.table(paste0("results/edta/", .x, ".mod.EDTA.TEanno.sum"), skip = 6, fill = TRUE, nrows = 12, col.names = c("te_classification", "te_count", "te_masked_bp", "percentage")) %>% filter(te_count != "--") %>%  mutate(genome = .x, percentage = as.numeric(str_remove(percentage, "%")))) %>%
  bind_rows()

## Orthofinder orthogroups
orthogroups <- read.table("results/orthofinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", sep = "\t", header = TRUE) %>%
  pivot_longer(!1:3, names_to = "genome", values_to = "gene") %>%
  separate_rows(gene, sep = ", ")

## Refplantnlr blast results

refplantnlr <- map(names(genomes), ~ read.table(paste0("results/refplantnlr/", .x, ".txt"), sep = "\t") %>% mutate(genome = .x)) %>%
  bind_rows()

## Helixer bed

helixer_bed <- map(names(genomes), ~ read.table(paste0("results/helixer/", .x, ".bed")) %>% mutate(genome = .x)) %>%
  bind_rows()

## genome size

genome_size <- data.frame(genome = character(), Total_Length = numeric(), stringsAsFactors = FALSE)
for (genome in names(genomes)) {
  # Read the fourth line from the file
  total_length_line <- readLines(paste0("results/edta/", genome, ".mod.EDTA.TEanno.sum"), n = 4)[4]

  # Extract the numeric value using regular expression
  total_length <- as.numeric(regmatches(total_length_line, regexpr("\\d+", total_length_line)))

  # Add the results to the data frame
  genome_size <- rbind(genome_size, data.frame(genome = genome, Total_Length = total_length, stringsAsFactors = FALSE))
}

## Overlaps

te_overlap <- map(names(genomes), ~ read.table(paste0("results/intersect/", .x, ".bed")) %>% mutate(genome = .x)) %>%
  bind_rows()

te_overlap <- te_overlap %>%
  mutate(te_id = str_extract(V21, "Name=(([^;]+))", group = 1)) %>%
  mutate(te_type = V15, gene = V4) %>%
  select(genome, te_id, te_type, gene)

## TE terminal motifs

te_terminal_motifs <- map(names(genomes), ~ read.table(paste0("results/edta/", .x, "_intact_terminal.tsv"), sep = "\t", comment.char = "") %>% mutate(genome = .x)) %>%
  bind_rows()

te_terminal_motifs <- te_terminal_motifs %>%
  mutate(te_id = str_remove(V1, "\\|.*"))

## TE intact gff

te_intact_gff <- map(names(genomes), ~ read.table(paste0("results/edta/", .x, ".mod.EDTA.intact.gff3"), sep = "\t", comment.char = "") %>% mutate(genome = .x)) %>%
  bind_rows()

# Analysis

genome_summary <- helixer_bed %>%
  group_by(genome) %>%
  summarise(count = n()) %>%
  left_join(genome_size, by = "genome")

genome_summary <- te_content %>%
  group_by(genome) %>%
  summarise(te_percentage = sum(percentage)) %>%
  left_join(genome_summary, by = "genome")

genome_summary <- te_intact_gff %>%
  filter(V3 %in% names(transposable_elements)) %>%
  group_by(genome) %>%
  summarise(intact_te_count = n()) %>%
  left_join(genome_summary, by = "genome")

genome_summary <- genome_summary %>%
  mutate(tuberising = tuberising[match(genome, names(tuberising))])

write.table(genome_summary, "results/genome_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

summary(lm(intact_te_count ~ tuberising, data = genome_summary))

summarise(group_by(resistify, genome), count = n()) %>%
  mutate(tuberising = tuberising[match(genome, names(tuberising))]) %>%
  ggplot(aes(x = tuberising, y = count)) +
  geom_boxplot(fill = nice_colours[1]) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black"))

te_association_plot <- inner_join(
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

ggsave("results/te_association_plot.pdf", te_association_plot, width = 4, height = 4, units = "in")

## How many overlapping TEs are there in each genome?
test <- te_overlap %>%
  filter(te_type %in% names(transposable_elements)) %>%
  group_by(genome) %>%
  summarise("Total TE overlaps (# NLR)" = paste(n(), " (", sum(gene %in% resistify$Sequence), ")", sep = ""))

test <- te_content %>%
  group_by(genome) %>%
  summarise("% TE" = sum(percentage)) %>%
  left_join(test, by = "genome")

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

# how many orthogroups are classified as N or NL?
resistify %>%
  left_join(orthogroups, by = c("Sequence" = "gene")) %>%
  group_by(Classification) %>%
  summarise(count = n_distinct(OG))
# N: 183, NL: 188, total: 771

# how many N and NL were identified by resistify?
resistify %>%
  filter(Classification %in% c("N", "NL")) %>%
  summarise(count = n_distinct(Sequence))

# how many in total were identified?
resistify %>%
  summarise(count = n_distinct(Sequence))



orthogroup_classification <- resistify %>%
  left_join(orthogroups, by = c("Sequence" = "gene")) %>%
  group_by(Classification) %>%
  summarise(count = n_distinct(OG)) %>%
  ggplot(aes(x = Classification, y = count, fill = Classification)) +
  geom_bar(stat = "identity") +
    theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Classification", y = "Number of orthogroups") +
  scale_fill_manual(values = paired_colours)

ggsave("results/orthogroups.pdf", plot_grid(orthogroup_histogram_plot, orthogroup_classification, ncol = 2, rel_widths = c(1, 2)), width = 6, height = 3, units = "in")

filtered_bed <- helixer_bed %>%
  filter(V4 %in% resistify$Sequence)

calculate_cluster_counts <- function(genome_name, maxgap_value) {
  # Subset the data for the current genome
  subsetted <- filtered_bed %>% filter(genome == genome_name)

  # Create a GRanges object
  subsetted_granges <- GRanges(seqnames = subsetted$V1, ranges = IRanges(start = subsetted$V2, end = subsetted$V3), strand = subsetted$V6)

  # Find overlaps within the specified maxgap
  overlaps <- findOverlaps(subsetted_granges, subsetted_granges, maxgap = maxgap_value, ignore.strand = TRUE)

  # Extract clusters from the overlaps
  cluster <- split(queryHits(overlaps), subjectHits(overlaps))

  # Count the clusters
  clusters_gt2 <- sum(lengths(cluster) > 2)
  clusters_eq2 <- sum(lengths(cluster) == 2)
  clusters_eq1 <- sum(lengths(cluster) == 1)

  # Create a data frame with the results
  result <- data.frame(genome = genome_name, maxgap = maxgap_value, cluster_status = c("cluster", "pair", "singleton"), count = c(clusters_gt2, clusters_eq2, clusters_eq1), stringsAsFactors = FALSE)

  return(result)
}

result_df <- data.frame(genome = character(), maxgap = numeric(), cluster_status = character(), count = numeric(), stringsAsFactors = FALSE)

for (genome_name in names(genomes)) {
  for (maxgap_value in seq(0, 200000, by = 10000)) {
    # Calculate cluster counts for the current genome and maxgap
    counts <- calculate_cluster_counts(genome_name, maxgap_value)

    # Append the results to the main data frame
    result_df <- bind_rows(result_df, counts)
  }
}

ggplot(result_df, aes(x = maxgap, y = count, fill = cluster_status)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(~ genome)

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

ggsave("results/plot.df", merged_plot, width = 6, height = 6, units = "in")

## Overlap analysis

overlapping_plot <- te_overlap %>%
  filter(gene %in% resistify$Sequence) %>%
  filter(te_type != "repeat_region") %>%
  group_by(genome, te_type) %>%
  summarise(count = n()) %>%
  ggplot(aes(y = genome, x = count, fill = te_type)) +
  geom_bar(stat = "identity") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black")) +
  labs(y = "Genome", x = "Number of overlapping TEs") +
  scale_fill_manual(values = nice_colours)

helitron_resistify <- resistify %>%
  filter(Sequence %in% te_overlap$gene)

ggsave("results/overlap.pdf", overlapping_plot, width = 6, height = 3, units = "in")

overlapping_helitron_motifs <- te_overlap %>%
  filter(gene %in% resistify$Sequence) %>%
  filter(te_type == "helitron") %>%
  select(te_id, genome) %>%
  distinct() %>%
  left_join(te_terminal_motifs, by = join_by(genome, te_id))

helitron_start <- ggseqlogo(overlapping_helitron_motifs$V2)
helitron_end <- ggseqlogo(overlapping_helitron_motifs$V3)

ggsave("motifs.png", plot_grid(helitron_start, helitron_end, nrow = 2), width = 4, height = 2, units = "in")

## known NLR analysis

refplantnlr_genomes[match("Bs4", names(refplantnlr_genomes))]

## Filter to only the best hit for each genome
filtered_refplantnlr <- refplantnlr %>%
  filter(V1 %in% names(refplantnlr_genomes)) %>%
  filter(genome == refplantnlr_genomes[V1]) %>%
  group_by(V1) %>%
  filter(V12 == max(V12)) %>%
  select(V1, V2, V3, V12, genome) %>%
  rename(refplantnlr_gene = V1, gene = V2, percent_identity = V3, bitscore = V12, genome = genome)

filtered_refplantnlr <- filtered_refplantnlr %>%
  left_join(orthogroups, by = join_by(gene)) %>%
  select(-HOG, -Gene.Tree.Parent.Clade, -genome.y) %>%
  rename(genome = genome.x)

test <- orthogroups %>%
  # get the orthogroups we're interested in
  filter(OG %in% filtered_refplantnlr$OG) %>%
  # for each orthogroup, count the number of members in each genome
  group_by(OG, genome) %>%
  summarise(count = n()) %>%
  select(OG, genome, count) %>%
  pivot_wider(names_from = genome, values_from = count, values_fill = 0) %>%
  right_join(filtered_refplantnlr, by = join_by(OG))

write.table(test, "results/refplantnlr_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

test %>%
  pivot_longer(!c(OG, refplantnlr_gene, gene, percent_identity, bitscore, genome), names_to = "genome", values_to = "count")
