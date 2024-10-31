library(tidyverse)
library(ggtree)
library(ggthemes)
library(cowplot)
library(ggtext)

# Configuration

genomes <- c(
  "bulbocastanum" = "italic(S. bulbocastanum)",
  "annuum" = "italic(C. annuum)",
  "candolleanum" = "italic(S. candolleanum)",
  "etuberosum" = "italic(S. etuberosum)",
  "lycopersicum" = "italic(S. lycopersicum)",
  "tuberosum_phureja_E4_63" = "italic(S. tuberosum) (Phureja E4-63)",
  "tuberosum_tuberosum_RH10_15" = "italic(S. tuberosum) (tuberosum RH10-15)",
  "chilense" = "italic(S. chilense)",
  "galapagense" = "italic(S. galapagense)",
  "neorickii" = "italic(S. neorickii)",
  "tuberosum_phureja_E86_69" = "italic(S. tuberosum) (Phureja E86-69)",
  "tuberosum_tuberosum_RH" = "italic(S. tuberosum) (tuberosum RH)",
  "chmielewskii" = "italic(S. chmielewskii)",
  "habrochaites" = "italic(S. habrochaites)",
  "peruvianum" = "italic(S. peruvianum)",
  "tuberosum_stenotomum_A6_26" = "italic(S. tuberosum) (Stenotomum A6-26)",
  "verrucosum" = "italic(S. verrucosum)",
  "corneliomulleri" = "italic(S. corneliomulleri)",
  "lycopersicoides" = "italic(S. lycopersicoides)",
  "pimpinellifolium" = "italic(S. pimpinellifolium)",
  "tuberosum_stenotomum_PG6359" = "italic(S. tuberosum) (Stenotomum PG6359)"
)

genomes <- c(
  "bulbocastanum" = "*S. bulbocastanum*",
  "annuum" = "*C. annuum*",
  "candolleanum" = "*S. candolleanum*",
  "etuberosum" = "*S. etuberosum*",
  "lycopersicum" = "*S. lycopersicum*",
  "tuberosum_phureja_E4_63" = "*S. tuberosum* (Phureja E4-63)",
  "tuberosum_tuberosum_RH10_15" = "*S. tuberosum* (tuberosum RH10-15)",
  "chilense" = "*S. chilense*",
  "galapagense" = "*S. galapagense*",
  "neorickii" = "*S. neorickii*",
  "tuberosum_phureja_E86_69" = "*S. tuberosum* (Phureja E86-69)",
  "tuberosum_tuberosum_RH" = "*S. tuberosum* (tuberosum RH)",
  "chmielewskii" = "*S. chmielewskii*",
  "habrochaites" = "*S. habrochaites*",
  "peruvianum" = "*S. peruvianum*",
  "tuberosum_stenotomum_A6_26" = "*S. tuberosum* (Stenotomum A6-26)",
  "verrucosum" = "*S. verrucosum*",
  "corneliomulleri" = "*S. corneliomulleri*",
  "lycopersicoides" = "*S. lycopersicoides*",
  "pimpinellifolium" = "*S. pimpinellifolium*",
  "tuberosum_stenotomum_PG6359" = "*S. tuberosum* (Stenotomum PG6359)"
)

tuberising <- c(
    "bulbocastanum" = TRUE,
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

paired_colours = c(
  "CN" = "#a6cee3", 
  "CNL" = "#1f78b4",
  "N" = "#b2df8a",
  "NL" = "#33a02c",
  "RN" = "#fb9a99",
  "RNL" = "#e31a1c",
  "TN" = "#fdbf6f",
  "TNL" = "#ff7f00"
)

nice_theme <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank(),

  strip.text = element_blank())

# Import data

## Merge all resistify results into a single dataframe
resistify <- map(
  names(genomes),
    read.table(paste0("output/", .x, "/results.tsv"), header = TRUE, sep = "\t") |> mutate(genome = .x)
) |>
  bind_rows()

## Orthofinder tree
genome_tree <- read.tree("output/SpeciesTree_rooted.txt")

## EDTA calculated % of TE content
te_content <- map(
  names(genomes),
    read.table(
    paste0("output/", .x, ".fa.mod.EDTA.TEanno.sum"),
    skip = 6,
    fill = TRUE,
    nrows = 12,
    col.names = c("te_classification", "te_count", "te_masked_bp", "percentage")
  ) |> filter(te_count != "--") |> mutate(genome = .x, percentage = as.numeric(str_remove(percentage, "%")))
) |>
  bind_rows()

## Orthofinder orthogroups
orthogroups <- read.table("output/N0.tsv", sep = "\t", header = TRUE) |>
  pivot_longer(!1:3, names_to = "genome", values_to = "gene") |>
  separate_rows(gene, sep = ", ")

## Refplantnlr blast results

refplantnlr <- map(
  names(genomes),
    read.table(paste0("output/", .x, ".pep.blast.txt"), sep = "\t") |> mutate(genome = .x)
) |>
  bind_rows()

## genome size

genome_size <- data.frame(genome = character(), Total_Length = numeric(), stringsAsFactors = FALSE)
for (genome in names(genomes)) {
  # Read the fourth line from the file
  total_length_line <- readLines(paste0("output/", genome, ".fa.mod.EDTA.TEanno.sum"), n = 4)[4]

  # Extract the numeric value using regular expression
  total_length <- as.numeric(regmatches(total_length_line, regexpr("\\d+", total_length_line)))

  # Add the results to the data frame
  genome_size <- rbind(genome_size, data.frame(genome = genome, Total_Length = total_length, stringsAsFactors = FALSE))
}

## Overlaps

te_overlap <- map(
  names(genomes),
    read.table(paste0("output/", .x, ".intersect.gff"))|>
    filter(V3 == "gene") |>
    mutate(overlap_percentage = V19/(V5 + 1 - V4)) |>
    filter(overlap_percentage > 0.9) |>
    mutate(genome = .x) 
) |>
  bind_rows()

## TE terminal motifs

te_terminal_motifs <- map(names(genomes),   read.table(paste0("results/edta/", .x, "_intact_terminal.tsv"), sep = "\t", comment.char = "") %>% mutate(genome = .x)) %>%
  bind_rows()

te_terminal_motifs <- te_terminal_motifs %>%
  mutate(te_id = str_remove(V1, "\\|.*"))

## TE intact gff

te_intact_gff <- map(names(genomes),   read.table(paste0("results/edta/", .x, ".mod.EDTA.intact.gff3"), sep = "\t", comment.char = "") %>% mutate(genome = .x)) %>%
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

summary(lm(intact_te_count   tuberising, data = genome_summary))

summarise(group_by(resistify, genome), count = n()) |>
  mutate(tuberising = tuberising[match(genome, names(tuberising))]) |>
  ggplot(aes(x = tuberising, y = count)) +
  geom_boxplot() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black")
  )

te_association_plot <- inner_join(
  te_content,
  summarise(group_by(resistify, genome), count = n()),
  by = "genome"
) |>
  filter(te_classification != "unknown") |>
  mutate(tuberising = tuberising[match(genome, names(tuberising))]) |>
  ggplot(aes(x = percentage, y = count, colour = tuberising)) +
  geom_point() +
  facet_wrap(  te_classification, scales = "free_x") +
  labs(x = "TE content (%)", y = "Number of resistance genes")

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

## how many orthogroups are classified as N or NL?
#resistify %>%
#  left_join(orthogroups, by = c("Sequence" = "gene")) %>%
#  group_by(Classification) %>%
#  summarise(count = n_distinct(OG))
## N: 183, NL: 188, total: 771
#
## how many N and NL were identified by resistify?
#resistify %>%
#  filter(Classification %in% c("N", "NL")) %>%
#  summarise(count = n_distinct(Sequence))
#
## how many in total were identified?
#resistify %>%
#  summarise(count = n_distinct(Sequence))

orthogroup_classification <- resistify |>
  left_join(orthogroups, by = c("Sequence" = "gene")) |>
  group_by(Classification) |>
  summarise(count = n_distinct(HOG)) |>
  ggplot(aes(x = Classification, y = count, fill = Classification)) +
  geom_bar(stat = "identity") +
  labs(x = "Classification", y = "Number of orthogroups") +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = paired_colours)

# find majority classification for each OG
OG_ids <- orthogroups |>
  inner_join(select(resistify, Sequence, Classification), by = join_by(gene == Sequence)) |>
  group_by(HOG) |>
  summarise(Classification = Classification |> table() |> which.max() |> names())

orthogroup_histogram_plot <- orthogroups |>
  inner_join(select(resistify, Sequence, Classification), by = join_by(gene == Sequence)) |>
  group_by(genome, HOG) |>
  summarise(count = n()) |>
  group_by(HOG) |>
  summarise(count = n()) |>
  left_join(OG_ids) |>
  ggplot(aes(x = count, fill = Classification)) +
  geom_bar() +
  labs(x = "Number of genomes", y = "Number of orthogroups") +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = paired_colours)

## Filter to only the best hit for each genome
filtered_refplantnlr <- refplantnlr |>
  filter(V1 %in% names(refplantnlr_genomes)) |>
  filter(genome == refplantnlr_genomes[V1]) |>
  group_by(V1) |>
  filter(V12 == max(V12)) |>
  select(V1, V2, V3, V12, genome) |>
  rename(refplantnlr_gene = V1, gene = V2, percent_identity = V3, bitscore = V12, genome = genome)

filtered_refplantnlr <- filtered_refplantnlr |>
  left_join(orthogroups, by = join_by(gene)) |>
  select(-OG, -Gene.Tree.Parent.Clade, -genome.y) |>
  rename(genome = genome.x)

refplantnlr_counts <- orthogroups |>
  # get the orthogroups we're interested in
  filter(HOG %in% filtered_refplantnlr$HOG) |>
  # for each orthogroup, count the number of members in each genome
  group_by(HOG, genome) |>
  summarise(count = n()) |>
  select(HOG, genome, count) |>
  right_join(filtered_refplantnlr, by = join_by(HOG))

refplantnlr_plot <- ggplot(refplantnlr_counts, aes(x = genomes[genome.x], y = refplantnlr_gene, fill = count)) +
  geom_tile() +
  geom_text(aes(label = count)) +
  theme_bw(base_size = 8) +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_markdown(angle = 90, hjust = 0.95)
  ) +
  scale_fill_gradient(low = "white", high = "#1f78b4")

merged_plot <- ggarrange(
  ggarrange(orthogroup_classification, orthogroup_histogram_plot, ncol = 2, labels = c("A", "B")),
  refplantnlr_plot,
  nrow = 2,
  labels = c("", "C")
)

ggsave("orthogroup_plot.png", width )

# Phylogenetic tree ------------------------------------------------------------
resistify_summary <- resistify |>
    group_by(genome, Classification) |>
    summarise(count = n())

overlap_summary <- te_overlap |>
  mutate(gene = paste0(gsub("ID=", "", V9), ".1")) |>
  filter(gene %in% resistify$Sequence) |>
  filter(str_detect(V18, "Method=structural")) |>
  filter(!str_detect(V12, "repeat_region")) |>
  mutate(te_type = transposable_elements[V12]) |>
  group_by(genome, te_type) |>
  summarise(count = n())

tree_plot <- ggtree(genome_tree) +
  geom_tiplab(aes(label = genomes[label]), parse = TRUE) +
  xlim(0, 0.4)

nice_theme <- theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank()
  )

resistify_plot <- ggplot(resistify_summary, aes(y = genome, x = count, fill = Classification)) +
  geom_bar(stat = "identity") +
  labs(x = "# NLRs") +
  scale_fill_manual(name = "NLR Classification", values = paired_colours) +
  nice_theme

overlap_plot <- ggplot(overlap_summary, aes(y = genome, x = count, fill = te_type)) +
  geom_bar(stat = "identity") +
  labs(x = "# TE-embedded NLRs") +
  scale_fill_tableau(name = "TE type") +
  nice_theme

final <- resistify_plot |> insert_left(tree_plot, width = 1.5) |> insert_right(overlap_plot)
ggsave2("phylogeny.pdf", final, width = 12, height = 6)

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





write.table(test, "results/refplantnlr_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

test %>%
  pivot_longer(!c(OG, refplantnlr_gene, gene, percent_identity, bitscore, genome), names_to = "genome", values_to = "count")
