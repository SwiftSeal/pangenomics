library(tidyverse)
library(ggthemes)

motif_translation = c(
  "extEDVID" = "CC",
  "bA" = "TIR",
  "aA" = "TIR",
  "bC" = "TIR",
  "aC" = "TIR",
  "bDaD1" = "TIR",
  "aD3" = "TIR",
  "VG" = "NB-ARC",
  "P-loop" = "NB-ARC",
  "RNSB-A" = "NB-ARC",
  "Walker-B" = "NB-ARC",
  "RNSB-B" = "NB-ARC",
  "RNSB-C" = "NB-ARC",
  "RNSB-D" = "NB-ARC",
  "GLPL" = "NB-ARC",
  "MHD" = "NB-ARC",
  "LxxLxL" = "LRR"
)

results <- read_tsv("refplantnlr_benchmark/resistify_output/results.tsv")
coconat <- read_tsv("refplantnlr_benchmark/resistify_output/coconat.tsv", col_names = FALSE)
annotation <- read_tsv("refplantnlr_benchmark/resistify_output/annotations.tsv")
motif <- read_tsv("refplantnlr_benchmark/resistify_output/motifs.tsv")
domain <- read_tsv("refplantnlr_benchmark/resistify_output/domains.tsv")
refplantnlr <- readxl::read_xlsx("refplantnlr_benchmark/journal.pbio.3001124.s008.xlsx")

View(full_join(results, refplantnlr, by = join_by("Sequence" == "ID")))

SEQUENCE = "Vat"

ggplot() +
  geom_rect(data = domain |> filter(Sequence == SEQUENCE), aes(xmin = Start, xmax = End, ymin = 0, ymax = 1, fill = Domain), alpha = 0.5) +
  geom_segment(data = annotation |> filter(Sequence == SEQUENCE), aes(x = Start, xend = End, y = 0, yend = 0, colour = Domain)) +
  geom_col(data = motif |> filter(Sequence == SEQUENCE), aes(x = Position, y = Probability, fill = motif_translation[Motif])) +
  geom_line(data = coconat |> filter(X1 == SEQUENCE), aes(x = X2, y = (1 - X3), colour = "CC")) +
  scale_fill_tableau() +
  scale_colour_tableau()

            