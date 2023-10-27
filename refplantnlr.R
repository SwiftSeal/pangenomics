library(tidyverse)

refplantnlr <- read.table("3936022/refplantnlr.csv", header = TRUE)
resistify <- read.table("3936022/output/results.tsv", header = TRUE)

# join by ID and Sequence column
result <- refplantnlr %>%
    left_join(resistify, by = join_by("ID" == "Sequence"))

# split Subclass by ; and set as the second value
result <- result %>%
    mutate(Subclass = str_split(Subclass, ";")) %>%
    mutate(Subclass = map_chr(Subclass, ~ .x[2]))

# if Subclass is GT, set it as TNL. If GR then RNL. If G10 then CNL
result <- result %>%
    mutate(Subclass = case_when(
        Subclass == "GT" ~ "TNL",
        Subclass == "GR" ~ "RNL",
        Subclass == "G10" ~ "CNL",
        TRUE ~ Subclass
    ))

# identify rows where Subclass is not the same as Classification
result %>%
    filter(Subclass != Classification) %>%
    select(ID, Subclass, Classification) %>%
    write.table("3936022/output/incorrect_classifications.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
