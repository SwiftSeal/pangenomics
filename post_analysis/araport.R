library(ggplot2)
library(dplyr)

results <- read.table("results/results.tsv", header=TRUE, sep="\t")

head(results)

ggplot(results, aes(x = Classification)) +
  geom_bar()

summary_table <- results %>%
  group_by(Classification) %>%
  # calculate percentage which has TRUE in MADA column
  summarise(count = n(), 
            percentage_mada = 100 * sum(MADA == "True") / count,
            percentage_CJID = 100 * sum(CJID == "True") / count,
            percentage_functional = 100 * sum(Functionality == "Likely") / count)

write.table(summary_table, "results/summary_table.tsv", sep="\t", quote=FALSE, row.names=FALSE)




