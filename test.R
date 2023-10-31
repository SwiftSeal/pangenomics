library("data.table")
library("ape")
tree <- read.tree("results/orthofinder/Species_Tree/SpeciesTree_rooted.txt")

if(is.ultrametric(tree)) {
    utree <- tree
} else{
    utree <- chronos(tree)
}

write.tree(utree, "results/cafe/SpeciesTree_rooted_ultrametric.txt")

hog <- fread("results/orthofinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")
hog[, OG := NULL]
hog[, `Gene Tree Parent Clade` := NULL]
hog <- melt(hog, id.vars='HOG', variable.name='species', value.name='pid')
hog <- hog[pid != '']
hog$n <- sapply(hog$pid, function(x) length(strsplit(x, ', ')[[1]]))

# Exclude HOGs with lots of genes in a one or more species. 
# See also cafe tutorial about filtering gene families
keep <- hog[, list(n_max=max(n)), HOG][n_max < 100]$HOG
hog <- hog[HOG %in% keep]

# Exclude HOGs present in only 1 species
keep <- hog[, .N, HOG][N > 1]$HOG
hog <- hog[HOG %in% keep]

counts <- dcast(hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')

fwrite(counts, 'results/cafe/hog_gene_counts.tsv', sep='\t')