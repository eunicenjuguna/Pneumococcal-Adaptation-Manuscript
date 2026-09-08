
library(phytools)

library(ape)

# 1. Load your tree
tree <- read.tree("core.full.final_tree.tre")

# 2. Calculate cophenetic distance (genetic distance between all pairs of tips)
dist_matrix <- cophenetic.phylo(tree)

# 3. Perform hierarchical clustering on the distance matrix
# This does NOT require the tree to be ultrametric
hc <- hclust(as.dist(dist_matrix), method = "average")

# 4. Visualize the tree and the height to choose your 'h'
plot(hc, labels = FALSE)
abline(h = 170, col = "red") # Adjust this height (h) as needed

# 5. Cut the tree into clusters
# Adjust 'h' until you get the number of sub-lineages you expect
clusters <- cutree(hc, h = 170)

# 6. Save to data frame
sublineage_df <- data.frame(
  sample_id = names(clusters), 
  sub_cluster = clusters
)

# Preview
head(sublineage_df)

# Check the number of clusters for different heights
for(h_val in c(150,  200)) {
  num_clusters <- length(unique(cutree(hc, h = h_val)))
  message(paste("At height", h_val, "you get", num_clusters, "clusters"))
}

library(dplyr)
library(readxl)

# 1. Load your metadata
# Ensure the ID column in your excel file matches the ID column in sublineage_df
metadata <- read.csv("GPSC22-metadata.csv")

# 2. Link clusters to metadata
# Replace 'sample_id' and 'ID_in_metadata' with the actual column names in your files
final_metadata <- metadata %>%
  left_join(sublineage_df, by = c("lanes" = "sample_id"))

# 3. Verify the merge
# Check if any samples failed to match (where sub_cluster is NA)
missing_matches <- final_metadata %>% filter(is.na(sub_cluster))

if(nrow(missing_matches) > 0) {
  warning("Some samples did not match! Check your ID naming consistency.")
} else {
  message("Success! All samples linked to clusters.")
}

# 4. Save the combined file
write.csv(final_metadata, "metadata_with_sublineages.csv", row.names = FALSE)

