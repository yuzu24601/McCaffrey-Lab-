##Friday, December 20th- making a heatmap showing metabolic pathway expression 
###Adapting "5.0 trying ggalign" code to make heatmap
##Pathway scores generated using KEGG, then added them as columns to "data_merged_all"

library(ggplot2)
library(tibble)

generate_scaled_ggheatmap_for_pathways <- function(data_merged_all, pathway_scores, cluster_numbers) {
  # Fetch the pathway scores and cluster information
  expression_data_all_clusters <- FetchData(data_merged_all, vars = pathway_scores)
  
  # Get the cluster IDs of all cells in the Seurat object
  cluster_ids <- Idents(data_merged_all)
  
  # Add cluster information to the pathway score data
  expression_data_all_clusters$cluster <- cluster_ids
  
  # Aggregate pathway scores by cluster (mean score per cluster)
  expression_data_cluster_avg <- expression_data_all_clusters %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    column_to_rownames(var = "cluster")
  
  # Scale the columns (pathways) to have mean = 0 and standard deviation = 1
  scaled_expression_data <- scale(expression_data_cluster_avg)
  
  # Generate the heatmap using ggalign::ggheatmap
  ggheatmap(scaled_expression_data) + 
    ggtitle("Pathway Scores Across Clusters") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Define the pathway scores
pathway_scores <- c("Glycolysis_score1", 
                    "Pentose_Phosphate_score1", 
                    "Oxidative_Phosphorylation_score1", 
                    "Glutamate_Glutamine_score1", 
                    "Hypoxia_score1")

# Example usage for all clusters
print(generate_scaled_ggheatmap_for_pathways(data_merged_all, pathway_scores, cluster_numbers = 0:12))





