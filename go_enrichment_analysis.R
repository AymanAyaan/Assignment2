#We read the CSV file
deg_data <- read.csv("deg_results_GSE9750.csv", stringsAsFactors = FALSE)
cat("Dataset Columns:\n")
print(colnames(deg_data))
# 1. Set the column name for your gene symbols
gene_col <- "X"
#cleaned Data
cleaned_data <- deg_data[!is.na(deg_data[[gene_col]]) & deg_data[[gene_col]] != "" & deg_data[[gene_col]] != "---", ]
# Unique gene
unique_genes <- unique(cleaned_data[[gene_col]])

write.table(unique_genes, 
            file = "GSE9750_string_input.txt", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

node_table <- read.csv("gse9750_cytoscape_node_table.csv", stringsAsFactors = FALSE)


# Load the verified biological libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)


module1_genes <- c("NUSAP1", "CDC20", "CCNB1", "NEK2", "PRC1", "UBE2C", "MELK", 
                   "RRM2", "KIF4A", "CENPF", "CEP55", "TOP2A", "DLGAP5", "KIF20A", "CDK1")

module2_genes <- c("RFC4", "GMNN", "MCM6", "MCM2", "MCM5", "CDC7", "GINS1")

run_go_enrichment <- function(gene_list, module_name) {
  message("Executing pathway analysis for: ", module_name)
  
  entrez_ids <- bitr(gene_list, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
  
  go_results <- enrichGO(gene          = entrez_ids$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "BP",           
                         pAdjustMethod = "BH",           
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)  
  
  return(go_results)
}

m1_analysis <- run_go_enrichment(module1_genes, "Module1")
m2_analysis <- run_go_enrichment(module2_genes, "Module2")

write.csv(as.data.frame(m1_analysis), file = "GO_enrichment_Module1_results.csv", row.names = FALSE)
write.csv(as.data.frame(m2_analysis), file = "GO_enrichment_Module2_results.csv", row.names = FALSE)


