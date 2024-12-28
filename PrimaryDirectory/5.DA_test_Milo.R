library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(S4Vectors)
library(ggplot2)
set.seed(123)


# Load data
#______________
combined_data=read.csv("FinalProcessed_CD4+_CD8+/Final_Processed_Data_sub.csv"
                       ,row.names = 1,header = TRUE)

# Milo of CD8+ T cells
####################################################################################
# Load data
week3=combined_data %>%
  filter(Group %in% c("3-week CD8+ T cells" ))
rownames(week3)=paste0("W3_",1:nrow(week3))


week0=combined_data %>%
  filter(Group %in% c("Baseline CD8+ T cells" ))
rownames(week0)=paste0("W0_",1:nrow(week0))



# Combine data
combined_data <- rbind(week0, week3)
metabolic_space <- c("Glut1", "MCT1",  "G6PD", "IDH2", "ATP5A", "CD98", "CD36", "CPT1A","GAPDH", "ACAC")
combined_data_metaolic_space <- combined_data[,metabolic_space]
combined_data_metaolic_space_t <- t(combined_data_metaolic_space)


# Categorize the cells into Exhausted and Functional
labels_week0 = ifelse(week0$Exhausted, "Exhausted", "Functional")
labels_week3 <- ifelse(week3$Exhausted, "Exhausted", "Functional")

combined_labels <- c(labels_week0, labels_week3)

# ANOTHER HACK: Milo needs replicates, we simply produce those synthetically and make up random values for 3 replicates for each condition
#create synthetic replicates for statistical testing
x <- sample(1:3, nrow(week0) , replace=T) # For baseline cells
y <- sample(4:6, nrow(week3), replace=T) # For 3-week cells


# Create a SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(intensities = as.matrix(combined_data_metaolic_space_t)))
colData_sce <- DataFrame(row.names=colnames(combined_data_metaolic_space_t),cell_id = colnames(combined_data_metaolic_space_t), Group = c(rep("Baseline",nrow(week0)),rep("Induced cells",nrow(week3))), replicate_id = c(x,y),State=combined_labels)
colData(sce) = colData_sce

#WE USE UMAP ONLY FOR VISUALIZATION: FOR THE SPACE IN WHICH WE APPLY MILO WE WANT TO USE ETHER PCA(and use all dimensons)
# or we use a small hack here and actually input the origional data, acting like it is some reduced dimension...
reducedDim(sce, "TRUE_VALUES") <- as.matrix(combined_data_metaolic_space)
head(reducedDim(sce, "TRUE_VALUES"))


# Add UMAP data: ONLY FOR VISUALIZATION, look at parameter d in functions below, we want to use ALL DIMENSIONS there and not only the 2 UMAP dimensions
sce <- runUMAP(sce,exprs_values = "intensities",name = "UMAP")


# Create a Milo object
#______________
sce <- Milo(sce)
# Construct KNN graph
sce <- buildGraph(sce, k=30, d=length(metabolic_space), reduced.dim = "TRUE_VALUES")
# Defining representative neighbourhoods on the KNN graph
sce <- makeNhoods(sce, prop = 1, k=30, d=length(metabolic_space), refined = TRUE, reduced_dims = "TRUE_VALUES")
plotNhoodSizeHist(sce)

# Counting cells in neighbourhoods
sce <- countCells(sce, meta.data = data.frame(colData(sce)), sample="replicate_id")
head(nhoodCounts(sce))

# Defining experimental design(
sce_design <- data.frame(colData(sce))[,c("Group","replicate_id"),drop=FALSE]
# Convert batch info from integer to factor
sce_design$Group <- as.factor(sce_design$Group)
sce_design <- distinct(sce_design)
rownames(sce_design) <- sce_design$replicate_id
sce_design


#COMPUTING MILO RESULTS
#______________

# Computing neighbourhood connectivity
sce <- calcNhoodDistance(sce, d=10, reduced.dim = "TRUE_VALUES")
# Testing
da_results <- testNhoods(sce, design = ~ Group, design.df = sce_design, reduced.dim = "TRUE_VALUES")

# Check the distribution of p-value
ggplot(da_results, aes(PValue))+ geom_histogram(bins=50)


# Visualize the FDR results with a volcano plot
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point(alpha=0.7) +
  geom_hline(yintercept = 1)

# Visualize the p-value results with a volcano plot
ggplot(da_results, aes(logFC, -log10(PValue))) +
  geom_point(alpha=7) +
  geom_hline(yintercept = 1) 



#PLOT RESULTS
#______________
sce <- buildNhoodGraph(sce)
## Plot single-cell UMAP
umap_pl <- plotReducedDim(sce, dimred = "UMAP", colour_by="Group",
                          text_size = 3, point_size=0.5) +
  guides(fill="none")
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(sce, da_results, layout="UMAP",alpha=0.1) 
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")


# Assign a cell type label to each neighbourhood by finding the most abundant cell type within cells in each neighbourhood
da_results <- annotateNhoods(sce, da_results, coldata_col = "State")
head(da_results)

# Visualize the distribution of DA Fold Changes in different cell types
ggplot(da_results, aes(State_fraction)) + geom_histogram(bins=50)

# Exclude neighbourhoods that are a mix of cell types
da_results$celltype <- ifelse(da_results$State_fraction < 0.8, "Mixed", da_results$State)

# Visualize the distribution of DA Fold Changes in different cell types
plotDAbeeswarm(da_results, group.by = "State")


# Identifying signatures of DA subpopulations
dge_smp <- findNhoodMarkers(sce, da_results,
                            assay = "intensities", gene.offset = FALSE, da.fdr = 0.1,
                            aggregate.samples = FALSE, 
                            subset.nhoods = da_results$State %in% c("Exhausted","Functional")
)

head(dge_smp)


markers <- dge_smp[which(dge_smp$adj.P.Val_1 < 0.1 ), "GeneID"]


plotNhoodExpressionDA(sce, da_results, features = markers,
                      subset.nhoods = da_results$State %in% c("Exhausted","Functional"),
                      assay="intensities",
                      scale_to_1 = TRUE, cluster_features = TRUE
)
# Milo of CD4+ T cells
####################################################################################

# Load data
combined_data=read.csv("FinalProcessed_CD4+_CD8+/Final_Processed_Data_sub.csv"
                       ,row.names = 1,header = TRUE)

week3=combined_data %>%
  filter(Group %in% c("3-week CD4+ T cells" ))
rownames(week3)=paste0("W3_",1:nrow(week3))


week0=combined_data %>%
  filter(Group %in% c("Baseline CD4+ T cells" ))
rownames(week0)=paste0("W0_",1:nrow(week0))



# Combine data
combined_data <- rbind(week0, week3)
metabolic_space <- c("Glut1", "MCT1",  "G6PD", "IDH2", "ATP5A", "CD98", "CD36", "CPT1A","GAPDH", "ACAC")
combined_data_metaolic_space <- combined_data[,metabolic_space]
combined_data_metaolic_space_t <- t(combined_data_metaolic_space)


# Categorize the cells into Exhausted and Functional
labels_week0 = ifelse(week0$Exhausted, "Exhausted", "Functional")
labels_week3 <- ifelse(week3$Exhausted, "Exhausted", "Functional")

combined_labels <- c(labels_week0, labels_week3)

# ANOTHER HACK: Milo needs replicates, we simply produce those synthetically and make up random values for 3 replicates for each condition
#create synthetic replicates for statistical testing
x <- sample(1:3, nrow(week0) , replace=T) # For baseline cells
y <- sample(4:6, nrow(week3), replace=T) # For 3-week cells


# Create a SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(intensities = as.matrix(combined_data_metaolic_space_t)))
colData_sce <- DataFrame(row.names=colnames(combined_data_metaolic_space_t),cell_id = colnames(combined_data_metaolic_space_t), Group = c(rep("Baseline",nrow(week0)),rep("Induced cells",nrow(week3))), replicate_id = c(x,y),State=combined_labels)
colData(sce) = colData_sce

#WE USE UMAP ONLY FOR VISUALIZATION: FOR THE SPACE IN WHICH WE APPLY MILO WE WANT TO USE ETHER PCA(and use all dimensons)
# or we use a small hack here and actually input the origional data, acting like it is some reduced dimension...
reducedDim(sce, "TRUE_VALUES") <- as.matrix(combined_data_metaolic_space)
head(reducedDim(sce, "TRUE_VALUES"))


# Add UMAP data: ONLY FOR VISUALIZATION, look at parameter d in functions below, we want to use ALL DIMENSIONS there and not only the 2 UMAP dimensions
sce <- runUMAP(sce,exprs_values = "intensities",name = "UMAP")


# Create a Milo object
#______________
sce <- Milo(sce)
# Construct KNN graph
sce <- buildGraph(sce, k=30, d=length(metabolic_space), reduced.dim = "TRUE_VALUES")
# Defining representative neighbourhoods on the KNN graph
sce <- makeNhoods(sce, prop = 1, k=30, d=length(metabolic_space), refined = TRUE, reduced_dims = "TRUE_VALUES")
plotNhoodSizeHist(sce)

# Counting cells in neighbourhoods
sce <- countCells(sce, meta.data = data.frame(colData(sce)), sample="replicate_id")
head(nhoodCounts(sce))

# Defining experimental design(
sce_design <- data.frame(colData(sce))[,c("Group","replicate_id"),drop=FALSE]
# Convert batch info from integer to factor
sce_design$Group <- as.factor(sce_design$Group)
sce_design <- distinct(sce_design)
rownames(sce_design) <- sce_design$replicate_id
sce_design


#COMPUTING MILO RESULTS
#______________

# Computing neighbourhood connectivity
sce <- calcNhoodDistance(sce, d=10, reduced.dim = "TRUE_VALUES")
# Testing
da_results <- testNhoods(sce, design = ~ Group, design.df = sce_design, reduced.dim = "TRUE_VALUES")

# Check the distribution of p-value
ggplot(da_results, aes(PValue))+ geom_histogram(bins=50)


# Visualize the FDR results with a volcano plot
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point(alpha=0.7) +
  geom_hline(yintercept = 1)

# Visualize the p-value results with a volcano plot
ggplot(da_results, aes(logFC, -log10(PValue))) +
  geom_point(alpha=7) +
  geom_hline(yintercept = 1) 



#PLOT RESULTS
#______________
sce <- buildNhoodGraph(sce)
## Plot single-cell UMAP
umap_pl <- plotReducedDim(sce, dimred = "UMAP", colour_by="Group",
                          text_size = 3, point_size=0.5) +
  guides(fill="none")
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(sce, da_results, layout="UMAP",alpha=0.1) 
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")


# Assign a cell type label to each neighbourhood by finding the most abundant cell type within cells in each neighbourhood
da_results <- annotateNhoods(sce, da_results, coldata_col = "State")
head(da_results)

# Visualize the distribution of DA Fold Changes in different cell types
ggplot(da_results, aes(State_fraction)) + geom_histogram(bins=50)

# Exclude neighbourhoods that are a mix of cell types
da_results$celltype <- ifelse(da_results$State_fraction < 0.8, "Mixed", da_results$State)

# Visualize the distribution of DA Fold Changes in different cell types
plotDAbeeswarm(da_results, group.by = "State")


# Identifying signatures of DA subpopulations
dge_smp <- findNhoodMarkers(sce, da_results,
                            assay = "intensities", gene.offset = FALSE, da.fdr = 0.1,
                            aggregate.samples = FALSE, 
                            subset.nhoods = da_results$State %in% c("Exhausted","Functional")
)

head(dge_smp)


markers <- dge_smp[which(dge_smp$adj.P.Val_1 < 0.1 ), "GeneID"]


plotNhoodExpressionDA(sce, da_results, features = markers,
                      subset.nhoods = da_results$State %in% c("Exhausted","Functional"),
                      assay="intensities",
                      scale_to_1 = TRUE, cluster_features = TRUE
)
