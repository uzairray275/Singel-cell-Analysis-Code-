Integrated Single Cell Analysis and Machine Learning for Precision Diagnosis in Human Papillomavirus Associated Cancers 
##############################################################
Code for Single Cell Analysis
##############################################################

STEP1: Install R along with R studio and run the following code in environment

#!/usr/bin/env Rscript

# R script to install requirements for exercises -------------------------------
dir.exists()
dir_exists <- sapply(datadirs, function(dir) dir.exists(dir))
print(dir_exists)
## a vector of packages to install (edit in this section) ----------------------
### packages could be either on CRAN or bioconductor
pkgs <- c("ggplot2", "BiocManager", "sctransform",
          "devtools", "cowplot", "matrixStats",
          "ggbeeswarm", "ggnewscale", "msigdbr",
          "Seurat", "bit64", "scater",
          "AnnotationDbi",
          "SingleR", "clusterProfiler", "celldex",
          "dittoSeq", "DelayedArray",
          "DelayedMatrixStats",
          "limma", "SingleCellExperiment",
          "SummarizedExperiment",
          "slingshot", "batchelor",
          "clustree", "edgeR")
for (pkg in pkgs) {
  library(pkg, character.only = TRUE)
}
install.packages("celldex")
BiocManager::install("celldex")
library(celldex)
install.packages(BiocManager)
setwd("C:/Users/Uzair/OneDrive - National University of Sciences & Technology/Desktop/uzair thesis material/CERVICAL CANCER")
sampleinfo <- read.csv("sample_info.csv")
datadirs <- file.path(".", sampleinfo$ID)
names(datadirs) <- gsub("_", "-", sampleinfo$ID)
datadirs
library(Seurat)
sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
seu <- Seurat::CreateSeuratObject(counts = sparse_matrix,
                                  project = "Primary",
                                  min.cells = 3,
                                  min.features = 100)
seu
View(seu)
seu@meta.data
metadata_table <- seu@meta.data
write.csv(metadata_table, "Meta_Data_Table.csv", row.names=TRUE)
table(seu@active.ident)
Cells <- table(seu@active.ident)
write.csv(Cells, "Cells.csv", row.names=TRUE)
head(seu@meta.data)
summary(seu@meta.data$nCount_RNA)
summary(seu@meta.data$nFeature_RNA)
mypalette <- c("#FF0000", "#00796B", "#27CED7", "#FF00CC", "#00FF00","purple","pink","orange","yellow")
mypalette2 <- c("#27CED7", "#FF00CC", "#00FF00","#FF2A00", "#00468B","purple","pink","orange","yellow")
mypalette3 <- c("#27CED7", "#FF00CC", "#00FF00","#FF2A00", "#00468B","purple","pink","orange","yellow","#FF0000", "#00796B","brown")
hist(seu$nCount_RNA, col = mypalette, main = paste0("Histogram of RNA counts per cell"))
hist(seu$nFeature_RNA, col = mypalette, main = paste0("Histogram of Gene counts per cell"))
Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", col = mypalette2)
Seurat::VlnPlot(seu, features = "nCount_RNA", cols = mypalette2)
Seurat::VlnPlot(seu, features = "nFeature_RNA", cols = mypalette)
Seurat::VlnPlot(seu, features = c("nCount_RNA",
                                  "nFeature_RNA"))
# mitochondrial genes
seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^MT-",
                                    col.name = "percent.mito")
# ribosomal genes
seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^RP[SL]",
                                    col.name = "percent.ribo")
# hemoglobin genes (but not HBP)
seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^HB[^(P)]",
                                    col.name = "percent.globin")
head(seu@meta.data)
Seurat::VlnPlot(seu, features = "percent.mito", cols = mypalette2)
Seurat::VlnPlot(seu, features = "percent.mito", cols = mypalette2, pt.size = 0)
Seurat::VlnPlot(seu, features = "percent.ribo", cols = mypalette2)
Seurat::VlnPlot(seu, features = "percent.ribo", cols = mypalette2, pt.size = 0)
Seurat::VlnPlot(seu, features = "percent.globin", cols = mypalette2)
Seurat::VlnPlot(seu, features = "percent.globin", cols = mypalette2, pt.size = 0)
Seurat::VlnPlot(seu, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))
Seurat::FeatureScatter(seu,
                      feature1 = "percent.mito",
                       feature2 = "percent.ribo", cols = mypalette2)
library(ggplot2)
library(Matrix)
library(Seurat)
most_expressed_boxplot <- function(object, ngenes = 20){
   # matrix of raw counts
  cts <- Seurat::GetAssayData(object, assay = "RNA", slot = "counts")
    # get percentage/cell
  cts <- t(cts)/colSums(cts)*100
  medians <- apply(cts, 2, median)
   # get top n genes
  most_expressed <- order(medians, decreasing = T)[ngenes:1]
  most_exp_matrix <- as.matrix((cts[,most_expressed]))
   # prepare for plotting
  most_exp_df <- stack(as.data.frame(most_exp_matrix))
  colnames(most_exp_df) <- c("perc_total", "gene")
    # boxplot with ggplot2
  boxplot <- ggplot(most_exp_df, aes(x=gene, y=perc_total)) +
    geom_boxplot() +
    coord_flip()
  return(boxplot)
}
most_expressed_boxplot(seu, 20)

seu <- subset(seu, subset = nFeature_RNA > 200 &
                nFeature_RNA < 5000 &
                percent.mito < 8)

Seurat::VlnPlot(seu, features = c("nFeature_RNA",
                                  "percent.mito"))
Seurat::GetAssayData(seu)[1:30,1:30] 
Assay_Data_Before_Normalization <- Seurat::GetAssayData(seu)[1:30,1:30]
write.csv(Assay_Data_Before_Normalization, "Assay_Data_Before_Normalization.csv", row.names=TRUE)

############Normalization##############

seu <- Seurat::NormalizeData(seu,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)
Seurat::GetAssayData(seu)[1:30,1:30]
Assay_Data_After_Normalization <- Seurat::GetAssayData(seu)[1:30,1:30]
write.csv(Assay_Data_After_Normalization, "Assay_Data_After_Normalization.csv", row.names=TRUE)

#############variable feature############

seu <- Seurat::FindVariableFeatures(seu,
                                    selection.method = "vst",
                                    nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(seu), 10)
top10
vf_plot <- Seurat::VariableFeaturePlot(seu)
Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE)

############Scaling###########
seu <- Seurat::ScaleData(seu,
                         features = rownames(seu))
plots <- VlnPlot(seu, features = c("percent.mito",
                                   "percent.ribo",
                                   "percent.globin"),
                 pt.size =0,
                 combine = FALSE, cols = mypalette2)
for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot() + theme(legend.position = 'none')
}
CombinePlots(plots)

##################PCA############

seu <- Seurat::RunPCA(seu)
Seurat::DimPlot(seu, reduction = "pca")
Seurat::DimPlot(seu, reduction = "pca", cols = mypalette2, pt.size = 1)

##############Heatmap############

Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)
Seurat::ElbowPlot(seu, ndims = 40)
seu <- Seurat::RunUMAP(seu, dims = 1:25)
Seurat::DimPlot(seu, reduction = "umap", cols = mypalette2)
# The default number of neighbours is 30. If your dataset is small, a decrease in the number of neighbors can be considered
seu <- Seurat::RunUMAP(seu, dims = 1:25, n.neighbors = 30) 
Seurat::DimPlot(seu, reduction = "umap", cols = mypalette2)
# Taking too few PCs we see everything looks connected
seu <- Seurat::RunUMAP(seu, dims = 1:5)      
Seurat::DimPlot(seu, reduction = "umap", cols = mypalette2)
#if more precision makes sense, for instance, if the genes that is of interest for your study is not present when the RunPCA was calculated, then an increase in the number of components calculated at start might be interesting to be considered
seu <- Seurat::RunUMAP(seu, dims = 1:50)     
Seurat::DimPlot(seu, reduction = "umap", cols = mypalette2)

#changing back to 30 clusters

seu <- Seurat::RunUMAP(seu, dims = 1:30)   
Seurat::DimPlot(seu, reduction = "umap", cols = mypalette2)
seu <- Seurat::RunTSNE(seu, dims = 1:30)  
Seurat::DimPlot(seu, reduction = "tsne", group.by = 'orig.ident', cols = mypalette2)
BiocManager::install("harmony")
library(Rcpp) 
seu <- seu %>% RunHarmony('orig.ident', plot_convergence = F)

###############Integration#################


library(ggplot2)
library(tidyverse)
seu@reductions
seu.embed <- Embeddings(seu, "harmony")
seu.embed

############Clustering############

library(ggraph)
library(clustree)
seu_clusters_UMAP <- seu %>%
  RunUMAP(reduction = "harmony", dims = 1:25) %>%
  FindNeighbors(reduction = "harmony", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 0.8, by=0.1) )
seu_clusters_TSNE <- seu %>%
  RunTSNE(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 0.8, by=0.1) ) 
clustree::clustree(seu_clusters_TSNE@meta.data[,grep("RNA_snn_res", colnames(seu_clusters_TSNE@meta.data))],
                   prefix = "RNA_snn_res.")
Clusters_tSNE <- DimPlot(seu_clusters_TSNE, reduction = 'tsne', group.by = 'RNA_snn_res.0.3', raster = FALSE, cols = mypalette3 ) 
Clusters_tSNE
BiocManager::install("SingleR")
detach("package:Matrix", unload = TRUE)
library(celldex)
library(SingleR)
seu_clusters_TSNE <- Seurat::SetIdent(seu_clusters_TSNE, value = seu_clusters_TSNE$RNA_snn_res.0.3)
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = 'tsne', "JCHAIN", label = TRUE)
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = "tsne", "IGHA1", label = TRUE)
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = "tsne", "IGHA2", label = TRUE)
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = "tsne", "S100A9", label = TRUE)
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = "tsne", "IGHG1", label = TRUE)
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = "tsne", "IGKV3-20", label = TRUE)
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = "tsne", "IGHV3-23", label = TRUE)
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = "tsne", "S100A8", label = TRUE)
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = "tsne", "TIMP1", label = TRUE)
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = 'tsne', tcell_genes, ncol=2, label = TRUE)
Seurat::VlnPlot(seu_clusters_TSNE,
                features = tcell_genes,
                ncol = 2, pt.size = 0, cols = mypalette3)
Seurat::VlnPlot(seu_clusters_TSNE,
                features = monocyte_genes,
                ncol = 2, pt.size = 0, cols = mypalette3)
seu_clusters_TSNE <- Seurat::AddModuleScore(seu_clusters_TSNE,
                                            features = list(tcell_genes),
                                            name = "tcell_genes")
Seurat::FeaturePlot(seu_clusters_TSNE, reduction = 'tsne', "tcell_genes1", label = TRUE)
Seurat::VlnPlot(seu_clusters_TSNE,
                "tcell_genes1",
                pt.size = 0, cols = mypalette3)
mypalette4 <- c("#27CED7", "#FF0000", "#FF9900")
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
seu_clusters_TSNE <- Seurat::CellCycleScoring(seu_clusters_TSNE,
                                              s.features = s.genes,
                                              g2m.features = g2m.genes)
Seurat::DimPlot(seu_clusters_TSNE, reduction = 'tsne', group.by = "Phase",
                cols = mypalette4, label = TRUE)

ref <- celldex::HumanPrimaryCellAtlasData()
class(ref)
table(ref$label.main)
seu_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu_clusters_TSNE, slot = "data"),
                                ref = ref,
                       
    BiocManager::install("dittoSeq")
SingleR::plotScoreHeatmap(seu_SingleR)
SingleR::plotDeltaDistribution(seu_SingleR)
singleR_labels <- seu_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"
seu_clusters_TSNE$SingleR_annot <- singleR_labels
install.packages(dittoseq)
mypalette0 <- c("#00796B", "#27CED7", "#FF00CC", "#00FF00", "#00FFFF", "#FF0000", "#00468B", "#FDAF91", "#5050FF", "#350E20", "#999999", "#7B4173", "#FF9900", "#358000", "#0000CC", "#99CCFF", "#FFCCCC", "#004C00", "#CCFFFF", "#CC99FF", "#9900CC", "#996600", "#666600", "#CCFF00", "#FFCC00", "#000000", "#FF420E", "#79CC3D", "#7E0021", "#FFF7F3", "#6699FF", "#CCCC99")
dittoSeq::dittoDimPlot(seu_clusters_TSNE, reduction = 'tsne', "SingleR_annot", size = 0.7)
dittoSeq::dittoBarPlot(seu_clusters_TSNE, var = "SingleR_annot", group.by = "orig.ident")
dittoSeq::dittoBarPlot(seu_clusters_TSNE, 
                       var = "SingleR_annot", 
                       group.by = "RNA_snn_res.0.3")
Seurat::VlnPlot(seu_clusters_TSNE,
                features = "percent.ribo 
 pt.size = 0, cols = mypalette3)

#########Differential gene expression###########
BiocManager::install("edgeR")
library(edgeR)
library(limma)
de_genes <- Seurat::FindAllMarkers(seu_clusters_TSNE,  min.pct = 0.25,
                                   only.pos = TRUE)
de_genes <- subset(de_genes, de_genes$p_val_adj<0.05)
write.csv(de_genes, "de_genes_FindAllMarkers.csv", row.names = F, quote = F)
write.table(de_genes, "de_genes_FindAllMarkers.txt", row.names = FALSE, quote = FALSE)
library(dplyr)
top_specific_markers <- de_genes %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)
top_specific_markers <- de_genes %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)
dittoSeq::dittoDotPlot(seu_clusters_TSNE, vars = unique(top_specific_markers$gene), 
                       group.by = "RNA_snn_res.0.3")
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
de_genes[de_genes$gene %in% tcell_genes,]
seu_clusters_TSNE <- Seurat::SetIdent(seu_clusters_TSNE, value = "SingleR_annot")
deg_im_cells <- Seurat::FindMarkers(seu_clusters_TSNE,
                                    ident.1 = "Keratinocytes",
                                    ident.2 = "Epithelial_cells",
                                    group.by = seu_clusters_TSNE$SingleR_annot,
                                    test.use = "wilcox")
deg_im_cells <- subset(deg_im_cells, deg_im_cells$p_val_adj<0.05)
View(deg_im_cells)
write.csv(deg_im_cells, "deg_im_cells.csv")

############ Matrics File Generation From Single cell analysis################ 

data <- seu@assays$RNA@layers$counts
metadata <- seu@meta.data
data1 <- (data)
str(data1)
str(metadata)
# Assuming 'data' is a sparse matrix (dgCMatrix) and 'metadata' is a data frame
# Check if 'SingleR_annot' column is present in metadata
if ('SingleR_annot' %in% colnames(metadata)) {
  # Extract the SingleR_annot column from metadata
  singleR_annot_values <- metadata$SingleR_annot
   # Identify the common row names between metadata and data
  common_row_names <- intersect(rownames(data1), rownames(metadata))
   # Replace row names in the sparse matrix with values from 'SingleR_annot'
  new_row_names <- rep("", length(rownames(data1)))
  new_row_names[match(common_row_names, rownames(data1))] <- singleR_annot_values[match(common_row_names, rownames(metadata))]
  rownames(data1) <- new_row_names
}

# Now 'data' should have row names replaced with values from 'SingleR_annot' column in 'metadata'
install.packages("scrattch.io")
remotes::install_github("AllenInstitute/scrattch.io")
devtools::install_github("AllenInstitute/scrattch.io")
library(rhdf5)
library(scrattch.io)
BiocManager::install("rhdf5")
install.packages("rhdf5")
str(data1)
write_dgCMatrix_csv(data1, "MLupdated(data).csv")



                                                       


































 Code of Machine Learning Models for Cancer Cancer Cell Identification on the base of Gene Expression
                   #############################################
                Step 3:Install Spyder from the Annaconda and Run That code 
import pandas as pd
import tensorflow as tf
import numpy as np
# Define the paths to the datasets
folder_path = "E:/ML UZAIR DATASET/"
cer_dataset_path = folder_path + "cer_dataset1.csv"
hd_dataset_path = folder_path + "HD_dataset2.csv"
oral_dataset_path = folder_path + "oral_dataset3.csv"
# Read the datasets
cer_dataset=pd.read_csv(cer_dataset_path)
hd_dataset = pd.read_csv(hd_dataset_path)
oral_dataset = pd.read_csv(oral_dataset_path)
cer_dataset.columns.values[0] = 'gene'
# Add prefix to the first column where cells are present
cer_dataset.iloc[:, 0] = cer_dataset.iloc[:, 0].apply(lambda x: 'D1_' + str(x) if pd.notnull(x) else x)
hd_dataset.iloc[:, 0] = hd_dataset.iloc[:, 0].apply(lambda x: 'D2_' + str(x) if pd.notnull(x) else x)
oral_dataset.iloc[:, 0] = oral_dataset.iloc[:, 0].apply(lambda x: 'D3_' + str(x) if pd.notnull(x) else x)
# Sort remaining columns alphabetically for each DataFrame
sorted_cer = cer_dataset.iloc[:, 1:].sort_index(axis=1)
sorted_hd = hd_dataset.iloc[:, 1:].sort_index(axis=1)
sorted_oral = oral_dataset.iloc[:, 1:].sort_index(axis=1)
# Merge the extracted first column with the sorted DataFrames
merged_cer = pd.concat([cer_dataset.iloc[:, 0], sorted_cer], axis=1)
merged_hd = pd.concat([hd_dataset.iloc[:, 0], sorted_hd], axis=1)
merged_oral = pd.concat([oral_dataset.iloc[:, 0], sorted_oral], axis=1)
print(merged_cer)
# Extract common columns
common_columns = set(merged_cer.columns[1:]).intersection(merged_hd.columns[1:], merged_oral.columns[1:])
common_columns_list = list(common_columns)
# Create new DataFrames containing only the common columns and the first column
common_cer = merged_cer[[merged_cer.columns[0]] + common_columns_list]
common_hd = merged_hd[[merged_hd.columns[0]] + common_columns_list]
common_oral = merged_oral[[merged_oral.columns[0]] + common_columns_list]
print(common_oral)
# Concatenate row by row along axis=0
merged_combined = pd.concat([common_cer, common_hd, common_oral], axis=0)
print(merged_combined)
merged_combined.to_csv('merged_combined.csv', index=False)



######### SVM  model Building and evaluation###########
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
# Split the dataset into features (X) and labels (y)
X = merged_combined.iloc[:, 1:]  # Features (all columns except the first one)
y = merged_combined.iloc[:, 0]   # Labels (first column)
# Encode the labels using LabelEncoder
label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)
# Split the data into training and testing sets (80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.2, random_state=42)
# Initialize the SVM classifier
svm_classifier = SVC(kernel='linear')
# Train the SVM classifier
svm_classifier.fit(X_train, y_train)
# Predict the labels for the test set
y_pred = svm_classifier.predict(X_test)
# Decode the predicted labels
y_pred_decoded = label_encoder.inverse_transform(y_pred)
# Calculate the accuracy of the model
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)
print("Precision:", precision)
print("Recall:", recall)
print("F1-score:", f1)
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix
# Calculate the confusion matrix
conf_matrix = confusion_matrix(y_test, y_pred)
# Normalize the confusion matrix by row (i.e., by the number of samples in each class)
conf_matrix_normalized = conf_matrix.astype('float') / conf_matrix.sum(axis=1)[:, np.newaxis]
# Plot the normalized confusion matrix with improved readability
plt.figure(figsize=(25, 25))  # Increase the figure size
ax = sns.heatmap(conf_matrix_normalized, annot=True, fmt='.2f', cmap='Blues', 
                 xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_,
                 annot_kws={"size": 8}, linewidths=.5, linecolor='black')  # Adjust font size for annotations
plt.xlabel('Predicted Labels', fontsize=14)  # Adjust font size for x-axis label
plt.ylabel('True Labels', fontsize=14)  # Adjust font size for y-axis label
plt.title('Normalized Confusion Matrix', fontsize=16)  # Adjust font size for title
plt.xticks(fontsize=10, rotation=90)  # Rotate x-axis labels for better readability
plt.yticks(fontsize=10)  # Adjust font size for y-axis labels
plt.show()


############### Random Forest Model for Cancer Classification ################

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
# Initialize the Random Forest classifier
rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)
# Train the Random Forest classifier
rf_classifier.fit(X_train, y_train)
# Predict the labels for the test set
y_pred_rf = rf_classifier.predict(X_test)
# Calculate the accuracy of the Random Forest model
accuracy_rf = accuracy_score(y_test, y_pred_rf)
print("Random Forest Accuracy:", accuracy_rf)
print("Precision:", precision)
print("Recall:", recall)
print("F1-score:", f1)
# Calculate the confusion matrix
conf_matrix_rf = confusion_matrix(y_test, y_pred_rf)
# Normalize the confusion matrix by row (i.e., by the number of samples in each class)
conf_matrix_rf_normalized = conf_matrix_rf.astype('float') / conf_matrix_rf.sum(axis=1)[:, np.newaxis]

# Plot the normalized confusion matrix with improved readability
plt.figure(figsize=(25, 25))  # Increase the figure size
ax = sns.heatmap(conf_matrix_rf_normalized, annot=True, fmt='.2f', cmap='Blues', 
                 xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_,
                 annot_kws={"size": 8}, linewidths=.5, linecolor='black')  # Adjust font size for annotations
plt.xlabel('Predicted Labels', fontsize=14)  # Adjust font size for x-axis label
plt.ylabel('True Labels', fontsize=14)  # Adjust font size for y-axis label
plt.title('Normalized Confusion Matrix for Random Forest Model', fontsize=16)  # Adjust font size for title
plt.xticks(fontsize=10, rotation=90)  # Rotate x-axis labels for better readability
plt.yticks(fontsize=10)  # Adjust font size for y-axis labels
plt.show()


# Artifical Neural Network model Build and evaluation for Cell identification on the base of Gene Expression
# ching the layers 
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.utils import to_categorical
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report
# Features and Labels
X = merged_combined.iloc[:, 1:].values  # Features
y = merged_combined.iloc[:, 0].values    # Labels
# Convert labels to categorical
label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)
y_categorical = to_categorical(y_encoded)
# Number of classes
num_classes = len(label_encoder.classes_)
# Split the data into training and testing sets (80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y_categorical, test_size=0.2, random_state=42)
# Define the model architecture with modifications
model = Sequential()
model.add(Dense(128, activation='relu', input_dim=X_train.shape[1]))  # Increase the number of neurons
model.add(Dropout(0.2))  # Add Dropout for regularization
model.add(Dense(64, activation='relu'))  # Add another hidden layer
model.add(Dense(num_classes, activation='softmax'))  # Output layer
# Compile the model with a different optimizer and learning rate
model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
# Train the model with a different number of epochs
history = model.fit(X_train, y_train, epochs=10, batch_size=32, validation_data=(X_test, y_test))
# Evaluate the model
loss, accuracy = model.evaluate(X_test, y_test)
print("Test Loss:", loss)
print("Test Accuracy:", accuracy)
# Predictions on test set
y_pred = model.predict(X_test)
# Convert predictions from categorical to labels
y_pred_labels = label_encoder.inverse_transform(y_pred.argmax(axis=1))
y_test_labels = label_encoder.inverse_transform(y_test.argmax(axis=1))
# Generate classification report
report = classification_report(y_test_labels, y_pred_labels)
print(report)
import matplotlib.pyplot as plt
# Plot training loss and validation loss
plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.title('Training and Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.show()

# Plot training accuracy and validation accuracy
plt.plot(history.history['accuracy'], label='Training Accuracy')
plt.plot(history.history['val_accuracy'], label='Validation Accuracy')
plt.title('Training and Validation Accuracy')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.legend()
plt.show()

# Xgboost Model building and implementation for cancer cells identification on the base of gene expression
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import seaborn as sns
# Assuming 'merged_combined' is already prepared
# Split X (features) and y (target)
X = merged_combined.iloc[:, 1:]  # All columns except the first one
y = merged_combined.iloc[:, 0]   # First column
# Convert categorical labels to numeric labels
label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)
# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.3, random_state=42)
# Create the XGBoost classifier
xgb_clf = xgb.XGBClassifier(use_label_encoder=False, eval_metric='mlogloss')
# Train the model
xgb_clf.fit(X_train, y_train)
# Make predictions
y_pred = xgb_clf.predict(X_test)

# Convert numeric predictions back to original labels
y_pred_labels = label_encoder.inverse_transform(y_pred)
y_test_labels = label_encoder.inverse_transform(y_test)
# Evaluate the model
accuracy = accuracy_score(y_test_labels, y_pred_labels)
print(f"Accuracy: {accuracy * 100:.2f}%")
# Confusion Matrix
conf_matrix = confusion_matrix(y_test_labels, y_pred_labels)
print("Confusion Matrix:")
print(conf_matrix)
# Plot Confusion Matrix
plt.figure(figsize=(10, 8))
sns.heatmap(conf_matrix, annot=True, fmt="d", cmap="Blues", xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_)
plt.title("Confusion Matrix")
plt.xlabel("Predicted Labels")
plt.ylabel("True Labels")
plt.show()
# Classification Report
print("Classification Report:")
class_report = classification_report(y_test_labels, y_pred_labels)
print(class_report)


