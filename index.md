## RMTL based cell type detection on CBMC data

Regularized multi-task learning(RMTL) based cell type detection for single cell RNA sequence. 

### Summary

The proposed method is a supervised multi-tasking learning model to identify the marker genes from the data sequence. The result gives a clear overview of the framework. It handles the challenging goal of learning the cluster structure and treats them as a separate task in the model.

### How to use this Model
CBMC dataset used for manifestation purposes. Data can be download from https://www.ncbi.nlm.nih.gov/geo/ under accession no. GSE100866.
The column of the input matrix contains cells and genes represented in rows.

**Data Loading and Pre-processing**
Single cell RNA count matrix is loaded and then for preprocessing Seurat v3 and Limma package is used. This popular R packages supports gene and cell filtering, normalization and QC(quality control) analysis of single cell data.   
```
#load RNA UMI data
cbmc_data <- read.csv("Data/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv",header=FALSE)

#Before going into preprocess steps, prefix word 'HUMAN_' and 'MOUSE_' are removed from sequence. 

library(lattice)
library(ggplot2)
library(caret)

#Preprocessed data and cell annotations matrix used for model demontration  
data = as.matrix(read.csv("cbmc_data.csv",header=FALSE))
cell = as.matrix(read.csv("cbmc_cellannotation.csv",header=FALSE))
gene = as.matrix(read.csv("cbmc_gene.csv",header=FALSE))
```
After preprocessing step, the remaining dataset contains 2000 genes and 7895 cells. For visualization of the cell types, t-sne plot used.

![cbmc_fulldataset](https://user-images.githubusercontent.com/86721570/128569137-676b8f79-186b-4a5f-a241-776077d5e2d8.jpeg)


### Data Preparation
**1. Spliting the dataset into training and testing** 

The dataset is ramdomly split into two parts: training(80%) and testing(20%) 

**2. Construction of Classifier**

For building a classifier, we are using RMTL(Regularized Multi-Task Learning) concept. The advantage of this model is that it can take multiple cell clusters as input and simultaneously learn from the features. Before giving training data as input to the classifier, we have to arrange them into groups according to annotation matrix. 
```
# Performs Cross-Validation with parallel computing
train_cvfitc <- cvMTL(data_X, data_Y, type="Classification", Regularization="L21", Lam1_seq=10^seq(1,-4, -1),  
                Lam2=0, opts=list(init=0,  tol=10^-6, maxIter=1500), nfolds=5, stratify=FALSE, parallel=TRUE)

# Train the model
train_model=MTL(data_X, data_Y, type = "Classification", Regularization = "L21",Lam1 = train_cvfitc$Lam1.min, 
            Lam1_seq = NULL, Lam2 = 0, opts = list(init = 0, tol= 10^-3, maxIter = 100), G = NULL, k = 2)
```

**3. Validation of Classifier**
The preserved test datasets are used for evaluation of our model which is optimized by cross-validation technique. Accuracy of the model can be measured by confusion matrix and the table also shows the rate of cell type prediction.

| Cell Type	| Precision	|	Recall |	Accuracy | F1 |
|-----------|:---:|:---:|:---:|:---:|
| Eryth	    | 1	        | 1      |	1        |	1 |
| NK	      | 1         |	0.9992631 |	0.9994 |	0.9996314 |
| CD14+ Mono |	1 |	1 |	1 |	1 |
| Mk |	1 |	0.9948849 |	0.9949  |	0.9974359 |
| CD34+	| 1 |	1 |	1 |	1 |
| DC |	1 |	1 |	1	| 1 |
| Memory CD4 T |	1 |	0.9991797 |	0.9994 |	0.9995897 |
| CD8 T |	1 |	1 |	1 |	1 |
| CD16+ Mono |	1 |	0.9993455 |	0.9994 |	0.9996727 |
| B cell |	1 |	0.9993364 |	0.9994	| 0.9996681 |
| T/Mono doublets	| 1	| 0.996134 |	0.9962	| 0.9980633 |
| pDCs |	1	| 1	| 1 |	1 |
| Naive CD4 T	| 1 |	0.9992526 |	0.9994 |	0.9996262 |

This table shows, accuracy rate is 0.99, means our classifier has identified almost every cell types.

```
# Predict and Error Calculation
# Calculating error for training and testing datasets  
training_error=calcError(cbmc_train_model_cvfitr, data_cbmc_X_mtl, data_cbmc_Y_mtl)
test_error=calcError(cbmc_train_model_cvfitr, data_test_X_mtl, data_test_Y_mtl)
# Predict for classification
predicted_set_t_cbmc_classification1=predict(cbmc_train_model,data_test_X_mtl)

``````

![cbmc_web](https://user-images.githubusercontent.com/86721570/128607281-7491a4c4-85fd-490f-859c-ac252882947c.jpg)

This above plotting shows the confusion matrix heatmap plot for predicted v/s reference or ture cell types.


