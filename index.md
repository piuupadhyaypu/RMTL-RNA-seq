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
train_cvfitc <- cvMTL(data_X_training, data_Y_training, type="Classification", Regularization="L21", Lam1_seq=10^seq(1,-4, -1),  
                Lam2=0, opts=list(init=0,  tol=10^-6, maxIter=1500), nfolds=5, stratify=FALSE, parallel=TRUE)

# Train the model
train_model=MTL(data_X_training, data_Y_training, type = "Classification", Regularization = "L21",Lam1 = train_cvfitc$Lam1.min, 
            Lam1_seq = NULL, Lam2 = 0, opts = list(init = 0, tol= 10^-3, maxIter = 100), G = NULL, k = 2)
```

**3. Validation of Classifier**
The preserved test datasets are used for evaluation of our model which is optimized by cross-validation technique. 

```
# Predictive analysis
predicted_cbmc = predict(train_model,data_X_test)
```

The figure below give us a better visualization of original with respect to predicted cell types, with the help of individual cell type plotting.

**CD14 cell type:**

![cbmc_ori_pre_CD14_cell3](https://user-images.githubusercontent.com/86721570/129381090-55a3af86-2971-4ad7-983f-5753671466ed.jpeg) 

**Mk cell type:**

![cbmc_ori_pre_Mk_cell4](https://user-images.githubusercontent.com/86721570/129381265-b7c2d987-9b2e-401c-b325-89c9b1243fe3.jpeg) 


The table shows the percentage of correct prediction and recall in each cell type of CBMC 

| Cell Type	| Sample present in data	|	Recall |	Prediction |
|-----------|:---:|:---:|:---:|:---:|
| Eryth	    | 105	        | 94      |	93.1       |	
| NK	      | 1089         |	87.77 |	94.8 |	
| CD14+ Mono |	2293 |	97.7|	99.1 |	
| Mk |	96 |	92.1 |	89.6  |	
| CD34+	| 119 |	89.04 |	88.8 |	
| DC |	70 |	91.1 |	90.8	| 
| Memory CD4 T |	1781 |	97 |	95.1 |	
| CD8 T |	273 |	90.2 |	89.7 |	
| CD16+ Mono |	230 |	87.7 |	88.5 |	
| B cell |	350 |	93.3 |	91.7	| 
| T/Mono doublets	| 182	| 92.7 |	91.5	| 
| pDCs |	49	| 91.8	| 90 |	
| Naive CD4 T	| 1248 |	98.2 |	93.6 |	



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


