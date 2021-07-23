## RMTL based cell type detection on CBMC data

Regularized multi-task learning(RMTL) based cell type detection for single cell RNA sequence. 

### Summary

The proposed method is a supervised multi-tasking learning model to identify the marker genes from the data sequence. The result gives a clear overview of the framework. It handles the challenging goal of learning the cluster structure and treats them as a separate task in the model.

### How to use this Model
CBMC dataset used for manifestation purposes. Data can be download from https://www.ncbi.nlm.nih.gov/geo/ under accession no. GSE100866.
The column of the input matrix contains cells and genes represented in rows.

**Data Loading and Pre-processing**
Single cell RNA count matrix is loaded and then cells and genes are filtered. On the filtered matrix, logarithm normalization is applied.   
```
library(lattice)
library(ggplot2)
library(caret)

```
### Data Preparation
**1. Spliting the dataset into training and testing** 

The dataset is ramdomly split into two parts: training(80%) and testing(20%) 

**2. Construction of Classifier**

For building a classifier, we are using RMTL(Regularized Multi-Task Learning) concept. The advantage of this model is that it can take multiple cell clusters as input and simultaneously learn from the features. Before giving training data as input to the classifier, we have to arrange them into groups according to annotation matrix. 
```
train_cvfitc <- cvMTL(data_X, data_Y, type="Classification", Regularization="L21", Lam1_seq=10^seq(1,-4, -1),  Lam2=0, opts=list(init=0,  tol=10^-6, maxIter=1500), nfolds=5, stratify=FALSE, parallel=TRUE)
train_model=MTL(data_X, data_Y, type = "Classification", Regularization = "L21",Lam1 = train_cvfitc$Lam1.min, Lam1_seq = NULL, Lam2 = 0, opts = list(init = 0, tol= 10^-3, maxIter = 100), G = NULL, k = 2)
```

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/piuupadhyaypu/RMTL-method-used-for-cell-type-detection-in-sc-RNA-seq/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
