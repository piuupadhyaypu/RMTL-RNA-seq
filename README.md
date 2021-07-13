# RMTL-method-used-for-cell-type-detection-in-sc-RNA-seq
Proposed method is based on regularized multi-task learning(RMTL) approach which allow us to acquire the knowledge of relationship present with a particular cell type, by identifying common features. This determines the marker genes from each cell clusters which is a tedious task. 

## Load require R package 
```
library(SingleCellExperiment)
library(RMTL)
library(Rtsne)
library(lattice)
library(caret)
```
## Preprocessing the data
goolam<-readRDS('Data/goolam.rds')
