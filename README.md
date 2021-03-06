# RMTL-RNA-seq
Regularized multi-task learning(RMTL) based cell type detection for single cell RNA sequence. 

## Introduction
Proposed method is based on regularized multi-task learning(RMTL) approach which allow us to acquire the knowledge of relationship present with a particular cell type, by identifying common features. This determines the marker genes from each cell clusters which is a tedious task, if we do it manually.  

## Load require R package 
```
library(Seurat)
library(Limma)
library(SingleCellExperiment)
library(RMTL)
library(Rtsne)
library(lattice)
library(caret)
```
## Preprocessing the data
The data used here is preprocessed by Seurat v3 and Limma package for prior analysis, quality control, preprocessing like gene and cell filtering, normalization. 

## Dry run on CBMC data
[Demo run of RMTL based cell type detection on CBMC data](https://github.com/piuupadhyaypu/RMTL-RNA-seq/blob/fded1f9e8bc524229122f2c64854ed4993f09910/index.md)

