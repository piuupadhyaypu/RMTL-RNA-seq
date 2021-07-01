####### Processing Goolam data ###########

library(kernlab)
goolam<-readRDS('C:/Users/Piu Upadhyay/Documents/goolam.rds')
write.csv(goolam@assays@.xData$data$logcounts,"goolam_data.csv") # same melanoma_datafitfinal_mat
goolam_data<-read.csv('C:/Users/Piu Upadhyay/Documents/goolam_data.csv')

rownames(goolam_data)<-goolam_data[,1]
goolam_data<-subset(goolam_data,select=-c(X)) #to remove 1st col
goolam_data_matrix=as(goolam_data,"matrix")

goolam_celltype_read<-read.csv('C:/Users/Piu Upadhyay/Documents/cell_type_goolam.csv') #cell_type of melanoma with gene_name&cell_type
cell_type_goolam=goolam_celltype_read[2]
rownames(cell_type_goolam)<-colnames(goolam_data)
Unique_cell_goolam=unique(cell_type_goolam)
goolam_cell=c("2cell","4cell","8cell","16cell","blast")

library(Rtsne)
tsne_out_goolam <- Rtsne(t(goolam_data_matrix),check_duplicates=TRUE, pca=TRUE, perplexity=30, theta=0.5, dims=2,max_iter = 5000)
d_tsne_2_goolam= as.data.frame(tsne_out_goolam$Y)
#goolam_celltype=cell_type[!is.na(match(cell_type_goolam,goolam_cell))]
true_id_goolam = as.matrix(cell_type_goolam) #same as melanoma_celltype(containing cell name and goolam_cell)
library(ggplot2)
colours<-c("2cell"="red","4cell"="blue","8cell"="brown","16cell"="green","blast"="black")
ggplot(d_tsne_2_goolam)+geom_point(aes(x=V1, y=V2,color= factor(true_id_goolam)))+ scale_color_manual(values = colours)

#######Divide training data and testing data ######################

test_size_goolam=round(0.2*ncol(goolam_data_matrix))
train_size_goolam=round(0.8*ncol(goolam_data_matrix))
test_coordinate_goolam = sample(1:ncol(goolam_data_matrix), test_size_goolam, replace=F)
test_data_goolam <- goolam_data_matrix[,c(test_coordinate_goolam)]
training_data_goolam <- goolam_data_matrix[,-c(test_coordinate_goolam)]
library(Rtsne)
tsne_training_goolam <- Rtsne(t(training_data_goolam),check_duplicates=TRUE, pca=TRUE, perplexity=30, theta=0.5, dims=2,max_iter = 5000)
d_tsne_2_training_goolam= as.data.frame(tsne_training_goolam$Y)
true_id_training_goolam = as.matrix(cell_type_goolam[-c(test_coordinate_goolam),]) 
library(ggplot2)
colours_goolam<-c("2cell"="red","4cell"="blue","8cell"="brown","16cell"="green","blast"="black")
ggplot(d_tsne_2_training_goolam)+geom_point(aes(x=V1, y=V2,color= factor(true_id_training_goolam)))+ scale_color_manual(values = colours_goolam)

cell_type_training_goolam <- as.matrix(cell_type_goolam[-c(test_coordinate_goolam),])
rownames(cell_type_training_goolam) = colnames(training_data_goolam)
cell_type_testing_goolam <- as.matrix(cell_type_goolam[c(test_coordinate_goolam),])
rownames(cell_type_testing_goolam) = colnames(test_data_goolam)

######### processing training and testing data#######       

goolam_X_mtl=list()
goolam_Y_mtl=list()
goolam_X_test=list()
goolam_Y_test=list()

for(i in 1:nrow(unique(cell_type_training_goolam)))
{
  goolam_Y_mtl[[i]]=matrix(-1,1,ncol(training_data_goolam))  
  goolam_Y_mtl[[i]][which(cell_type_training_goolam==Unique_cell_goolam$cell_type1[i])]=1
}

for(i in 1:nrow(unique(cell_type_testing_goolam)))
{
  goolam_Y_test[[i]]=matrix(-1,1,ncol(test_data_goolam))  
  goolam_Y_test[[i]][which(cell_type_testing_goolam==Unique_cell_goolam$cell_type1[i])]=1
}

#for(i in 1:nrow(unique(cell_type_goolam)))
#{
#  goolam_Y_mtl[[i]]=matrix(-1,1,ncol(goolam_data))
#  goolam_Y_mtl[[i]][which(cell_type_goolam==Unique_cell_goolam$cell_type1[i])]=1
#}

for(i in 1:nrow(unique(cell_type_training_goolam)))
{
  goolam_X_mtl[[i]]=t(training_data_goolam)
  row_no1=which(cell_type_training_goolam!=Unique_cell_goolam$cell_type1[i])
  goolam_X_mtl[[i]][row_no1,]=0
}

for(i in 1:nrow(unique(cell_type_testing_goolam)))
{
  goolam_X_test[[i]]=t(test_data_goolam)
  row_no2=which(cell_type_testing_goolam!=Unique_cell_goolam$cell_type1[i])
  goolam_X_test[[i]][row_no2,]=0
}

#for(i in 1:nrow(unique(cell_type_goolam)))
#{
#  goolam_X_mtl[[i]]=t(goolam_data)
#  row_no=which(cell_type_goolam!=Unique_cell_goolam$cell_type1[i])
#  goolam_X_mtl[[i]][row_no,]=0
#}

for(i in 1:nrow(unique(cell_type_training_goolam)))
{
  goolam_Y_mtl[[i]]=t(goolam_Y_mtl[[i]])
}

for(i in 1:nrow(unique(cell_type_testing_goolam)))
{
  goolam_Y_test[[i]]=t(goolam_Y_test[[i]])
}

# training data
library(RMTL)
Sys.time()

goolam_train_cvfitc <- cvMTL(goolam_X_mtl, goolam_Y_mtl, type="Classification", Regularization="L21", Lam1_seq=10^seq(1,-4, -1),  Lam2=0, opts=list(init=0,  tol=10^-6, maxIter=1500), nfolds=5, stratify=FALSE, parallel=TRUE)
goolam_train_model=MTL(goolam_X_mtl, goolam_Y_mtl, type = "Classification", Regularization = "L21",Lam1 = goolam_train_cvfitc$Lam1.min, Lam1_seq = NULL, Lam2 = 0, opts = list(init = 0, tol= 10^-3, maxIter = 100), G = NULL, k = 2)

############
str(goolam_train_cvfitc)
plot(goolam_train_cvfitc)
training_error_goolam=calcError(goolam_train_model, goolam_X_mtl, goolam_Y_mtl)
test_error_goolam=calcError(goolam_train_model, goolam_X_test, goolam_Y_test)
paste0("training error: ", calcError(goolam_train_model, goolam_X_mtl, goolam_Y_mtl))
paste0("test error: ", calcError(goolam_train_model, goolam_X_test, goolam_Y_test))

#################
library(Rtsne)
tsne_testing_goolam <- Rtsne(t(test_data_goolam),check_duplicates=TRUE, pca=TRUE, perplexity=8, theta=0.5, dims=2,max_iter = 5000)
d_tsne_2_testing_goolam= as.data.frame(tsne_testing_goolam$Y)
true_id_testing_goolam = as.matrix(cell_type_goolam[c(test_coordinate_goolam),]) 
library(ggplot2)
ggplot(d_tsne_2_testing_goolam)+geom_point(aes(x=V1, y=V2,color= factor(true_id_testing_goolam)))+ scale_color_manual(values = colours_goolam)

#######prediction ############
predicted_set_t_goolam_classification=predict(goolam_train_model,goolam_X_test)

plot(density(predicted_set_t_goolam_classification[[4]]))

#overall threshold 0.7 on classification
for(i in 1:4)
{
  for(j in 1:nrow(predicted_set_t_goolam_classification[[i]]))
  {
    if((predicted_set_t_goolam_classification[[i]][j])>=0.7)
    {
      predicted_set_t_goolam_classification[[i]][j]=1
    }
    else
      predicted_set_t_goolam_classification[[i]][j]=0
  }
}
##########################
# Original vs Predicted plot
i=4
trueid_t_goolam=matrix(,nrow(cell_type_testing_goolam),1)
trueid_t_goolam[cell_type_testing_goolam==goolam_cell[i]]=goolam_cell[i]
predicted_id_t_goolam=matrix(,nrow(cell_type_testing_goolam),1)
order_predicted_goolam=which(predicted_set_t_goolam_classification[[i]]==1)
predicted_id_t_goolam[order_predicted_goolam]=goolam_cell[i]
library(cowplot)
par(mfrow=c(2,1))
plot1<-ggplot(d_tsne_2_testing_goolam, aes(x=V1, y=V2))+geom_point(col="grey")+geom_point(aes(col=predicted_id_t_goolam))+  scale_colour_manual(values = colours_goolam)+xlab("t-SNE dim-1")+ylab("t-SNE dim-2")+ theme_bw()+theme(legend.position = "none")
plot2<-ggplot(d_tsne_2_testing_goolam, aes(x=V1, y=V2))+geom_point(col="grey")+geom_point(aes(col=trueid_t_goolam))+  scale_colour_manual(values = colours_goolam)+xlab("t-SNE dim-1")+ylab("t-SNE dim-2")+ theme_bw()+theme(legend.position = "none")
plot_grid(plot2,plot1)

##########################
#Data preparation for confusion matrix

confusion_goolam_Y_test=list()
confusion_goolam_Y_test=goolam_Y_test
for (i in 1:4) 
{
  for(j in 1:nrow(confusion_goolam_Y_test[[i]]))
  {
    if((confusion_goolam_Y_test[[i]][j])== -1) 
    {
      confusion_goolam_Y_test[[i]][j]=0
    }
  }
  
}
library(lattice)
library(ggplot2)
library(caret)

con_mat_goolam<-confusionMatrix(factor(predicted_id_t_goolam),factor(cell_type_testing_goolam))

ggplot(as.data.frame(con_mat_goolam$table), aes(Prediction,Reference,fill= Freq))+geom_tile() + geom_text(aes(label=Freq)) + scale_fill_gradient(low="green", high="#009194") + labs(x = "Reference",y = "Prediction") 

ggplot(as.data.frame(con_mat_goolam$table), aes(Prediction,Reference,fill= Freq))+geom_tile() + scale_fill_gradient(low="light blue", high="purple") +labs(x = "Reference",y = "Prediction")+theme(axis.text.x = element_text(angle = 90))

###########################
F1_goolam=array()
TP_goolam=array()
precision_goolam=array()
recall_goolam=array()
predicted_idc=matrix(,nrow(cell_type_goolam),1)
for(i in 1:4){
  
  o=order(predicted_set_t_goolam_classification[[i]],decreasing = TRUE)
  uniq_cell=unique(cell_type_goolam)
  predicted_idc[o[1:(sum(cell_type_goolam==uniq_cell[i,]))]]=uniq_cell[i,]
  
  TP_goolam[i]=length(intersect(which(predicted_idc==uniq_cell[i,]),which(cell_type_goolam==uniq_cell[i,])))
  
  precision_goolam[i]=TP_goolam[i]/sum(!is.na(predicted_idc))
  
  recall_goolam[i]=TP_goolam[i]/sum(cell_type_goolam==uniq_cell[i,])
  F1_goolam[i] <- (2 * precision_goolam[i] * recall_goolam[i]) / (precision_goolam[i] + recall_goolam[i])
}

#########donuts chart for all true cell and all predicted cell#############

##for all true cell Goolam
library(ggpubr)
p=array()
for( i in 1:length(goolam_cell)){
  p[i]=length(which(cell_type_testing_goolam==goolam_cell[i]))/length(cell_type_testing_goolam)
}
df1=data.frame(val=100*p,cell=goolam_cell)
df1$cell=paste(df1$cell, " ",round(df1$val,digits =2), "%",sep="")

v=match(sort(df1$cell),unique(cell_type_goolam))
myPalette=getColors_goolam(v)
#df_order=df[order(df$cell),]
#labs=paste(df_order$cell, " (", df_order$val, "%)")
ggdonutchart(df1, "val", label = "cell",lab.font = "white",lab.pos="out", fill = "cell",palette = getColors_goolam(6))+theme(axis.text.x = element_text(face = "bold", color = "#993333", size = 10, angle = 60),legend.position ="top",legend.direction = "horizontal")

##for prdicted cell

predicted_id=matrix(,nrow(cell_type_testing_goolam),1)
for(i in 1:4){
  #predicted_id_t_cbmc=matrix(,nrow(cell_type_testing_cbmc),1)
  order1=which(predicted_set_t_goolam_classification[[i]]==1)
  predicted_id[order1]=goolam_cell[i]
}
ggplot(d_tsne_2_testing_goolam, aes(x=V1, y=V2))+geom_point(col="grey")+geom_point(aes(col=predicted_id))+  scale_colour_manual(values = getColors_goolam(6))+xlab("t-SNE dim-1")+ylab("t-SNE dim-2")+ theme_bw()+theme(legend.position = "none")

# predicted donut chart
p1=array()
for( i in 1:length(goolam_cell)){
  p1[i]=length(which(predicted_id==goolam_cell[i]))/length(predicted_id)
}
df=data.frame(val=100*p1,cell=goolam_cell)
df$cell=paste(df$cell, " ",round(df$val,digits =2), "%",sep="")
ggdonutchart(df, "val", label = "cell",lab.font = "white",lab.pos="out", fill = "cell",palette = getColors_goolam(6))+theme(axis.text.x = element_text(face = "bold", color = "#993333", size = 10, angle = 60),legend.position ="top",legend.direction = "horizontal")

###################################
#########color function Goolam ###########
getColors_goolam<-function(n){
  
  mycolors = c("#19FFFF", "#FFFF33", "#40E0D0" , "#FFA500" , "#000000" , "#964B00" , "#FF0000" , "#FF69B4"  , "#66FF00" , "#A020F0" , "#0000FF" , "#006A4E","#FF00FF")
  if(n>20){
    cat("Too many colors...Using fallback color scheme.\n")
    return(getPalette(n))
  }
  return(mycolors[1:n])
  
}