####### Processing CBMC data ###########

rna_scaled<-read.csv('C:/Users/Piu Upadhyay/Downloads/rna_scaled (1)/rna_scaled.csv')
cbmc_read<-read.csv('C:/Users/Piu Upadhyay/Downloads/cell_type_cbmc.csv')
cbmc_read<-subset(cbmc_read, select = -c(X.1))
cell_type_cbmc= cbmc_read[2]
rownames(cell_type_cbmc)<-cbmc_read[,1]
Unique_cell=unique(cell_type_cbmc)
cbmc_cell=c("Eryth","NK","CD14+ Mono","Mk","CD34+","DC","Memory CD4 T","CD8 T","CD16+ Mono","B","T/Mono doublets","pDCs","Naive CD4 T")

rownames(rna_scaled)<- rna_scaled[,1]
rna_scaled<-subset(rna_scaled,select=-c(X)) #to remove 1st col
cbmc_data_matrix=as(rna_scaled,"matrix")

library(Rtsne)
tsne_out_cbmc <- Rtsne(t(cbmc_data_matrix),check_duplicates=TRUE, pca=TRUE, perplexity=30, theta=0.5, dims=2,max_iter = 5000)
d_tsne_2_cbmc= as.data.frame(tsne_out_cbmc$Y)
true_id_cbmc = as.matrix(cell_type_cbmc) 
library(ggplot2)
colours_cbmc<-c("Eryth"="red","NK"="blue","CD14+ Mono"="yellow","Mk"="green","CD34+"="orange","DC"="brown","Memory CD4 T"="hot pink","CD8 T"="black","CD16+ Mono"="turquoise","B"="cyan","T/Mono doublets"="magenta","pDCs"="dark green","Naive CD4 T"="purple")
ggplot(d_tsne_2_cbmc)+geom_point(aes(x=V1, y=V2,color= factor(true_id_cbmc)))+ scale_color_manual(values = colours_cbmc)

######## Divide training data and testing data #######

test_size=round(0.2*ncol(cbmc_data_matrix))
train_size=round(0.8*ncol(cbmc_data_matrix))
test_coordinate = sample(1:ncol(cbmc_data_matrix), test_size, replace=F)
test_data_cbmc <- cbmc_data_matrix[,c(test_coordinate)]
training_data_cbmc <- cbmc_data_matrix[,-c(test_coordinate)]
library(Rtsne)
tsne_training_cbmc <- Rtsne(t(training_data_cbmc),check_duplicates=TRUE, pca=TRUE, perplexity=30, theta=0.5, dims=2,max_iter = 5000)
d_tsne_2_training_cbmc= as.data.frame(tsne_training_cbmc$Y)
true_id_training_cbmc = as.matrix(cell_type_cbmc[-c(test_coordinate),]) 
library(ggplot2)
colours_cbmc<-c("Eryth"="red","NK"="blue","CD14+ Mono"="yellow","Mk"="green","CD34+"="orange","DC"="brown","Memory CD4 T"="hot pink","CD8 T"="black","CD16+ Mono"="turquoise","B"="cyan","T/Mono doublets"="magenta","pDCs"="dark green","Naive CD4 T"="purple")
ggplot(d_tsne_2_training_cbmc)+geom_point(aes(x=V1, y=V2,color= factor(true_id_training_cbmc)))+ scale_color_manual(values = colours_cbmc)

############################################

cell_type_training_cbmc <- as.matrix(cell_type_cbmc[-c(test_coordinate),])
rownames(cell_type_training_cbmc) = colnames(training_data_cbmc)
cell_type_testing_cbmc <- as.matrix(cell_type_cbmc[c(test_coordinate),])
rownames(cell_type_testing_cbmc) = colnames(test_data_cbmc)

#############################################
#processing training n testing data set
data_cbmc_X_mtl=list()
data_cbmc_Y_mtl=list()
data_test_X_mtl=list()
data_test_Y_mtl=list()
for(i in 1:nrow(unique(cell_type_training_cbmc)))
{
  data_cbmc_Y_mtl[[i]]=matrix(-1,1,ncol(training_data_cbmc))  
  data_cbmc_Y_mtl[[i]][which(cell_type_training_cbmc==Unique_cell$x[i])]=1
}

for(i in 1:nrow(unique(cell_type_testing_cbmc)))
{
  data_test_Y_mtl[[i]]=matrix(-1,1,ncol(test_data_cbmc))  
  data_test_Y_mtl[[i]][which(cell_type_testing_cbmc==Unique_cell$x[i])]=1
}

for(i in 1:nrow(unique(cell_type_training_cbmc)))
{
  data_cbmc_X_mtl[[i]]=t(training_data_cbmc)
  row_no1=which(cell_type_training_cbmc!=Unique_cell$x[i])
  data_cbmc_X_mtl[[i]][row_no1,]=0
}

for(i in 1:nrow(unique(cell_type_testing_cbmc)))
{
  data_test_X_mtl[[i]]=t(test_data_cbmc)
  row_no2=which(cell_type_testing_cbmc!=Unique_cell$x[i])
  data_test_X_mtl[[i]][row_no2,]=0
}

for(i in 1:nrow(unique(cell_type_training_cbmc)))
{
  data_cbmc_Y_mtl[[i]]=t(data_cbmc_Y_mtl[[i]])
}

for(i in 1:nrow(unique(cell_type_testing_cbmc)))
{
  data_test_Y_mtl[[i]]=t(data_test_Y_mtl[[i]])
}

# training data
Sys.time()

cbmc_train_cvfitc <- cvMTL(data_cbmc_X_mtl, data_cbmc_Y_mtl, type="Classification", Regularization="L21", Lam1_seq=10^seq(1,-4, -1),  Lam2=0, opts=list(init=0,  tol=10^-6, maxIter=1500), nfolds=5, stratify=FALSE, parallel=TRUE)
cbmc_train_model=MTL(data_cbmc_X_mtl, data_cbmc_Y_mtl, type = "Classification", Regularization = "L21",Lam1 = cbmc_train_cvfitc$Lam1.min, Lam1_seq = NULL, Lam2 = 0, opts = list(init = 0, tol= 10^-3, maxIter = 100), G = NULL, k = 2)

str(cbmc_train_cvfitr)
plot(cbmc_train_cvfitr)
training_error=calcError(cbmc_train_model_cvfitr, data_cbmc_X_mtl, data_cbmc_Y_mtl)
test_error=calcError(cbmc_train_model_cvfitr, data_test_X_mtl, data_test_Y_mtl)
paste0("training error: ", calcError(cbmc_train_model_cvfitr, data_cbmc_X_mtl, data_cbmc_Y_mtl))
paste0("test error: ", calcError(cbmc_train_model_cvfitr, data_test_X_mtl, data_test_Y_mtl))

paste0("training error: ", calcError(cbmc_train_model, data_cbmc_X_mtl, data_cbmc_Y_mtl))
paste0("test error: ", calcError(cbmc_train_model, data_test_X_mtl, data_test_Y_mtl))

library(Rtsne)
tsne_testing_cbmc <- Rtsne(t(test_data_cbmc),check_duplicates=TRUE, pca=TRUE, perplexity=30, theta=0.5, dims=2,max_iter = 5000)
d_tsne_2_testing_cbmc= as.data.frame(tsne_testing_cbmc$Y)
true_id_testing_cbmc = as.matrix(cell_type_cbmc[c(test_coordinate),]) 
library(ggplot2)
colours_cbmc<-c("Eryth"="red","NK"="blue","CD14+ Mono"="yellow","Mk"="green","CD34+"="orange","DC"="brown","Memory CD4 T"="hot pink","CD8 T"="black","CD16+ Mono"="turquoise","B"="cyan","T/Mono doublets"="magenta","pDCs"="dark green","Naive CD4 T"="purple")
ggplot(d_tsne_2_testing_cbmc)+geom_point(aes(x=V1, y=V2,color= factor(true_id_testing_cbmc)))+ scale_color_manual(values = colours_cbmc)
##############
#Predict
predicted_set_t_cbmc_classification1=predict(cbmc_train_model,data_test_X_mtl)

predicted_set_t_cbmc_classification=predict(cbmc_train_model,data_test_X_mtl)


#display density of each cell type
plot(density(predicted_set_t_cbmc_regression[[9]]))
plot(density(predicted_set_t_cbmc_classification[[13]]))


#overall threshold 0.8 on classification
for(i in 1:13)
{
  for(j in 1:nrow(predicted_set_t_cbmc_classification[[i]]))
  {
    if((predicted_set_t_cbmc_classification[[i]][j])>=0.8)
    {
      predicted_set_t_cbmc_classification[[i]][j]=1
    }
    else
      predicted_set_t_cbmc_classification[[i]][j]=0
  }
}
######################
#Data preparation for confusion matrix
confusion_cbmc_Y_test=list()
confusion_cbmc_Y_test=data_test_Y_mtl
for (i in 1:13) 
{
  for(j in 1:nrow(confusion_cbmc_Y_test[[i]]))
  {
    if((confusion_cbmc_Y_test[[i]][j])== -1) 
    {
      confusion_cbmc_Y_test[[i]][j]=0
    }
  }
  
}
library(lattice)
library(caret)
i=2
print(confusion_matrix<-confusionMatrix(factor(confusion_cbmc_Y_test[[i]]),factor(predicted_set_t_cbmc_classification[[i]])))
print(precision<-confusion_matrix$table[1]/(confusion_matrix$table[1]+confusion_matrix$table[3]))
print(recall<-confusion_matrix$table[1]/(confusion_matrix$table[1]+confusion_matrix$table[2]))
print(F1<-(2*confusion_matrix$table[1])/((2*confusion_matrix$table[1])+confusion_matrix$table[2]+confusion_matrix$table[3]))


con_mat<-confusionMatrix(factor(predicted_id_t_cbmc),factor(cell_type_testing_cbmc))

ggplot(as.data.frame(con_mat$table), aes(Prediction,Reference,fill= Freq))+geom_tile() + scale_fill_gradient(low="blue", high="purple") +labs(x = "Reference",y = "Prediction")+theme(axis.text.x = element_text(angle = 90))
table_con_mat = as.data.frame(con_mat$table)
table_con_mat$Y1 <- cut(table_con_mat$Freq,breaks = c(-Inf,0,50,100,150,200,250,300,350,400,Inf),right = FALSE)
for (i in 1:nrow(table_con_mat)) {
  if(table_con_mat$Freq[i] == 0)
  {table_con_mat$Y1[i]="[-Inf,0)"}
}
ggplot(table_con_mat, aes(Prediction,Reference))+ geom_tile(aes(fill = Y1)) +labs(x = "Reference",y = "Prediction")+theme(axis.text.x = element_text(angle = 90))+ scale_fill_manual(breaks=c("[400, Inf)","[350,400)","[300,350)","[250,300)","[200,250)","[150,200)","[100,150)","[50,100)","[0,50)","[-Inf,0)"), values = c("grey","purple","dark blue", "yellow", "blue", "lightblue", "lightgreen", "green","darkgreen", "pink","dark blue","white"))

 ###########################
i=13
trueid_t_cbmc=matrix(,nrow(cell_type_testing_cbmc),1)
trueid_t_cbmc[cell_type_testing_cbmc==cbmc_cell[i]]=cbmc_cell[i]
predicted_id_t_cbmc=matrix(,nrow(cell_type_testing_cbmc),1)
order_predicted_cbmc=which(predicted_set_t_cbmc_classification[[i]]==1)
predicted_id_t_cbmc[order_predicted_cbmc]=cbmc_cell[i]
library(cowplot)
par(mfrow=c(2,1))
plot1<-ggplot(d_tsne_2_testing_cbmc, aes(x=V1, y=V2))+geom_point(col="grey")+geom_point(aes(col=predicted_id_t_cbmc))+  scale_colour_manual(values = colours_cbmc)+xlab("t-SNE dim-1")+ylab("t-SNE dim-2")+ theme_bw()+theme(legend.position = "none")+labs(title = "Predicted")
plot2<-ggplot(d_tsne_2_testing_cbmc, aes(x=V1, y=V2),color= factor(trueid_t_cbmc))+geom_point(col="grey")+geom_point(aes(col=trueid_t_cbmc))+  scale_colour_manual(values = colours_cbmc)+xlab("t-SNE dim-1")+ylab("t-SNE dim-2")+ theme_bw()+theme(legend.position = "none")+labs(title = "Original")
plot_grid(plot2,plot1)

###############heatmap
install.packages("remotes")
remotes::install_github("dviraran/SingleR")
####################
i=2
trueid_t_cbmc=matrix(,nrow(cell_type_testing_cbmc),1)
trueid_t_cbmc[cell_type_testing_cbmc==cbmc_cell[i]]=cbmc_cell[i]

predicted_id_t_cbmc=matrix(,nrow(cell_type_testing_cbmc),1)
order_predicted_cbmc=which(predicted_set_t_cbmc_classification[[i]]==1)
predicted_id_t_cbmc[order_predicted_cbmc]=cbmc_cell[i]

common_cbmc=intersect(which(predicted_id_t_cbmc==cbmc_cell[i]),which(trueid_t_cbmc==cbmc_cell[i]))
non_match_predicted_cbmc=setdiff(which(predicted_id_t_cbmc==cbmc_cell[i]),common_cbmc)
#non_match_predicted_cbmc=setdiff(which(trueid_t_cbmc==cbmc_cell[i]),which(predicted_id_t_cbmc==cbmc_cell[i]))

Acc_cbmc=length(common_cbmc)/sum(cell_type_testing_cbmc==cbmc_cell[i])
print(Acc_cbmc)

wrong_predicted_celltype_cbmc=cell_type_testing_cbmc[non_match_predicted_cbmc,]
uniq_wrong_predicted_celltype_cbmc= unique(wrong_predicted_celltype_cbmc)
percentage_wrong_predicted_celltype_cbmc=array()
for(j in 1:length(unique(wrong_predicted_celltype_cbmc))){
  percentage_wrong_predicted_celltype_cbmc[j]=length(which(wrong_predicted_celltype_cbmc==uniq_wrong_predicted_celltype_cbmc[j]))/sum(cell_type_testing_cbmc==cbmc_cell[j])
}
print(percentage_wrong_predicted_celltype_cbmc)

df_cbmc=data.frame(val=100*append(Acc_cbmc,percentage_wrong_predicted_celltype_cbmc),cell=append(cbmc_cell[i],uniq_wrong_predicted_celltype_cbmc))
##creating pi-chart

pie_cbmc=ggplot(df_cbmc, aes(x = "", y=val,fill=cell)) + geom_bar(width = 1,stat="identity")+coord_polar("y", start=0)

ggplot(df_cbmc, aes(x = 2, y = val, fill = cell)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = val, label = cell), color = "white")+
  theme_void()+
  xlim(0.5, 2.5)

v=match(sort(df_cbmc$cell),cbmc_cell)
myPalette=getColors(v)
ggdonutchart(df_cbmc, "val", label = "cell",lab.font = "white",lab.pos="out", fill = "cell",palette = myPalette)+theme(axis.text.x = element_text(face = "bold", color = "#993333", size = 10, angle = 30,hjust=0.95,vjust=0.9))

#########donuts chart for all true cell and all predicted cell#############

##for all true cell CBMC

p=array()
for( i in 1:length(cbmc_cell)){
  p[i]=length(which(cell_type_testing_cbmc==cbmc_cell[i]))/length(cell_type_testing_cbmc)
}
df1=data.frame(val=100*p,cell=cbmc_cell)
df1$cell=paste(df1$cell, " ",round(df1$val,digits =2), "%",sep="")
cbmc_cell=c("Eryth","NK","CD14+ Mono","Mk","CD34+","DC","Memory CD4 T","CD8 T","CD16+ Mono","B","T/Mono doublets","pDCs","Naive CD4 T")
v=match(sort(df1$cell),unique(cell_type_cbmc))
myPalette=getColors_cbmc(v)
#df_order=df[order(df$cell),]
#labs=paste(df_order$cell, " (", df_order$val, "%)")
ggdonutchart(df1, "val", label = "cell",lab.font = "white",lab.pos="out", fill = "cell",palette = getColors_cbmc(13))+theme(axis.text.x = element_text(face = "bold", color = "#993333", size = 10, angle = 60),legend.position ="top",legend.direction = "horizontal")

##for prdicted cell

predicted_id=matrix(,nrow(cell_type_testing_cbmc),1)
for(i in 1:13){
  #predicted_id_t_cbmc=matrix(,nrow(cell_type_testing_cbmc),1)
  order1=which(predicted_set_t_cbmc_classification[[i]]==1)
  predicted_id[order1]=cbmc_cell[i]
}
ggplot(d_tsne_2_testing_cbmc, aes(x=V1, y=V2))+geom_point(col="grey")+geom_point(aes(col=predicted_id))+  scale_colour_manual(values = getColors_cbmc(13))+xlab("t-SNE dim-1")+ylab("t-SNE dim-2")+ theme_bw()+theme(legend.position = "none")

#
p1=array()
for( i in 1:length(cbmc_cell)){
  p1[i]=length(which(predicted_id==cbmc_cell[i]))/length(predicted_id)
}
df=data.frame(val=100*p1,cell=cbmc_cell)
df$cell=paste(df$cell, " ",round(df$val,digits =2), "%",sep="")
ggdonutchart(df, "val", label = "cell",lab.font = "white",lab.pos="out", fill = "cell",palette = getColors_cbmc(13))+theme(axis.text.x = element_text(face = "bold", color = "#993333", size = 10, angle = 60),legend.position ="top",legend.direction = "horizontal")

###################################

#######################
###############analyzing precision and recall for each cell type

non_match_predicted_cbmc=setdiff(which(seperate_predictedid_precision==cbmc_cell[i]),common_cbmc)

Acc_cbmc=length(common_cbmc)/sum(cell_type_cbmc==cbmc_cell[i])
print(Acc_cbmc)
wrong_predicted_celltype_cbmc=cell_type_cbmc[non_match_predicted_cbmc,]
uniq_wrong_predicted_celltype_cbmc= unique(wrong_predicted_celltype_cbmc)
percentage_wrong_predicted_celltype_cbmc=array()
for(j in 1:length(unique(wrong_predicted_celltype_cbmc[[i]]))){
  percentage_wrong_predicted_celltype_cbmc[j]=length(which(wrong_predicted_celltype_cbmc==uniq_wrong_predicted_celltype_cbmc[j]))/sum(cell_type_cbmc==cbmc_cell[j])
}
print(percentage_wrong_predicted_celltype_cbmc)
percentage_wrong_p_cell= percentage_wrong_predicted_celltype_cbmc
df_cbmc=data.frame(val=100*append(Acc_cbmc,percentage_wrong_predicted_celltype_cbmc),cell=append(cbmc_cell[i],uniq_wrong_predicted_celltype_cbmc))

##creating pi-chart

library(ggplot2)
pie=ggplot(df_cbmc, aes(x = "", y=val,fill=cell)) + geom_bar(width = 1,stat="identity")+coord_polar("y", start=0)+geom_text(aes(y = val, label = cell), color = "white")


ggplot(df_cbmc, aes(x = 2, y = val, fill = cell)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = val, label = cell), color = "white")+
  theme_void()+
  xlim(0.5, 2.5)

F1=array()
TP_c=array()
precision_c=array()
recall_c=array()
predicted_idc=matrix(,nrow(cell_type_cbmc),1)
for(i in 1:13){
  o= which(predicted_set_t_cbmc_classification[[i]]==1)
  #o=order(predicted_set_t_cbmc_classification[[i]],decreasing = TRUE)
  uniq_cell=unique(cell_type_testing_cbmc)
  predicted_idc[o[1:(sum(cell_type_testing_cbmc==uniq_cell[i,]))]]=uniq_cell[i,]
  
  TP_c[i]=length(intersect(which(predicted_idc==uniq_cell[i,]),which(cell_type_testing_cbmc==uniq_cell[i,])))
  
  precision_c[i]=TP_c[i]/sum(!is.na(predicted_idc))
  
  recall_c[i]=TP_c[i]/sum(cell_type_testing_cbmc==uniq_cell[i,])
  F1[i] <- (2 * precision_c[i] * recall_c[i]) / (precision_c[i] + recall_c[i])
}

#########color function CBMC ###########
getColors_cbmc<-function(n){
  #mycolors = c("#00fe0a", "#794b05","#fc9a07" , "#000000", "#f6fc2a",
  #             "#0000ff", "#00ffff","#f87791" ,"#81b7dd" ,
  #             "#c3beb7", "#1e7309", "#625b51", "#6a09c3", "#189ff5",
  #             "#d19d00", "#0ebf06", "#88ffb3", "#bded1b", "#ff0000", "ff21d3")
mycolors = c("#19FFFF", "#FFFF33", "#40E0D0" , "#FFA500" , "#000000" , "#964B00" , "#FF0000" , "#FF69B4"  , "#66FF00" , "#A020F0" , "#0000FF" , "#006A4E","#FF00FF")
  if(n>20){
    cat("Too many colors...Using fallback color scheme.\n")
    return(getPalette(n))
  }
  return(mycolors[1:n])
  
}
############donout chart################ 
v=match(sort(df_cbmc$cell),cbmc_cell)
myPalette=getColors_cbmc(v)
library(ggpubr)
ggdonutchart(df_cbmc, "val", label = "cell",lab.font = "white",lab.pos="out", fill = "cell",palette = getColors_cbmc(13))+theme(axis.text.x = element_text(face = "bold", color = "#993333", angle = 30,hjust=0.95,vjust=0.9))
