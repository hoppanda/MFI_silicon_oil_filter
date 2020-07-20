########################################################################
###
###    Silicon Oil Image filter - Raw Image-based Filters
###
###    Gregory, Chen (May 2020)
###
########################################################################


#: ------------------------------------------------------
#:  Pre-setup
#: ------------------------------------------------------

# library(jpeg)
# library(imager)
# source("http://bioconductor.org/biocLite.R")
# biocLite("EBImage")
# library(tidyverse)
# library(DepthProc)
library(EBImage)
library(mvtnorm)
library(cluster)
library(factoextra)
library(ddalpha)

#: ------------------------------------------------------
#:  Functions
#: ------------------------------------------------------
chop_blank = function(tpic){
  ch1a=apply(tpic,1,mean)
  ch1b=apply(tpic,1,function(x) max(x)-min(x))
  ch2a=apply(tpic,2,mean)
  ch2b=apply(tpic,2,function(x) max(x)-min(x))
  row.cut=!(ch1a>0.78 & ch1b<0.28)
  col.cut=!(ch2a>0.78 & ch2b<0.28)
  tpic[row.cut,col.cut]
}

#: ------------------------------------------------------
#:  Loading Data
#: ------------------------------------------------------
IA_train_path = "H:/data/Projects/Research/MFI/Data/Filter Picture/"
img_files = list.files(IA_train_path)

# -> Read in files
raw_imgs = lapply(img_files,function(i) readImage(file=paste0(IA_train_path,i)))
# -> Info on the size of the files
raw_size = t(sapply(raw_imgs, function(x) dim(x)))
raw_size = as.data.frame(raw_size)
names(raw_size) = c("nrow","ncol")
raw_size$name = gsub('.jpg','',img_files)
raw_size$ID = paste(raw_size$name,paste(raw_size$nrow,raw_size$ncol,sep="*"),sep="_")
names(raw_imgs) = raw_size$ID
# -> Rescale them to 20 by 20, with augmentation of 90 degree rotation 3 times
imgs_20_a = lapply(raw_imgs,function(temp){
  temp=(temp-min(temp))/(max(temp-min(temp)))
  temp=chop_blank(temp)
  temp=(temp-min(temp))/(max(temp-min(temp)))
  temp=chop_blank(temp)
  temp=(temp-min(temp))/(max(temp-min(temp)))
  if(abs(nrow(temp)-ncol(temp))>1){
    if(ncol(temp)<nrow(temp)){
      pad_val=(max(temp[,1])+max(temp[,ncol(temp)]))/2
      nr_pad=nrow(temp)-ncol(temp)
      pad_mat=matrix(pad_val,nrow=nrow(temp),ncol=floor(nr_pad/2))
      temp=cbind(pad_mat,temp,pad_mat)
    }else{
      pad_val=(max(temp[1,])+max(temp[nrow(temp),]))/2
      nr_pad=ncol(temp)-nrow(temp)
      pad_mat=matrix(pad_val,nrow=floor(nr_pad/2),ncol=ncol(temp))
      temp=rbind(pad_mat,temp,pad_mat)
    }
  }
  resize(temp,20,20)
})
imgs_20_b = lapply(imgs_20_a,function(temp) rotate(temp,90))
names(imgs_20_b)=paste0(names(imgs_20_a),"_R90")
imgs_20_c = lapply(imgs_20_a,function(temp) rotate(temp,180))
names(imgs_20_c)=paste0(names(imgs_20_a),"_R180")
imgs_20_d = lapply(imgs_20_a,function(temp) rotate(temp,270))
names(imgs_20_d)=paste0(names(imgs_20_a),"_R270")

imgs_20 = c(imgs_20_a,imgs_20_b,imgs_20_c,imgs_20_d)
rm(imgs_20_a,imgs_20_b,imgs_20_c,imgs_20_d)
imgs_20 = lapply(imgs_20,function(x) x@.Data)
imgs_20 = sapply(imgs_20,as.vector)
imgs_20 = t(imgs_20)
dim(imgs_20)



# #: ------------------------------------------------------
# #:  Identify Very similar pictures
# #: ------------------------------------------------------
# # -> excluding duplicated pictures for later density estimation
# imgs_20_dist=dist(imgs_20,method="euclidean")
# windows()
# hist(imgs_20_dist)
# imgs_20_distMatrix=as.matrix(imgs_20_dist)
# ident_pics=which(imgs_20_distMatrix<=0.71,arr.in=T)  #sqrt(0.05^2*200)=0.7
# ident_pics=ident_pics[ident_pics[,1]!=ident_pics[,2],]
# id=paste(rowSums(ident_pics),apply(ident_pics,1,max))
# ident_pics=ident_pics[!duplicated(id),]
# all(ident_pics[,1]>ident_pics[,2])
# length(ident_pics)
# rm_pic_dens=split(ident_pics,f=ident_pics[,2])
# length(unique(as.vector(ident_pics)))-length(rm_pic_dens) # number of removed pic out of 3204


#: ------------------------------------------------------
#:  PCA for dimensioal reduction, 20 by 20
#: ------------------------------------------------------
pca_initial=prcomp(imgs_20)
pca_initial_check=summary(pca_initial)$importance
cumsum(summary(pca_initial)$importance[2,])

pc_imgs = list()
for(pc_id in 1:30){
  temp_pc_img=1-(pca_initial$rotation[,pc_id]-min(pca_initial$rotation[,pc_id]))/(max(pca_initial$rotation[,pc_id])-min(pca_initial$rotation[,pc_id]))
  temp_pc_img=Image(matrix(temp_pc_img,20,20),dim=c(20,20),colormode="Grayscale")
  pc_imgs=c(pc_imgs,list(temp_pc_img))
}
windows()
par(mfrow=c(6,5))
for(i in 1:30) plot(pc_imgs[[i]])
for(i in 1:20) writeImage(pc_imgs[[i]],files=paste0("h:/data/Projects/Research/MFI/Documents/Final Figures/pc20/pc",i,".png"),type="png")

eigvec_20x20=pca_initial$rotation[,1:20]
pc_scores= imgs_20%*%eigvec_20x20
app_pics= pc_scores%*%t(eigvec_20x20)

app_pics=sweep(app_pics,1,rowMeans(imgs_20)-rowMeans(app_pics),FUN="+")
sst_values = apply((imgs_20-app_pics)^2,1,sum)
# windows()
# hist(sst_values)
summary(sst_values)
quantile(sst_values,0.995)

check_id=which(sst_values>quantile(sst_values,0.95))
windows()
par(mfrow=c(20,10))
for(i in check_id) plot(as.Image(matrix(imgs_20[i,],20,20)))


windows()
windows()
par(mfrow=c(20,10))
for(i in 1:200) plot(as.Image(matrix(app_pics[i,],20,20)))

### check the worst approximated pic
###   Checking the PC approximation, seems good
# which.min(sst_values)
# org_pic=imgs_20[281,]
# app_pic=app_pics[281,]
# #app_pic=app_pic+mean(org_pic)-mean(app_pic)
# sum((org_pic-app_pic)^2)
# display(as.Image(matrix(org_pic,20,20)))
# display(as.Image(matrix(app_pic,20,20)))
# writeImage(as.Image(matrix(org_pic,20,20)),files="h:/data/Projects/Research/MFI/Results/Publication/approx img/SO_raw_good.png",type="png")
# writeImage(as.Image(matrix(app_pic,20,20)),files="h:/data/Projects/Research/MFI/Results/Publication/approx img/SO_approx_good.png",type="png")


# -> Data Depth based on first 30 loading
use_dat=pc_scores
# use_dat=cbind(AvgLight=use_dat,rowMeans(imgs_20)-rowMeans(use_dat%*%t(pc_scores[,1:10])))

# windows()
# cl.res=fviz_nbclust(use_dat, kmeans, method = "silhouette",k.max=50)
# # str(cl.res)
# # cl.res$data
# 
# windows()
# par(mar=c(5,4,2,1),tcl=-0.15,mgp=c(2.5,0.5,0),font.lab=6,font.axis=6,cex.axis=0.9)
# plot(cl.res$data$y,pch=16,cex=0.8,col="skyblue",ylab="Average Sihouetter Metric",xlab="Number of Clusters",main="",
#      xxat="n",yaxt="n")
# lines(cl.res$data$y,lty=1,col="skyblue")
# axis(2,las=1)
# abline(v=40,col="gray",lty=2)

set.seed(12345)
kid=kmeans(use_dat,centers=30)$cluster
range(summary(as.factor(kid)))

gr_imgs_raw = lapply(1:max(kid),function(x) imgs_20[kid==x,])
gr_imgs = lapply(1:max(kid),function(x) use_dat[kid==x,])
gr_imgs_avg = lapply(gr_imgs_raw,function(x){
  values=apply(x,2,mean,trim=0.05)
  values=(values-min(values))/(max(values)-min(values))
  as.Image(matrix(values,20,20))
})

windows()
par(mfrow=c(5,6))
for(i in seq(gr_imgs_avg)) plot(gr_imgs_avg[[i]])

# id=3
# set.seed(202003)
# simlate_pics = rmvnorm(100,colMeans(gr_imgs[[id]]),sigma=cov(gr_imgs[[id]]))
# simlate_pics = lapply(1:nrow(simlate_pics),function(i){
#   values=simlate_pics[i,]%*%t(pc_scores)
#   values=(values-min(values))/(max(values)-min(values))
#   as.Image(matrix(values,20,20))
# })
# windows()
# par(mfrow=c(10,10))
# for(i in seq(simlate_pics)) plot(simlate_pics[[i]])

set.seed(2020)
dd_res=sapply(seq(gr_imgs), function(i){
  org_dd = depth.Mahalanobis(gr_imgs[[i]],gr_imgs[[i]])
  org_dd_cut= min(org_dd)
  new_pics = rmvnorm(5000,colMeans(gr_imgs[[i]]),sigma=cov(gr_imgs[[i]]))
  sim_dd = depth.Mahalanobis(new_pics,data=gr_imgs[[i]])
  sim_dd_99 = quantile(sim_dd,0.01)
  sim_dd_995 = quantile(sim_dd,0.005)
  org_dd_cut_p = sum(sim_dd<round(org_dd_cut,4))/5000
  nr_org_f99 = sum(org_dd<round(sim_dd_99,4))
  nr_org_f995 = sum(org_dd<round(sim_dd_995,4))
  c(org_dd_cut=org_dd_cut,
    org_dd_cut_p=org_dd_cut_p,
    sim_dd_99=sim_dd_99,
    sim_dd_995=sim_dd_995,
    nr_org_f99=nr_org_f99,
    nr_org_f995=nr_org_f995)
})
dd_res=t(dd_res)
range(dd_res[,1])
range(dd_res[,4])
mean(dd_res[,1])
median(dd_res[,1])

train_test_res=sapply(seq(gr_imgs), function(i){
  org_dd = depth.Mahalanobis(use_dat,gr_imgs[[i]])
  org_dd
})
summary(apply(train_test_res,1,max))
sum(apply(train_test_res,1,max)<0.015)

check=cbind(rownames(imgs_20),ifelse(apply(train_test_res,1,sum)==0,"Fail","Pass"))
check=check[check[,2]=="Fail",]

windows()
par(mfrow=c(8,10))
for(i in which(apply(train_test_res,1,max)<0.015)) plot(as.Image(matrix(imgs_20[i,],20,20)))

#: -> Make reference pictures 20 by 20
set.seed(12345)
kid_50=kmeans(use_dat,centers=30)$cluster
summary(as.factor(kid_50))
gr_imgs_raw_50 = lapply(1:max(kid_50),function(x) imgs_20[kid_50==x,])
gr_imgs_avg_50 = lapply(gr_imgs_raw_50,function(x){
  values=apply(x,2,mean,trim=0.05)
  values=(values-min(values))/(max(values)-min(values))
  as.Image(matrix(values,20,20))
})
windows()
par(mfrow=c(5,10))
for(i in seq(gr_imgs_avg_50)) plot(gr_imgs_avg_50[[i]])
reference=t(sapply(gr_imgs_avg_50,as.vector))
dif_fil_50=apply(imgs_20,1,function(x){
  temp=sweep(reference,2,x,FUN="-")
  temp=sweep(temp,1,rowMeans(temp),FUN="-")
  min(apply(temp,1,function(y) sum(abs(y))))
})
summary(dif_fil_50)
windows()
plot(dif_fil_50)
quantile(dif_fil_50,probs=0.99)


reference=rbind(t(sapply(gr_imgs_avg_50,as.vector)),imgs_20[dif_fil_50>38,])
windows()
par(mfrow=c(15,10))
for(i in 1:nrow(reference)) plot(as.Image(matrix(reference[i,],20,20)))
#reference=reference[-c(55,58,59,62,63,65,67,73,74,77:80,82,85),]
reference=reference[-c(32,38,40,44:49,51:54,56,59:63,67:70,72,74:77,80:81,83,85:88,90,91,93:96,98:103,105:110,112,116,119:121,124:126,129:140,142),]
windows()
par(mfrow=c(10,7))
for(i in 1:nrow(reference)) plot(as.Image(matrix(reference[i,],20,20)))
dif_fil_final=apply(imgs_20,1,function(x){
  temp=sweep(reference,2,x,FUN="-")
  temp=sweep(temp,1,rowMeans(temp),FUN="-")
  min(apply(temp,1,function(y) sum(abs(y))))
})
summary(dif_fil_final)
windows()
plot(dif_fil_final)
quantile(dif_fil_final,probs=0.99)

#: ------------------------------------------------------
#:  Small picture Index
#: ------------------------------------------------------
raw_size$small = pmax(raw_size$nrow,raw_size$ncol)<=9
sum(raw_size$small)

raw_small_imgs=raw_imgs[names(raw_imgs)%in%(raw_size$ID[raw_size$small])]

imgs_9_a = lapply(raw_small_imgs,function(temp){
  temp=(temp-min(temp))/(max(temp-min(temp)))
  temp=chop_blank(temp)
  temp=(temp-min(temp))/(max(temp-min(temp)))
  temp=chop_blank(temp)
  temp=(temp-min(temp))/(max(temp-min(temp)))
  if(abs(nrow(temp)-ncol(temp))>1){
    if(ncol(temp)<nrow(temp)){
      pad_val=(max(temp[,1])+max(temp[,ncol(temp)]))/2
      nr_pad=nrow(temp)-ncol(temp)
      pad_mat=matrix(pad_val,nrow=nrow(temp),ncol=floor(nr_pad/2))
      temp=cbind(pad_mat,temp,pad_mat)
    }else{
      pad_val=(max(temp[1,])+max(temp[nrow(temp),]))/2
      nr_pad=ncol(temp)-nrow(temp)
      pad_mat=matrix(pad_val,nrow=floor(nr_pad/2),ncol=ncol(temp))
      temp=rbind(pad_mat,temp,pad_mat)
    }
  }
  resize(temp,9,9)
})
imgs_9_b = lapply(imgs_9_a,function(temp) rotate(temp,90))
names(imgs_9_b)=paste0(names(imgs_9_a),"_R90")
imgs_9_c = lapply(imgs_9_a,function(temp) rotate(temp,180))
names(imgs_9_c)=paste0(names(imgs_9_a),"_R180")
imgs_9_d = lapply(imgs_9_a,function(temp) rotate(temp,270))
names(imgs_9_d)=paste0(names(imgs_9_a),"_R270")

imgs_9 = c(imgs_9_a,imgs_9_b,imgs_9_c,imgs_9_d)
rm(imgs_9_a,imgs_9_b,imgs_9_c,imgs_9_d)
imgs_9 = lapply(imgs_9,function(x) x@.Data)
imgs_9 = sapply(imgs_9,as.vector)
imgs_9 = t(imgs_9)
dim(imgs_9)

pca_9=prcomp(imgs_9)
(pca_9_check=summary(pca_9)$importance)

pc_imgs_9 = list()
for(pc_id in 1:22){
  temp_pc_img=1-(pca_9$rotation[,pc_id]-min(pca_9$rotation[,pc_id]))/(max(pca_9$rotation[,pc_id])-min(pca_9$rotation[,pc_id]))
  temp_pc_img=Image(matrix(temp_pc_img,9,9),dim=c(9,9),colormode="Grayscale")
  pc_imgs_9=c(pc_imgs_9,list(temp_pc_img))
}
windows()
par(mfrow=c(5,5))
for(i in 1:22) plot(pc_imgs_9[[i]])

eigvec_9x9=pca_9$rotation[,1:22]
pc_sc_9= imgs_9%*%eigvec_9x9
app_pics_9= pc_sc_9%*%t(eigvec_9x9)

sst_values_9 = apply((imgs_9-app_pics_9)^2,1,sum)
# windows()
# hist(sst_values)
summary(sst_values_9)
quantile(sst_values_9,0.995)


use_dat_9=pc_sc_9

set.seed(2020)
kid_9=kmeans(use_dat_9,centers=30)$cluster
summary(as.factor(kid_9))

gr_imgs_raw_9 = lapply(1:max(kid_9),function(x) imgs_9[kid_9==x,])
gr_imgs_9 = lapply(1:max(kid_9),function(x) use_dat_9[kid_9==x,])
# gr_imgs_avg_9 = lapply(seq(gr_imgs_raw_9),function(i){
#   x=gr_imgs_9[[i]]
#   dd = depth.Mahalanobis(x,x)
#   values=gr_imgs_raw_9[[i]][which.max(dd),]
#   as.Image(matrix(values,9,9))
# })
gr_imgs_avg_9 = lapply(gr_imgs_raw_9,function(x){
  values=apply(x,2,mean,trim=0.1)
  values=(values-min(values))/(max(values)-min(values))
  as.Image(matrix(values,9,9))
})


windows()
par(mfrow=c(5,6))
for(i in seq(gr_imgs_avg_9)) plot(gr_imgs_avg_9[[i]])

reference_9=t(sapply(gr_imgs_avg_9,as.vector))
dif_fil_9=apply(imgs_9,1,function(x){
  temp=sweep(reference_9,2,x,FUN="-")
  temp=sweep(temp,1,rowMeans(temp),FUN="-")
  min(apply(temp,1,function(y) sum(abs(y))))
})
summary(dif_fil_9)
windows()
plot(dif_fil_9)
quantile(dif_fil_9,0.95)

reference_9=rbind(t(sapply(gr_imgs_avg_9,as.vector)),imgs_9[dif_fil_9>6.7,])
windows()
par(mfrow=c(13,10))
for(i in 1:nrow(reference_9)) plot(as.Image(matrix(reference_9[i,],9,9)))
#reference_9=reference_9[c(1:7,9:30,51,58,74,78,85,98,104,119),]
#reference_9=reference_9[c(1,2,5:7,9:14,17,21:28,30,51,58,74,78,85,98,104,119),]
reference_9=reference_9[c(1:13,15:24,26:31,52,75,78,123),]
windows()
par(mfrow=c(7,5))
for(i in 1:nrow(reference_9)) plot(as.Image(matrix(reference_9[i,],9,9)))
dif_fil_9=apply(imgs_9,1,function(x){
  temp=sweep(reference_9,2,x,FUN="-")
  temp=sweep(temp,1,rowMeans(temp),FUN="-")
  min(apply(temp,1,function(y) sum(abs(y))))
})
summary(dif_fil_9)
windows()
plot(dif_fil_9)
quantile(dif_fil_9,probs=0.90)

# 
# #-----------------------------------------------------
# test_data_path="h:/data/Projects/Research/MFI/Back Up Data/Real Examples/filter test/"  # testing,image and metadata
# pool=c(2,3,4)
# raw_test_pic=list()
# for(j in pool){
#   temp_files=list.files(paste0(test_data_path,j,"/check Oil/"))
#   temp=lapply(temp_files,function(y)readImage(file=paste0(test_data_path,j,"/check Oil/",y)))
#   raw_test_pic=c(raw_test_pic,temp)
# }
# length(raw_test_pic)
# small_ind= sapply(raw_test_pic,function(y) max(dim(y)))<=9
# raw_test_pic_vec=sapply(raw_test_pic,function(temp) {
#   temp=(temp-min(temp))/(max(temp-min(temp)))
#   temp=chop_blank(temp)
#   temp=(temp-min(temp))/(max(temp-min(temp)))
#   temp=chop_blank(temp)
#   temp=(temp-min(temp))/(max(temp-min(temp)))
#   if(abs(nrow(temp)-ncol(temp))>1){
#     if(ncol(temp)<nrow(temp)){
#       pad_val=(max(temp[,1])+max(temp[,ncol(temp)]))/2
#       nr_pad=nrow(temp)-ncol(temp)
#       pad_mat=matrix(pad_val,nrow=nrow(temp),ncol=floor(nr_pad/2))
#       temp=cbind(pad_mat,temp,pad_mat)
#     }else{
#       pad_val=(max(temp[1,])+max(temp[nrow(temp),]))/2
#       nr_pad=ncol(temp)-nrow(temp)
#       pad_mat=matrix(pad_val,nrow=floor(nr_pad/2),ncol=ncol(temp))
#       temp=rbind(pad_mat,temp,pad_mat)
#     }
#   }
#   temp=resize(temp,20,20)
#   as.vector(temp@.Data)
# })
# raw_test_pic_vec_9=sapply(raw_test_pic,function(temp) {
#   temp=(temp-min(temp))/(max(temp-min(temp)))
#   temp=chop_blank(temp)
#   temp=(temp-min(temp))/(max(temp-min(temp)))
#   temp=chop_blank(temp)
#   temp=(temp-min(temp))/(max(temp-min(temp)))
#   if(abs(nrow(temp)-ncol(temp))>1){
#     if(ncol(temp)<nrow(temp)){
#       pad_val=(max(temp[,1])+max(temp[,ncol(temp)]))/2
#       nr_pad=nrow(temp)-ncol(temp)
#       pad_mat=matrix(pad_val,nrow=nrow(temp),ncol=floor(nr_pad/2))
#       temp=cbind(pad_mat,temp,pad_mat)
#     }else{
#       pad_val=(max(temp[1,])+max(temp[nrow(temp),]))/2
#       nr_pad=ncol(temp)-nrow(temp)
#       pad_mat=matrix(pad_val,nrow=floor(nr_pad/2),ncol=ncol(temp))
#       temp=rbind(pad_mat,temp,pad_mat)
#     }
#   }
#   temp=resize(temp,9,9)
#   as.vector(temp@.Data)
# })
# raw_test_pic_vec=t(raw_test_pic_vec)
# raw_test_pic_vec_9=t(raw_test_pic_vec_9)
# 
# loading_test=solve(t(pc_scores)%*%pc_scores)%*%t(pc_scores)%*%t(raw_test_pic_vec)
# loading_test=t(loading_test)
# # loading_test_9=solve(t(pc_sc_9)%*%pc_sc_9)%*%t(pc_sc_9)%*%t(raw_test_pic_vec_9)
# # loading_test_9=t(loading_test_9)
# 
# app_pics_test = loading_test%*%t(pc_scores)
# app_pics_test = sweep(app_pics_test,1,rowMeans(raw_test_pic_vec)-rowMeans(app_pics_test),FUN="+")
# sst_values_test = apply((raw_test_pic_vec-app_pics_test)^2,1,sum)
# sum(sst_values_test>1.1)
# # app_pics_test_9 = loading_test_9%*%t(pc_sc_9)
# # app_pics_test_9 = sweep(app_pics_test_9,1,rowMeans(raw_test_pic_vec_9)-rowMeans(app_pics_test_9),FUN="+")
# # sst_values_test_9 = apply((raw_test_pic_vec_9-app_pics_test_9)^2,1,sum)
# # sum(sst_values_test_9>0.18)
# # sum(sst_values_test_9[small_ind]>0.18)
# 
# test_IA_res=sapply(seq(gr_imgs), function(i){
#   org_dd = depth.Mahalanobis(loading_test,gr_imgs[[i]])
#   org_dd 
# })
# # test_IA_res_9=sapply(seq(gr_imgs_9), function(i){
# #   org_dd = depth.Mahalanobis(loading_test_9,gr_imgs_9[[i]])
# #   org_dd 
# # })
# 
# mean(sst_values_test<1.1)
# layer1 = sst_values_test<1.1 & apply(test_IA_res,1,max)>0.01
# mean(layer1)
# # layer2 = layer1
# # layer2[(!layer1)&small_ind] = (apply(test_IA_res_9,1,max)>0.02)[(!layer1)&small_ind]
# # mean(layer2)
# 
# test_dif_fil=apply(raw_test_pic_vec,1,function(x){
#   temp=sweep(reference,2,x,FUN="-")
#   temp=sweep(temp,1,rowMeans(temp),FUN="-")
#   min(apply(temp,1,function(y) sum(abs(y))))
# })
# test_dif_fil_9=apply(raw_test_pic_vec_9,1,function(x){
#   temp=sweep(reference_9,2,x,FUN="-")
#   temp=sweep(temp,1,rowMeans(temp),FUN="-")
#   min(apply(temp,1,function(y) sum(abs(y))))
# })
# layer2 = test_dif_fil<38
# layer3 = test_dif_fil_9<7
# pic_test=layer2
# pic_test[small_ind]=layer3[small_ind]
# 
# mean(layer1)
# mean(layer2)
# mean(layer3)
# 
# (sum(layer1) + sum(test_dif_fil[!layer1]<38))/length(sst_values_test)
# (sum(layer1) + sum(pic_test[!layer1]))/length(sst_values_test)
# 
# #sum((test_dif_fil_20<=45)|(apply(test_IA_res_oil,1,sum)>0))
# check=raw_test_pic[layer1]
# 
# windows()
# par(mfrow=c(10,10))
# for(i in 1:length(check)) plot(check[[i]])
# 
# 
# #save(list=ls(all=T),file="H:/data/Projects/Research/MFI/Programs/IA_Train.RData")
# # ---------------------------------------------------------------
# #  Appendix: Exercise with individuals to tune some para
# # ---------------------------------------------------------------
# 
# # test image: 163534.jpg, 28155.jpg, 30223.jpg, 329235.jpg, 274720.jpg, 23651.jpg(good one for backgroud)
# temp=readImage(file=paste0(IA_train_path,"23651.jpg"))
# dim(temp)
# display(temp)
# dim(chop_blank(temp))
# display(chop_blank(temp))
# display(resize(chop_blank(temp),30,30))
# 
# # contrast preprocessing, 40109.jpg 
# temp=readImage(file=paste0(IA_train_path,"40109.jpg"))
# temp=(temp-min(temp))/(max(temp-min(temp)))
# dim(temp)
# display(temp)
# dim(chop_blank(temp))
# display(chop_blank(temp))
# display(resize(chop_blank(temp),30,30))
# 
# # contrast preprocessing 
# temp=readImage(file=paste0(IA_train_path,"195212.jpg"))
# temp=(temp-min(temp))/(max(temp-min(temp)))
# temp=chop_blank(temp)
# temp=(temp-min(temp))/(max(temp-min(temp)))
# temp=chop_blank(temp)
# temp=(temp-min(temp))/(max(temp-min(temp)))
# display(resize(temp,20,20))
# 
# # contrast preprocessing 
# temp=readImage(file=paste0(IA_train_path,"31942.jpg"))
# temp=(temp-min(temp))/(max(temp-min(temp)))
# temp=chop_blank(temp)
# temp=(temp-min(temp))/(max(temp-min(temp)))
# temp=chop_blank(temp)
# temp=(temp-min(temp))/(max(temp-min(temp)))
# display(resize(temp,20,20))
# 
# # contrast preprocessing 
# temp=readImage(file=paste0(IA_train_path,"333869.jpg"))
# temp=(temp-min(temp))/(max(temp-min(temp)))
# temp=chop_blank(temp)
# temp=(temp-min(temp))/(max(temp-min(temp)))
# temp=chop_blank(temp)
# temp=(temp-min(temp))/(max(temp-min(temp)))
# temp2=resize(temp,20,20)
# display(temp2)
# display(rotate(temp2,180))
# 












