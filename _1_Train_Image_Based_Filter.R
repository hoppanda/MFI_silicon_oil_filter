#############################################################################
#
#      Reproduce the results from Paper - Train IBF
# 
#
#   ------------------------------------------------------------------------
###    Correspondence: X. Gregory, Chen 
###                    stat1013@gmail.com
###
###    Version: Oct 2020
###
#############################################################################

#:=>  Dependencies: (+ Requireed R package, tested version number)
#       + EBImage, 4.22.1
#       + mvtnorm, 1.1-0
#       + cluster, 2.0.7-1
#       + factorextra, 1.0.6
#       + ddalpha, 1.3.11

#: ------------------------------------------------------
#:  Pre-setup
#: ------------------------------------------------------

# *** To install EBImage package for the first time
# source("http://bioconductor.org/biocLite.R")
# biocLite("EBImage")

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
  row.ind= range(which(row.cut))
  row.cut[row.ind[1]:row.ind[2]] = TRUE
  col.cut=!(ch2a>0.78 & ch2b<0.28)
  col.ind= range(which(col.cut))
  col.cut[col.ind[1]:col.ind[2]] = TRUE
  tpic[row.cut,col.cut,drop=F]
}

#: ------------------------------------------------------
#:  Loading Data, Augmentation, Preprocess
#: ------------------------------------------------------
setwd("h:/data/Projects/Research/MFI/Programs/GitHub/RData Files/") # change as per your case
load("Raw_Imgs_and_Meta_Data.RData") # download this from "RData Files" in Github

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

#: ------------------------------------------------------
#:  PCA for dimensioal reduction, 20 by 20, Layer 1
#: ------------------------------------------------------
pca_initial=prcomp(imgs_20)
pca_initial_check=summary(pca_initial)$importance

pc_imgs = list()
for(pc_id in 1:30){
  temp_pc_img=1-(pca_initial$rotation[,pc_id]-min(pca_initial$rotation[,pc_id]))/(max(pca_initial$rotation[,pc_id])-min(pca_initial$rotation[,pc_id]))
  temp_pc_img=Image(matrix(temp_pc_img,20,20),dim=c(20,20),colormode="Grayscale")
  pc_imgs=c(pc_imgs,list(temp_pc_img))
}
windows()
par(mfrow=c(6,5))
for(i in 1:30) plot(pc_imgs[[i]])

# -> calculating dPCA
eigvec_20x20=pca_initial$rotation[,1:20]
pc_scores= imgs_20%*%eigvec_20x20
app_pics= pc_scores%*%t(eigvec_20x20)
app_pics=sweep(app_pics,1,rowMeans(imgs_20)-rowMeans(app_pics),FUN="+")
sst_values = apply((imgs_20-app_pics)^2,1,sum)
windows()
hist(sst_values)
summary(sst_values)
trunc(quantile(sst_values,0.995)*10)/10  # conservative cut-off for dPCA

check_id=which(sst_values>quantile(sst_values,0.95)) # check gallery of most deviating pic in training set
windows()
par(mfrow=c(20,10))
for(i in check_id) plot(as.Image(matrix(imgs_20[i,],20,20)))
windows()
par(mfrow=c(20,10))
for(i in 1:200) plot(as.Image(matrix(app_pics[i,],20,20)))

# -> Data Depth based on first 30 loading
use_dat=pc_scores
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
summary(dd_res[,4]) # range of the 99.5% cutoff of dMD across 30 clusters

train_test_res=sapply(seq(gr_imgs), function(i){
  org_dd = depth.Mahalanobis(use_dat,gr_imgs[[i]])
  org_dd
})
summary(apply(train_test_res,1,max))
quantile(apply(train_test_res,1,max),0.005)
sum(apply(train_test_res,1,max)<0.02)

# check the gallery of the most deviating pics
check=cbind(rownames(imgs_20),ifelse(apply(train_test_res,1,sum)==0,"Fail","Pass"))
check=check[check[,2]=="Fail",]
windows()
par(mfrow=c(8,10))
for(i in which(apply(train_test_res,1,max)<0.02)) plot(as.Image(matrix(imgs_20[i,],20,20)))

#: ------------------------------------------------------
#:  Reference book, 20 by 20, Layer 2
#: ------------------------------------------------------
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
quantile(dif_fil_50,probs=0.99) # interiam cutoff, not the final one

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
quantile(dif_fil_final,probs=0.99) # final cutoff

#: ------------------------------------------------------
#:  Reference book, 9 by 9, Layer 2
#: ------------------------------------------------------
raw_size$small = pmax(raw_size$nrow,raw_size$ncol)<=9
sum(raw_size$small)

# resize
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

# pca for average images
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
use_dat_9=pc_sc_9

set.seed(2020)
kid_9=kmeans(use_dat_9,centers=30)$cluster
summary(as.factor(kid_9))
gr_imgs_raw_9 = lapply(1:max(kid_9),function(x) imgs_9[kid_9==x,])
gr_imgs_9 = lapply(1:max(kid_9),function(x) use_dat_9[kid_9==x,])
gr_imgs_avg_9 = lapply(gr_imgs_raw_9,function(x){
  values=apply(x,2,mean,trim=0.1)
  values=(values-min(values))/(max(values)-min(values))
  as.Image(matrix(values,9,9))
})
windows()
par(mfrow=c(5,6))
for(i in seq(gr_imgs_avg_9)) plot(gr_imgs_avg_9[[i]])

# compile the reference book 
reference_9=t(sapply(gr_imgs_avg_9,as.vector))
dif_fil_9=apply(imgs_9,1,function(x){
  temp=sweep(reference_9,2,x,FUN="-")
  temp=sweep(temp,1,rowMeans(temp),FUN="-")
  min(apply(temp,1,function(y) sum(abs(y))))
})
summary(dif_fil_9)
windows()
plot(dif_fil_9)
quantile(dif_fil_9,0.95) # interim cutoff, not the final one

reference_9=rbind(t(sapply(gr_imgs_avg_9,as.vector)),imgs_9[dif_fil_9>6.7,])
windows()
par(mfrow=c(13,10))
for(i in 1:nrow(reference_9)) plot(as.Image(matrix(reference_9[i,],9,9)))
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
quantile(dif_fil_9,probs=0.90) # final cutoff
