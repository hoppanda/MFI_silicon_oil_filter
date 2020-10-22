#############################################################################
#
#      Reproduce the results from Paper - Evaluation on the Test Set
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
#       + randomForest, 4.6-14

#:-> Run this script right after "_1_Train_Image_Based_Filter.R", "_2_Train_RF.R"
#     DO NOT clear up the previous session image


#: ------------------------------------------------------
#:  Functions
#: ------------------------------------------------------
#: Benchmark 1 ~ Aspect Ration >=0.85
#: Benckmark 2
strehl=function(use){
  sf=use$Circularity*use$Aspect.Ratio*use$Intensity.Max*use$Intensity.STD
  th=sapply(use$ECD.um.,function(x){
    if(x<=2.875){29732*log(x)+3913
    }else{
      if(x<=6.125){46609*log(x)-10464
      }else{3368*x+49124}
    }
  })
  ifelse(sf>=th,"Oil","Not Oil")
}

#: ------------------------------------------------------
#:  Some discriptive of test set
#: ------------------------------------------------------
quantile(test_oil_meta$ECD.um.,probs=c(0.025,0.05,0.95,0.975))
quantile(test_badoil_meta$ECD.um.,probs=c(0.025,0.05,0.95,0.975))
quantile(test_notoil_meta$ECD.um.,probs=c(0.025,0.05,0.95,0.975))
nrow(test_oil_meta)
nrow(test_badoil_meta)
nrow(test_notoil_meta)

sum(apply(test_oil_size[,c("nrow","ncol")],1,max)<=9)
sum(apply(test_oil_size[,c("nrow","ncol")],1,max)>9)
sum(apply(test_badoil_size[,c("nrow","ncol")],1,max)<=9)
sum(apply(test_badoil_size[,c("nrow","ncol")],1,max)>9)
sum(apply(test_notoil_size[,c("nrow","ncol")],1,max)<=9)
sum(apply(test_notoil_size[,c("nrow","ncol")],1,max)>9)
check=t(sapply(test_oil_imgs,dim))
all(row.names(check)==row.names(test_oil_size))
sum(apply(check,1,max)<=9)


#: ------------------------------------------------------
#:  Testing Outcome: Benchmark 1 (AR>=0.89)
#: ------------------------------------------------------
test_ben1_oil_res = factor(ifelse(test_oil_meta$Aspect.Ratio>=0.85,"Oil","Not Oil"))
test_ben1_badoil_res = factor(ifelse(test_badoil_meta$Aspect.Ratio>=0.85,"Oil","Not Oil"))
test_ben1_notoil_res = factor(ifelse(test_notoil_meta$Aspect.Ratio>=0.85,"Oil","Not Oil"))

summary(test_ben1_oil_res)
summary(test_ben1_badoil_res)
summary(test_ben1_notoil_res)

all(rownames(test_oil_meta)[match(test_oil_size$ID,rownames(test_oil_meta))]==test_oil_size$ID)

test_oil_size$ben_1 = test_ben1_oil_res[match(test_oil_size$ID,rownames(test_oil_meta))]
test_badoil_size$ben_1 = test_ben1_badoil_res[match(test_badoil_size$ID,rownames(test_badoil_meta))]
test_notoil_size$ben_1 = test_ben1_notoil_res[match(test_notoil_size$ID,rownames(test_notoil_meta))]

#: ------------------------------------------------------
#:  Testing Outcome: Benchmark 2 (S factor, Strehl et.al.)
#: ------------------------------------------------------
test_ben2_oil_res = factor(strehl(test_oil_meta),levels=c("Not Oil","Oil"))
test_ben2_badoil_res = factor(strehl(test_badoil_meta),levels=c("Not Oil","Oil"))
test_ben2_notoil_res = factor(strehl(test_notoil_meta),levels=c("Not Oil","Oil"))

summary(test_ben2_oil_res)
summary(test_ben2_badoil_res)
summary(test_ben2_notoil_res)

test_oil_size$ben_2 = test_ben2_oil_res[match(test_oil_size$ID,rownames(test_oil_meta))]
test_badoil_size$ben_2 = test_ben2_badoil_res[match(test_badoil_size$ID,rownames(test_badoil_meta))]
test_notoil_size$ben_2 = test_ben2_notoil_res[match(test_notoil_size$ID,rownames(test_notoil_meta))]


#: ------------------------------------------------------
#:  Testing Outcome: RF 
#: ------------------------------------------------------
# Testing Outcome: Random Forest (RF in the main text)
use_names=names(rft_dat)[names(rft_dat)!="res"]
sum(names(test_oil_meta)%in%use_names)
length(use_names)
test_rf_oil_res=predict(mod_rf,newdata=test_oil_meta[,names(test_oil_meta)%in%use_names])
test_rf_badoil_res=predict(mod_rf,newdata=test_badoil_meta[,names(test_badoil_meta)%in%use_names])
test_rf_notoil_res=predict(mod_rf,newdata=test_notoil_meta[,names(test_notoil_meta)%in%use_names])

summary(test_rf_oil_res)
summary(test_rf_badoil_res)
summary(test_rf_notoil_res)

test_oil_size$rf = test_rf_oil_res[match(test_oil_size$ID,rownames(test_oil_meta))]
test_badoil_size$rf = test_rf_badoil_res[match(test_badoil_size$ID,rownames(test_badoil_meta))]
test_notoil_size$rf = test_rf_notoil_res[match(test_notoil_size$ID,rownames(test_notoil_meta))]

# Testing Outcome: Random Forest_0
use_names=names(rft_dat)[names(rft_dat)!="res"]
sum(names(test_oil_meta)%in%use_names)
length(use_names)
test_rf0_oil_res=predict(mod_rf0,newdata=test_oil_meta[,names(test_oil_meta)%in%use_names])
test_rf0_badoil_res=predict(mod_rf0,newdata=test_badoil_meta[,names(test_badoil_meta)%in%use_names])
test_rf0_notoil_res=predict(mod_rf0,newdata=test_notoil_meta[,names(test_notoil_meta)%in%use_names])

summary(test_rf0_oil_res)
summary(test_rf0_badoil_res)
summary(test_rf0_notoil_res)

test_oil_size$rf0 = test_rf0_oil_res[match(test_oil_size$ID,rownames(test_oil_meta))]
test_badoil_size$rf0 = test_rf0_badoil_res[match(test_badoil_size$ID,rownames(test_badoil_meta))]
test_notoil_size$rf0 = test_rf0_notoil_res[match(test_notoil_size$ID,rownames(test_notoil_meta))]


#: ------------------------------------------------------
#:  Testing Outcome: Image based Filter (IBF)
#: ------------------------------------------------------
dd_filter=function(raw_pics){
  small_ind= sapply(raw_pics,function(y) max(dim(y)))<=9
  raw_pic_vec=sapply(raw_pics,function(temp) {
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
    temp=resize(temp,20,20)
    as.vector(temp@.Data)
  })
  raw_pic_vec_9=sapply(raw_pics,function(temp) {
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
    temp=resize(temp,9,9)
    as.vector(temp@.Data)
  })
  raw_pic_vec=t(raw_pic_vec)
  raw_pic_vec_9=t(raw_pic_vec_9)
  #loading_test=solve(t(pc_scores)%*%pc_scores)%*%t(pc_scores)%*%t(raw_pic_vec)
  #loading_test=t(loading_test)
  loading_test= raw_pic_vec%*%eigvec_20x20
  
  app_pics_test = loading_test%*%t(eigvec_20x20)
  app_pics_test = sweep(app_pics_test,1,rowMeans(raw_pic_vec)-rowMeans(app_pics_test),FUN="+")
  sst_values_test = apply((raw_pic_vec-app_pics_test)^2,1,sum)
  
  # display(as.Image(matrix(app_pics_test[33,],20,20)))
  # display(as.Image(matrix(raw_pic_vec[33,],20,20)))
  # sst_values_test[33]
  # writeImage(as.Image(matrix(raw_pic_vec[33,],20,20)),files="h:/data/Projects/Research/MFI/Results/Publication/approx img/NSO_raw_worst.png",type="png")
  # writeImage(as.Image(matrix(app_pics_test[33,],20,20)),files="h:/data/Projects/Research/MFI/Results/Publication/approx img/NSO_approx_worst.png",type="png")
  
  test_IA_res=sapply(seq(gr_imgs), function(i){
    org_dd = depth.Mahalanobis(loading_test,gr_imgs[[i]])
    org_dd 
  })
  
  test_dif_fil_20=apply(raw_pic_vec,1,function(x){
    temp=sweep(reference,2,x,FUN="-")
    temp=sweep(temp,1,rowMeans(temp),FUN="-")
    min(apply(temp,1,function(y) sum(abs(y))))
  })
  test_dif_fil_9=apply(raw_pic_vec_9,1,function(x){
    temp=sweep(reference_9,2,x,FUN="-")
    temp=sweep(temp,1,rowMeans(temp),FUN="-")
    min(apply(temp,1,function(y) sum(abs(y))))
  })
  
  layer1 = sst_values_test<1.4 & apply(test_IA_res,1,max)>0.02
  layer2 = test_dif_fil_20<38
  layer3 = test_dif_fil_9<6
  pic_test=layer2
  pic_test[small_ind]=layer3[small_ind]

  output=data.frame(ID=names(raw_pics),sst=sst_values_test,
                    dd_max=apply(test_IA_res,1,max),
                    dd_layer1=ifelse(layer1,"Oil","Not Oil"),
                    dd_layer2=ifelse(pic_test,"Oil","Not Oil"),
                    dd_full=ifelse(layer1,"Oil",ifelse(pic_test,"Oil","Not Oil")))
  rownames(output)=NULL
  return(output)
}

test_dd_oil = dd_filter(raw_pics=test_oil_imgs)
test_dd_badoil = dd_filter(raw_pics=test_badoil_imgs)
test_dd_notoil = dd_filter(raw_pics=test_notoil_imgs)

test_oil_size=merge(test_oil_size,test_dd_oil,by="ID")
test_badoil_size=merge(test_badoil_size,test_dd_badoil,by="ID")
test_notoil_size=merge(test_notoil_size,test_dd_notoil,by="ID")

test_oil_meta$ID=row.names(test_oil_meta)
test_badoil_meta$ID=row.names(test_badoil_meta)
test_notoil_meta$ID=row.names(test_notoil_meta)

all(row.names(test_notoil_size)==row.names(test_notoil_meta))


#: ------------------------------------------------------
#:  Summarize the evaluation outcome
#: ------------------------------------------------------
reportbox=function(xx){
  xx$small= pmax(xx$nrow,xx$ncol)<=9
  makeres=
  data.frame(Category=c("Small","Not-Small","Total"),
             Number=c(sum(xx$small),sum(!xx$small),nrow(xx)),
             Ben_1=c(sum(xx$ben_1=="Not Oil" & xx$small),sum(xx$ben_1=="Not Oil" & (!xx$small)), sum(xx$ben_1=="Not Oil")),
             Ben_2=c(sum(xx$ben_2=="Not Oil" & xx$small),sum(xx$ben_2=="Not Oil" & (!xx$small)), sum(xx$ben_2=="Not Oil")),
             RF=c(sum(xx$rf=="Not Oil" & xx$small),sum(xx$rf=="Not Oil" & (!xx$small)), sum(xx$rf=="Not Oil")),
             RF0=c(sum(xx$rf0=="Not Oil" & xx$small),sum(xx$rf0=="Not Oil" & (!xx$small)), sum(xx$rf0=="Not Oil")),
             DD_1=c(sum(xx$dd_layer1=="Not Oil" & xx$small),sum(xx$dd_layer1=="Not Oil" & (!xx$small)), sum(xx$dd_layer1=="Not Oil")),
             DD_2=c(sum(xx$dd_layer2=="Not Oil" & xx$small),sum(xx$dd_layer2=="Not Oil" & (!xx$small)), sum(xx$dd_layer2=="Not Oil")),
             DD_full=c(sum(xx$dd_full=="Not Oil" & xx$small),sum(xx$dd_full=="Not Oil" & (!xx$small)), sum(xx$dd_full=="Not Oil"))
             )
  for(i in 3:9) makeres[,i]=round(makeres[,i]/makeres[,2],3)
  makeres
}
reportbox2=function(xx){
  xx$small= pmax(xx$nrow,xx$ncol)<=9
  makeres=
    data.frame(Category=c("Small","Not-Small","Total"),
               Number=c(sum(xx$small),sum(!xx$small),nrow(xx)),
               Ben_1=c(sum(xx$ben_1=="Oil" & xx$small),sum(xx$ben_1=="Oil" & (!xx$small)), sum(xx$ben_1=="Oil")),
               Ben_2=c(sum(xx$ben_2=="Oil" & xx$small),sum(xx$ben_2=="Oil" & (!xx$small)), sum(xx$ben_2=="Oil")),
               RF=c(sum(xx$rf=="Oil" & xx$small),sum(xx$rf=="Oil" & (!xx$small)), sum(xx$rf=="Oil")),
               RF0=c(sum(xx$rf0=="Oil" & xx$small),sum(xx$rf0=="Oil" & (!xx$small)), sum(xx$rf0=="Oil")),
               DD_1=c(sum(xx$dd_layer1=="Oil" & xx$small),sum(xx$dd_layer1=="Oil" & (!xx$small)), sum(xx$dd_layer1=="Oil")),
               DD_2=c(sum(xx$dd_layer2=="Oil" & xx$small),sum(xx$dd_layer2=="Oil" & (!xx$small)), sum(xx$dd_layer2=="Oil")),
               DD_full=c(sum(xx$dd_full=="Oil" & xx$small),sum(xx$dd_full=="Oil" & (!xx$small)), sum(xx$dd_full=="Oil"))
    )
  for(i in 3:9) makeres[,i]=round(makeres[,i]/makeres[,2],3)
  makeres
}

reportbox(test_oil_size)
reportbox(test_badoil_size)
reportbox2(test_notoil_size)

