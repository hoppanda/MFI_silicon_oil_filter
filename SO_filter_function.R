#############################################################################
#
#      Silicon Oil(SO) Image filter 
# 
#
#   ------------------------------------------------------------------------
###    Correspondence: X. Gregory, Chen 
###                    stat1013@gmail.com
###
###    Version: Oct 2020
###
#############################################################################

#:=>  Description:
#       This R file wraps up the resulting SO filter in a function <dd_filter> that can be readily applied.
#       Required R-objects to run this function is in "SO_Filter.RData"
#       One should first set the working directory to the local direcotry where "SO_Filter.RData" is stored 
#       for example, setwd("C:/myPC/Folder_where_RDATA_is_in")

#:=>  Dependencies: (+ Requireed R package, tested version number)
#       + EBImage, 4.22.1
#       + ddalpha, 1.3.11

setwd("h:/data/Projects/Research/MFI/Programs")

load("SO_Filter.RData")
#load("Raw_Imgs_and_Meta_Data.RData")
library(EBImage)
library(ddalpha)

#' Function: Clopping
#' -------------------------------------
#' @param tpics a matrix or a single 'Image' object (definded by EBImage) 
#' @return a matrix or a 'Image' object, where some row and columns may be chopped 
#'          if the average is >0.78 and Max-Min is <0.28 within that row/column
#' @seealso N/A
#' @export
#' 
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

#' Function: SO filter
#' -------------------------------------
#' @param raw_pics a list of or a single 'Image' object (definded by EBImage) 
#' @return a data frame with five columns, 
#'         ''
#'         'dd_max' is the maximum depth across the 30 clusters in layer 1, 
#'         'dd_layer1', 'dd_layer2','dd_full' are class prediction 'Oil'/'Not Oil' 
#'          from layer 1, layer 2 and all-combined respectively.
#'          One can use <output_data_frame>$dd_full to get access to the binary predicted label from the full filter.
#' @seealso N/A
#' @export
#' 
dd_filter=function(raw_pics){
  #prelimiary check
  if(!is.list(raw_pics)) raw_pics=list('Single Pic'=raw_pics)
  if(!is.Image(raw_pics[[1]])) stop("!!!! Input 'raw_pics' should be a list of or a single 'Image' object  !!!! ")
  #preprocessing
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
  
  #depth
  loading_test= raw_pic_vec%*%eigvec_20x20
  app_pics_test = loading_test%*%t(eigvec_20x20)
  app_pics_test = sweep(app_pics_test,1,rowMeans(raw_pic_vec)-rowMeans(app_pics_test),FUN="+")
  sst_values_test = apply((raw_pic_vec-app_pics_test)^2,1,sum)

  test_IA_res=sapply(seq(gr_imgs), function(i){
    org_dd = depth.Mahalanobis(loading_test,gr_imgs[[i]])
    org_dd 
  })
  if(length(raw_pics)==1) test_IA_res=matrix(test_IA_res,nrow=1)
  
  #reference book
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
  
  #filter layers
  layer1 = sst_values_test<1.4 & apply(test_IA_res,1,max)>0.02
  layer2 = test_dif_fil_20<38
  layer3 = test_dif_fil_9<6
  pic_test=layer2
  pic_test[small_ind]=layer3[small_ind]
  
  #outputing
  output=data.frame(ID=names(raw_pics),sst=sst_values_test,
                    dd_max=apply(test_IA_res,1,max),
                    dd_layer1=ifelse(layer1,"Oil","Not Oil"),
                    dd_layer2=ifelse(pic_test,"Oil","Not Oil"),
                    dd_full=ifelse(layer1,"Oil",ifelse(pic_test,"Oil","Not Oil")))
  rownames(output)=NULL
  output
}


# # Example
# dd_filter(test_oil_imgs[1:10])
# dd_filter(test_oil_imgs[[1]])






















