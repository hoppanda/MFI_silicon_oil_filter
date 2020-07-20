
library(EBImage)
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


# Oil pic
IA_train_path = "H:/data/Projects/Research/MFI/Data/Filter Picture/"
temp=readImage(file=paste0(IA_train_path,"8893.jpg"))
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/SO_0_Raw.png",type="png")
temp=(temp-min(temp))/(max(temp-min(temp)))
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/SO_1.png",type="png")
temp=chop_blank(temp)
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/SO_2.png",type="png")
temp=(temp-min(temp))/(max(temp-min(temp)))
temp=chop_blank(temp)
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/SO_3.png",type="png")
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
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/SO_4.png",type="png")
temp=resize(temp,20,20)
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/SO_5.png",type="png")




IA_train_path = "H:/data/Projects/Research/MFI/Back Up Data/Real Examples/filter test/2/check non oil/"
temp=readImage(file=paste0(IA_train_path,"4326.jpg"))
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/NSO_0_Raw.png",type="png")
temp=(temp-min(temp))/(max(temp-min(temp)))
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/NSO_1.png",type="png")
temp=chop_blank(temp)
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/NSO_2.png",type="png")
temp=(temp-min(temp))/(max(temp-min(temp)))
temp=chop_blank(temp)
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/NSO_3.png",type="png")
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
temp=as.Image(temp)
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/NSO_4.png",type="png")
temp=resize(temp,20,20)
writeImage(temp,files="h:/data/Projects/Research/MFI/Results/Publication/preprocess/NSO_5.png",type="png")






















