#############################################################################
#
#      Reproduce the results from Paper - Train RF Filter
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

#:-> Run this script right after "_1_Train_Image_Based_Filter.R", 
#     DO NOT clear up the previous session image

summary(rft_dat)
library(randomForest)

set.seed(12345)
cad_treesize=seq(3,30,by=1)
test_acc2=matrix(NA,length(cad_treesize),50)
train_acc2=matrix(NA,length(cad_treesize),50)
for(i in cad_treesize){
  print(i)
  for(j in 1:50){
    ind=sample(1:nrow(rft_dat),size=0.7*nrow(rft_dat),replace=F)
    tr_dat=rft_dat[ind,]
    te_dat=rft_dat[setdiff(1:nrow(rft_dat),ind),]
    tmp_rf= randomForest(res~.,tr_dat,maxnodes=i)
    ty=predict(tmp_rf,newdata=te_dat)
    oy=te_dat$res
    train_acc2[i-2,j]=sum(diag(tmp_rf$confusion[,1:2]))/nrow(tr_dat)
    test_acc2[i-2,j]=sum(diag(xtabs(~ty+oy)))/nrow(te_dat)
  }
}

windows(width=8,height=5)
par(font.main=6,font.axis=6,font.lab=6,tcl=-0.15,mgp=c(2.5,0.5,0),bty="n",mar=c(5,4,1,1))
plot(0,0,type="n",xlim=range(cad_treesize),ylim=range(train_acc2,test_acc2,1),xlab="Max Number of Terminal Nodes (Tree Depth)",
     yaxt="n",xaxt="n",ylab="Classification Accuracy",main="")
axis(2,las=1)
axis(1,at=cad_treesize)
lines(apply(train_acc2,1,mean)~cad_treesize,type="b",lty=1,pch=15, col=rgb(135,206,235,200,max=255))
lines(apply(test_acc2,1,mean)~cad_treesize,type="b",lty=2,pch=16, col=rgb(255,127,80,200,max=255))
arrows(x0=cad_treesize,y0=apply(test_acc2,1,mean)-apply(test_acc2,1,sd),y1=apply(test_acc2,1,mean)+apply(test_acc2,1,sd),
       code=3,angle=90,length=0.07,col=rgb(255,127,80,180,max=255))
legend("bottomright",bty="n",col=c(rgb(135,206,235,200,max=255),rgb(255,127,80,200,max=255)),
       lty=c(1,2),legend=c("Validation-Training","Validation-Testing"))
tun_id=which.max(apply(test_acc2,1,mean))
abline(h=apply(test_acc2,1,mean)[tun_id]-apply(test_acc2,1,sd)[tun_id],lty=4,col="lightgray")
abline(v=4,lty=4,col="lightgray")


set.seed(12345)
mod_rf = randomForest(res~., rft_dat,maxnodes=4,ntree=501)
xtabs(~predict(mod_rf)+rft_dat$res)
windows()
varImpPlot(mod_rf,main="")

rm(list=c("te_dat","tr_dat","tmp_rf","i","j","ind","oy","tun_id","ty"))

load("RF_Train_Supplemental.RData") # download this image from Github, it contains another RF trained with more similar NSO as in the testing set




