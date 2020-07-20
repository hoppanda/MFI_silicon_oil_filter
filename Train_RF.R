########################################################################
###
###    Silicon Oil Image filter - Train RF
###
###    Gregory, Chen (March 2020)
###
########################################################################


#: ------------------------------------------------------
#:  Pre-setup
#: ------------------------------------------------------
prog_path="H:/data/Projects/Research/MFI/Programs/" 
load(paste0(prog_path,"RF_Train_0.RData"))

col_name=names(rft_dat)
alls=ls(all=T)
rm(list=setdiff(alls,c("RF_train_path","col_name")))
rm(list=c("alls"))
#: ------------------------------------------------------
#:  Train RF
#: ------------------------------------------------------
RF_train_path="H:/data/Projects/Research/MFI/Data/Train Set RF Final/"  
library(readxl)
rft_oil = read_excel(paste0(RF_train_path,"SO_differentation.xlsx"),sheet="SO data set",range="A1:L2317",col_names=T) 
rft_oil = as.data.frame(rft_oil)

no_files= list.files(RF_train_path,pattern="csv$")[2]
rft_notoil = read.csv(paste0(RF_train_path,no_files),header=T)
rft_notoil$dataset=no_files

names(rft_oil)=c("Particle..","ECD.um.","Aspect.Ratio","Circularity","Intensity.Mean",
                 "Intensity.STD","Intensity.Min","Intensity.Max","Area.pixels.","Perimeter",
                 "Max.Feret.Diameter.um.","Time.Stamp.s.")
rft_notoil$res="Not Oil"
rft_oil$res="Oil"

rft_notoil = rft_notoil[,names(rft_notoil)%in%names(rft_oil)]

# set.seed(12345)
# rft_oil=rft_oil[sample(size=1500,1:nrow(rft_oil),replace=F),]

rft_dat=rbind(rft_notoil,rft_oil[,match(names(rft_notoil),names(rft_oil))])
rft_dat=rft_dat[,-c(1:2)]
rft_dat=rft_dat[rft_dat$ECD.um.>5,]
rft_dat$res=factor(rft_dat$res)
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

rm(list=c("name_check","te_dat","tr_dat","tmp_rf","i","j","ind","oy","tun_id","ty"))





