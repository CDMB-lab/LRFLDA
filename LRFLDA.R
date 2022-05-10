library(openxlsx)
library(randomForest)
library(ROCR)
library(plyr)
library(parallel)
###################################read data###########################################
L <- read.xlsx("/home/ubuntu/Documents/R/01-lncRNAs-240.xlsx", sheet = 1, colNames = FALSE)
D <- read.xlsx("/home/ubuntu/Documents/R/02-diseases-412.xlsx", sheet = 1, colNames = FALSE)
M <- read.xlsx("/home/ubuntu/Documents/R/03-miRNAs-495.xlsx", sheet = 1, colNames = FALSE)
LL <- read.xlsx("/home/ubuntu/Documents/R/04-lncRNA-lncRNA.xlsx", sheet = 1, colNames = FALSE)
LD <- read.xlsx("/home/ubuntu/Documents/R/05-lncRNA-disease.xlsx", sheet = 1, colNames = FALSE)
MD <- read.xlsx("/home/ubuntu/Documents/R/06-miRNA-disease.xlsx", sheet = 1, colNames = FALSE)
DD <- read.xlsx("/home/ubuntu/Documents/R/07-disease-disease.xlsx", sheet = 1, colNames = FALSE)
LM <- read.xlsx("/home/ubuntu/Documents/R/08-lncRNA-miRNA.xlsx", sheet = 1, colNames = FALSE)
###################################feature matrix###########################################
L1147 <- cbind(LL, LM, LD)
L1148 <- cbind(L[,1], L1147)
DL <- t(LD)
DM <- t(MD)
D1147 <- cbind(DD, DM, DL)
D1148 <- cbind(D[,1], D1147)
LDALL <- merge(x = L1148, y = D1148, by = NULL)
d1 <- subset(LDALL,select=1)
d2 <- subset(LDALL,select=1149)
d3 <- subset(LDALL,select=c(-1,-1149))
LDALL <-data.frame(d1,d2,d3)
for(i in 1:ncol(LDALL)){
  colnames(LDALL)[i] <- c(paste("X",i,sep = ""))
}

###################################non-zero###########################################
d1 <- subset(LDALL,select=c(1,2))
d2 <- subset(LDALL,select=c(-1,-2))
d2 <- d2[,which(colSums(d2) > 0)] 
LDExcl0 <- cbind(d1, d2)


###################################samples###########################################
xy <- which(LD[,] == 1, arr.ind = TRUE)
LDA <- cbind(L[xy[,1],1],D[xy[,2],1])
label.fr <- data.frame("label" = c(0))
LDExcl0 <- cbind(label.fr, LDExcl0)
for (i in 1:nrow(LDExcl0))
{
  for (j in 1:nrow(LDA))
  {
    if((LDExcl0[i,2] == LDA[j,1]) && (LDExcl0[i,3] == LDA[j,2]))
    {
      LDExcl0[i,1] <- 1
    }
  }    
}

PositiveSample <- LDExcl0[LDExcl0[,1]==1, ]
UnlabeledSample <- LDExcl0[LDExcl0[,1]==0, ]
B11<-subset(UnlabeledSample[,], select=-(X2:X3))
sp <- sample(nrow(UnlabeledSample), nrow(PositiveSample), replace = FALSE, prob = NULL)
sps<-sort(sp)
NegativeSample <- UnlabeledSample[sps,]
TrainingSample <- rbind(PositiveSample, NegativeSample)
##############################################################################################################

#########################################de-redundancy ################################
A1 <- unlist(subset(TrainingSample[,], select=(X1)))
A2 <- (subset(TrainingSample[,], select=-(X1:X3)))

for(i in 1:ncol(A2)){
  colnames(A2)[i] <- c(paste("X",i,sep = ""))
  
}
{Lasso<-Lasso(A2,A1,lambda=0.0004,fix.lambda=TRUE,cv.OLS=FALSE)
  c<-colnames(A2[,which(Lasso$beta!=0)])
  xx_lasso<-A2[,which(Lasso$beta!=0)]   
  CC1<-xx_lasso
  CC1 <- cbind(A1, CC1)
  C<-colnames(xx_lasso)
  CC2<-subset(UnlabeledSample[,], select=(C))
  A3<-unlist(subset(UnlabeledSample[,], select=(X1)))
  label.fr <- data.frame("label" = c(0))
  CC2 <- cbind(label.fr, CC2)
  for(i in 1:ncol(CC1)){
    colnames(CC1)[i] <- c(paste("X",i,sep = ""))
    
  }
  for(i in 1:ncol(CC2)){
    colnames(CC2)[i] <- c(paste("X",i,sep = ""))
    
  }
  T1<-CC1
  T2<-CC2
}
write.xlsx(T1, "/home/ubuntu/Documents/R/TrainingSample_Lasso.xlsx",colNames = FALSE)
#################################### 5fold###############################    
{ 

TB <- T1
NB <- T2   

sumauc <- 0
sumap <- 0
PTB <- TB[1:2697,]
ind <- sample (2697, 2697, replace =FALSE)
PTB11 <- PTB[ind,]
NTB <- TB[2698:5394,] 	
ind <- sample (2697, 2697, replace =FALSE)
NTB11 <- NTB[ind,]
}
{
  num <- c(1:4)
  cl <- makeCluster(20) # 初始化四核心集群
  clusterExport(cl,"PTB11",envir = environment())
  clusterExport(cl,"NTB11",envir = environment())
  clusterExport(cl,"NB",envir = environment())
  clusterExport(cl,"sumauc",envir = environment())
  clusterExport(cl,"sumap",envir = environment())
  #print(system.time({
}  
Fun<-function(i)
{
  ### training sample set
  PTB1 <- PTB11[-(((540*(i-1))+1):(540*i)),]
  NTB1 <- NTB11[-(((540*(i-1))+1):(540*i)),]
  TrainB <- rbind(PTB1, NTB1)
  
  ### test sample set
  PTB2 <- PTB11[(((540*(i-1))+1):(540*i)),]
  TestB <- rbind(PTB2,NB)
  library(randomForest)
  library(ROCR)
  rf=randomForest(X1~.,data = TrainB, mtry=172, importance = TRUE, ntree=500, na.action=na.omit)
  pred <- predict(rf, TestB)
  pred1 <- prediction(pred, TestB$X1)
  roc <- performance(pred1, "tpr", "fpr")
  plot(roc, main = "ROC chart")
  auc <- performance(pred1, "auc")@y.values
  print(auc)
  sumauc <- sumauc + as.numeric(auc[[1]])
  perf1 <- performance(pred1, "prec", "rec")
  plot(perf1)
  prec <- performance(pred1, "prec")@y.values
  rec <- performance(pred1, "rec")@y.values
  ap <-0
  cur_rec <- rec[[1]][2]
  cur_prec <- prec[[1]][2]     
  for (j in 3:length(rec[[1]])) {
    if(prec[[1]][j] >= cur_prec)
    {
      cur_prec = prec[[1]][j]
    }
    if (abs(cur_rec - rec[[1]][j]) > 0) {
      ap = ap + cur_prec * abs(cur_rec - rec[[1]][j])
    }
    cur_rec = rec[[1]][j]
  }
  print(ap)
  sumap <- sumap + ap
  
  return(list(auc, ap))
}
results <- parLapply(cl,num,Fun)#调用parLapply并行计算平方函数
final <- do.call('c',results)#整合结果
final
stopCluster(cl)

{i <- 5
  PTB1 <- PTB11[-(((540*(i-1))+1):2697),]
  NTB1 <- NTB11[-(((540*(i-1))+1):2697),]
  TrainB <- rbind(PTB1, NTB1)
  
  ### test sample set
  PTB2 <- PTB11[(((540*(i-1))+1):2697),]
  TestB <- rbind(PTB2,NB)
  rf=randomForest(X1~.,data = TrainB, mtry=172, importance = TRUE, ntree=500, na.action=na.omit)
  ### predict using RandomForest
  pred <- predict(rf, TestB)
  pred1 <- prediction(pred, TestB$X1)
  roc <- performance(pred1, "tpr", "fpr")
  plot(roc, main = "ROC chart")
  
  auc <- performance(pred1, "auc")@y.values
  print(auc)
  sumauc <- sumauc + as.numeric(auc[[1]])
  sumauc <- sumauc/5
  #print(sumauc)
  
  perf1 <- performance(pred1, "prec", "rec")
  plot(perf1)
  
  prec <- performance(pred1, "prec")@y.values
  rec <- performance(pred1, "rec")@y.values
  ap <-0
  cur_rec <- rec[[1]][2]
  cur_prec <- prec[[1]][2]     
  for (j in 3:length(rec[[1]])) {
    if(prec[[1]][j] >= cur_prec)
    { 
      cur_prec = prec[[1]][j]
    }
    if (abs(cur_rec - rec[[1]][j]) > 0) {
      ap = ap + cur_prec * abs(cur_rec - rec[[1]][j])
    }
    cur_rec = rec[[1]][j]
  }
  print(ap)
  sumap <- sumap + ap
  sumap <- sumap/5
  # print(sumap) 
  print(n)
}
##############################################################################################################
################################### predict all lncRNA-disease samples #######################
TB <- T1
NB <-T2
rf=randomForest(X1~.,data = TB, mtry=172, importance = TRUE, ntree=500, na.action=na.omit)
pred <- predict(rf, NB)
NB1 <-UnlabeledSample  [,1:3]
UnlabeledSampleScore <- cbind(NB1, data.frame(pred))
write.xlsx(UnlabeledSampleScore, "/home/ubuntu/Documents/R/UnlabeledSampleScore-Lasso.xlsx",colNames = FALSE)
################################################################################################################
