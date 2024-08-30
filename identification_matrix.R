#identification matrix for each pair of fMRI conditions in a given dataset
#MIC from HCP-YA minimal pre-processed data is used here as an example; 
#to analyze other datasets, the tasks and file path should be replaced appropriately


library(data.table)#导入包，用于读取大文件


#——————————subjects for each condition——————————

thisDataset <- "HCP1200"

thisTaskTypes <- c("rfMRI_REST1_LR", "rfMRI_REST2_LR", "tfMRI_EMOTION_LR", "tfMRI_GAMBLING_LR", "tfMRI_LANGUAGE_LR",
                   "tfMRI_MOTOR_LR", "tfMRI_RELATIONAL_LR", "tfMRI_SOCIAL_LR", "tfMRI_WM_LR")
subNumlist <- list() 

for (i in 1:length(thisTaskTypes)) {#read the files to select subjects
  
  thisfilePath <- paste0("D:\\LWfiles\\liuwei\\minimalpre\\ICA\\", thisDataset, "\\comp\\LR\\FirstIC_", thisTaskTypes[i],".csv")
  
  thisdata <- fread(thisfilePath,header = TRUE, select = c(1,2))
  
  thisdata <- as.matrix(thisdata)
  
  subNumlist[[i]]<- as.numeric(thisdata[,1])
  
  assign(thisTaskTypes[i],thisdata)
}

#——————————the subjects for all conditions——————————

subs <- Reduce(intersect, subNumlist)

#——————————subjects with error FC or MIC——————————

ErrorSubs <- c()

for(i in 1:length(thisTaskTypes)){
  
  thistask <- get(thisTaskTypes[i])
  
  thisidx <- match(subs,as.numeric(thistask[,1]))
  
  ErrorSubs <- c(ErrorSubs,which(thistask[thisidx,2]=="Error"))
  
}
if (length(ErrorSubs)!=0){
  noErrorSubs <- subs[-ErrorSubs]
} else (noErrorSubs <- subs)

write.csv(noErrorSubs,file = "noErrorSubs.csv")
#——————————subjects without error in each condition——————————

rowidx <- list()

for(i in 1:length(thisTaskTypes)){
  
  thistask <- get(thisTaskTypes[i])
  
  thisidx <- match(noErrorSubs,as.numeric(thistask[,1]))
  
  rowidx[[i]] <- thisidx
  
}

names(rowidx) <- thisTaskTypes


#——————————calculating identification matrix————————

for (i in 1:length(thisTaskTypes)) {
  
  thisfilePath <- paste0("D:\\LWfiles\\liuwei\\minimalpre\\ICA\\", thisDataset, "\\comp\\LR\\FirstIC_", thisTaskTypes[i],".csv")
  
  thisdata <- fread(thisfilePath,header = TRUE)
  
  thisdata <- thisdata[rowidx[[i]],] 
  
  thisdata <- apply(thisdata, 2, as.numeric)
  
  assign(thisTaskTypes[i],thisdata)
}



for (i in 1:(length(thisTaskTypes)-1)) {
  
  for (j in (i+1):length(thisTaskTypes)) {
    
    FCtask1 <- get(thisTaskTypes[i])
    
    FCtask2 <- get(thisTaskTypes[j])
    
    thiscorr <- cor(t(FCtask1[,-1]), t(FCtask2[,-1]),use = "pairwise.complete.obs")
    
    write.csv(thiscorr,
              
              file =  paste0("D:\\LWfiles\\liuwei\\minimalpre\\ICA\\", thisDataset, "/corr/ICAcorr_", thisTaskTypes[i], "_", thisTaskTypes[j],".csv"))
  }
  
}