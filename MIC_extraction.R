#extract MIC for a given dataset
#HCP-YA minimal pre-processed data is used here as an example; 
#to analyze other datasets, the variables in lines 14-35 and 78-82 should be replaced appropriately

#————————————import packages——————————
library(ciftiTools)
library(stringr)
library(ica)
library(foreach)
library(doSNOW)
library(parallel)

#——————————set working folder——————————
setwd("E:\\liuwei\\minimalpre")

#——————————output folder——————————
thisDataset <- "HCP1200"

fileSaveDir <- paste0("results/ICA/",thisDataset,"/")

#——————————file path——————————

Tasks <- c("rfMRI_REST1_LR","rfMRI_REST2_LR","tfMRI_EMOTION_LR",
           "tfMRI_GAMBLING_LR","tfMRI_MOTOR_LR","tfMRI_RELATIONAL_LR",
           "tfMRI_SOCIAL_LR","tfMRI_WM_LR","tfMRI_LANGUAGE_LR")

Filedir <- c("H:\\HCP1200\\REST1_preproc\\",
             "H:\\HCP1200\\REST2_preproc\\",
             "H:\\HCP1200\\EMOTION_preproc\\",
             "H:\\HCP1200\\GAMBLING_preproc\\",
             "H:\\HCP1200\\MOTOR_preproc\\",
             "H:\\HCP1200\\RELATIONAL_preproc\\",
             "H:\\HCP1200\\SOCIAL_preproc\\",
             "H:\\HCP1200\\WM_preproc\\",
             "H:\\HCP1200\\LANGUAGE_preproc\\")

nExtractcomp <- 30 #number of components to be extracted
nvolume = 0#volume to be analyzed, 0 refers to all volumes

#———————————function of extracing MIC and corresponding explained variance for a single fMRI data——————————

ICA_vec <- function(xii_path, ncomp = 30, nvolume = 0 ){
  #read file
  xii <- read_xifti(xii_path)
  
  xii_mat <- as.matrix(xii)
  
  #volume to be analyzed
  if (nvolume == 0) {
    nv = ncol(xii_mat)
  }else{
    nv = nvolume
  }
  
  #perform ICA
  xii_ica <- ica(xii_mat[,1:nv],nc = ncomp, method = "fast")#进行独立成份分析
  
  #explained variance
  vec_vafs <- xii_ica$vafs
  
  return(list(comp_vec=xii_ica$S[,1],vafs_vec=vec_vafs))
}

#parallel analysis for the whole dataset
for (i in 1:length(Tasks)){
  
  thisTask <- Tasks[i]
  
  thisFiledir <- Filedir[i]
  
  #——————————file path should be changed for other dataset——————————
  File_Path <- list.dirs(thisFiledir,
                         full.names = TRUE,recursive = FALSE)
    
  subNum <- substr(File_Path,nchar(thisFiledir)+1,nchar(thisFiledir)+6)#subject number,other dataset may be different
    
  File_Path <- paste0(File_Path,"\\",subNum)
  
  fullPath <- paste0(File_Path,"\\MNINonLinear\\Results\\",thisTask,"\\",thisTask,"_Atlas_MSMAll.dtseries.nii")
  
  #——————————check if the file exist and remove the missing files————————
  exitIndx <- file.exists(fullPath)
  
  subNum <- subNum[exitIndx]
  
  fullPath <- fullPath[exitIndx]
  
  
  #~~~~~~~~parallel analysis~~~~~~~~
  ncores <- detectCores() 
  
  cl <- makeCluster(ncores-2)
  
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(fullPath),style = 3)
  
  progress <- function(n){setTxtProgressBar(pb, n)}
  
  opt <- list(progress = progress)
  
  #~~~~~~~~~calculating~~~~~~~
  
  
  thistaskICA <- foreach(j = 1:length(fullPath),
                         
                         .packages = c("ica","ciftiTools"), 
                         .combine = "rbind", 
                         .options.snow = opt)%dopar%{
                           
                           ciftiTools.setOption('wb_path', 'D:\\Program Files (x86)\\workbench-windows64-v1.5.0\\workbench\\bin_windows64')
                           
                           tryCatch({icaResults <- ICA_vec(fullPath[j],ncomp = nExtractcomp, nvolume = nvolume)#最大独立成分
                           
                           cv1 <- icaResults[[1]]
                           
                           cv2 <- icaResults[[2]]
                           
                           results <- c(cv2, cv1)},
                           
                           error = function(e){
                             
                             return("Error")
                             
                           }
                           )
                         }
  close(pb)
  stopCluster(cl)
  
  #~~~~~~~~save the results~~~~~~
  #explained variance
  vafs_all <- thistaskICA[,1:nExtractcomp]
  rownames(vafs_all) <- as.character(subNum)
  colnames(vafs_all) <- c(1:ncol(vafs_all))
  write.csv(vafs_all,file = paste0(fileSaveDir,"var/vafs_",thisTask,".csv"))
  
  #MIC
  ICA_all <- thistaskICA[,(nExtractcomp+1):ncol(thistaskICA)]
  rownames(ICA_all) <- as.character(subNum)
  colnames(ICA_all) <- c(1:ncol(ICA_all))
  write.csv(ICA_all,file = paste0(fileSaveDir,"comp/FirstIC_",thisTask,".csv"))
}


