#construct functional connectivity for each subject in a given dataset
#HCP-YA minimal pre-processed data is used here as an example; 
#to analyze other datasets, the variables in lines 57-74 and 85-89 should be replaced appropriately
#——————————导入所需包——————————
library(data.table)

library(foreach)

library(doParallel)

library(doSNOW)

library(ciftiTools)

#——————————————function of constructing FC for a single fMRI data—————————————————

FC_calculator <- function(xii_path,parcellation="Schaefer_400",savefile = FALSE, filePath = NULL){
  #————————————————————parcellation————————————————————
  parc <- load_parc(parcellation)
  parc_vec <- c(as.matrix(parc))#
  n_parc <- length(unique(parc_vec))-1
  
  #……………………fMRI file………………
  xii <- read_xifti(xii_path)
  xii_mwallto0 <- move_from_mwall(xii, NA)
  xii_mat <- as.matrix(xii_mwallto0)
  
  #mean time series in each parcel
  xii_pmean <- matrix(nrow = n_parc, ncol = ncol(xii_mat))
  for (p in 1:n_parc) { 
    data_p <- xii_mat[parc_vec == p,] 
    xii_pmean[p, ] <- colMeans(data_p, na.rm=TRUE)
  }
  
  #————————————————————correlation matrix————————————————————
  full_cor <- cor(t(xii_pmean))
  
  #————————————————————save the FC matrix——————————————————
  if (savefile) {
    write.csv(full_cor, file = filePath, row.names = FALSE)
  }
  
  #————————————————————flatten the matrix——————————————————
  thisFCmat <- as.matrix(full_cor)
  
  thisFCvec <- thisFCmat[upper.tri(thisFCmat)]
  
  thisFCvec <- 0.5*log((1+thisFCvec)/(1-thisFCvec))#fisher-z transformation
  
  return(thisFCvec)
}

#——————————working folder——————————
setwd("D:\\LWfiles\\liuwei\\minimalpre")

#——————————output folder——————————
thisDataset <- "HCP1200"

fileSaveDir <- paste0("FC/",thisDataset,"/FCvec/RL/")

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

for (i in 1:length(Tasks)){
  
  thisFiledir <- Filedir[i]
  
  thisTask <- Tasks[i]
  #——————————file path, which should be changed for another dataset——————————
  File_Path <- list.dirs(thisFiledir,
                         full.names = TRUE,recursive = FALSE)
  
  subNum <- substr(File_Path,nchar(thisFiledir)+1,nchar(thisFiledir)+6)#subject number,other dataset may be different
  
  File_Path <- paste0(File_Path,"\\",subNum)
  
  fullPath <- paste0(File_Path,"\\MNINonLinear\\Results\\",thisTask,"\\",thisTask,"_Atlas_MSMAll.dtseries.nii")
  #———————————check if the file exist and remove the missing files————————
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
  
  #__________________________________________calculating________________________________________________
  
  thistaskFC <- foreach(j = 1:length(fullPath),
                        .packages = c("ciftiTools"), 
                        .combine = "rbind", 
                        .options.snow = opt)%dopar%{

                          ciftiTools.setOption('wb_path', 'D:\\programs\\workbench\\bin_windows64')
                          
                          tryCatch({
                            
                            FCvec <- FC_calculator(fullPath[j])
                            
                            return(FCvec)},
                            
                            error = function(e){
                              
                              return("Error")
                              
                            }
                          )
                        }
  close(pb)
  stopCluster(cl)
  
  #~~~~~~~~save the result~~~~~~
  FCvec <- thistaskFC
  rownames(FCvec) <- as.character(subNum)
  colnames(FCvec) <- c(1:ncol(FCvec))
  write.csv(FCvec,file = paste0(fileSaveDir,"FCvec_",thisTask,".csv"))
  
}



