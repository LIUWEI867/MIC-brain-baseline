#plot the brain map of MIC

library(ciftiTools)
library(data.table)
#load a cifti file as template
xii <- read_xifti("D:\\HCPtest\\100307_3T_rfMRI_REST1_preproc\\100307\\MNINonLinear\\Results\\rfMRI_REST1_LR\\rfMRI_REST1_LR_Atlas_MSMAll.dtseries.nii")

ciftiTools.setOption('wb_path', 'D:\\programs\\workbench\\bin_windows64') 

HCPD_PA_MICmatDir <- "D:\\LWfiles\\liuwei\\resultsICAFIX\\ICA\\HCPD\\30\\comp\\PA"

HCPD_PA_MICpath <- list.files(HCPD_PA_MICmatDir,full.names = T)

picPath <- "D:\\LWfiles\\liuwei\\resultsICAFIX\\pic\\brainmap"

plotsubidx <- c(15,120,233,301,486)#choose participants

for(i in 1:length(HCPD_PA_MICpath)){
  
  MICmat <- fread(HCPD_PA_MICpath[i],header = T)#读取数据
  
  MICmat_c <- MICmat[,-1]
  
  subnum <- MICmat$V1
  
  plotsubNum <- subnum[plotsubidx]#需要绘制的被试编号
  
  for (j in 1:length(plotsubidx)) {
    
    thiscomp <- t(MICmat_c[plotsubidx[j],])
    
    thiscomp_xii <- newdata_xifti(xii,thiscomp)
    
    plot(thiscomp_xii,fname=paste0(picPath,plotsubNum[j],"_brain.png"))
  }
}



