#identification rate for each identification matrix in a given dataset
#MIC identification matrix from HCP-YA minimal pre-processed data is used here as an example; 
#to analyze other datasets, the tasks and file path should be replaced appropriately
library(data.table)
library(pheatmap)
library(corrplot)
library(ggplot2)
library(ggsignif)

#——————————设置基本参数————————
#path to identification matrix
filePath <- "D:\\LWfiles\\liuwei\\minimalpre\\ICA\\HCP1200\\corr\\LR"

#output path for pictures
picPath <- "D:\\LWfiles\\liuwei\\resultsMINIMAL\\pic\\ICA_HCP1200_LR\\pic\\test\\"

#task names in this dataset. note the order should be in line with the file order loaded in line 35
thistasks <- c("RS1","RS2","EMO","GAM","LAN","MOT","REL","SOC","WM") 

#name the pairs
taskPairs <- c()
k <- 1

for (i in 1:(length(thistasks)-1)) {
  
  for (j in (i+1):length(thistasks)) {
    
    taskPairs[k] <- paste0(thistasks[i]," v.s. ",thistasks[j])
    
    k <- k+1
  }
  
}

taskFiles <- list.files(path = filePath, full.names = TRUE)#获取文件完整路径

#———————————functions——————————————
#function of calculating the identification rate for a given identification matrix

IDrateCalculator <- function(filePath, fileName, savedir){
  
  #——————————read identification matrix file——————————
  thisCorrfile <- fread(filePath,header = TRUE)
  
  thisCorrfile <- as.matrix(thisCorrfile[,-1])
  #——————————heat map————————
  
  pheatmap(thisCorrfile, cluster_rows = FALSE,cluster_cols = FALSE,
           breaks = seq(-1,1,0.02),
           legend = FALSE,
           show_rownames = FALSE,show_colnames = FALSE,
           filename = paste0(savedir, fileName,".tiff"))
  
  #——————————identification rate——————————
  
  absthisCorrfile <- as.matrix(abs(thisCorrfile))
  
  withinSubCorr <- diag(absthisCorrfile)
  
  rowmax <- apply(absthisCorrfile, 1, max)
  
  colmax <- apply(absthisCorrfile, 2, max)
  
  rowID <- sum(withinSubCorr == rowmax)/length(withinSubCorr)
  
  colID <- sum(withinSubCorr == colmax)/length(withinSubCorr)
  
  combindedID <- sum((withinSubCorr == colmax) & (withinSubCorr == rowmax))/length(withinSubCorr)
  
  #——————————violin plot————————
  vioPlot(thisCorrfile,fileName,savedir)
  
  return(c(fileName,rowID, colID,combindedID))
}

#function of plotting identification power violin plot

vioPlot <- function(cormat, fileName, savedir){
  #——————————intra-sub correlation————————
  
  within_score <- diag(cormat)
  
  within_score_z <- abs(0.5*log((1+within_score)/(1-within_score)))#fisher z transform
  
  #——————————inter-sub correlation————————
  between_score=c()
  for (i in 1:ncol(cormat)) {
    row_betweenScore <- cormat[i,-i]
    col_betweenScore <- cormat[-i,i]
    between_score[i] <- max(abs(c(row_betweenScore,col_betweenScore)))
  }
  
  between_score_z <- abs(0.5*log((1+between_score)/(1-between_score)))#fisher z transform
  
  #data for plot
  idData <- data.frame(score = c(within_score_z,between_score_z),
                       group = rep(c("intra","inter"),each = ncol(cormat)))
  
  t.result <- t.test(within_score_z,between_score_z,paired = T)
  
  CI95 = round(t.result$conf.int,digits = 2)
  
  cohenD <- round(t.result$estimate/sd(within_score_z-between_score_z),2)
  
  ggplot(idData,aes(x = group, y = score, fill = group))+
    
    geom_violin(alpha=0.6)+
    
    scale_fill_manual(values = c("blue","red"))+
    
    geom_boxplot(alpha=1,size=0.3, width=0.2,fill="white")+
    
    stat_summary(fun="mean",geom="point",shape=21, size=2,fill="blue")+
    
    geom_signif(comparisons = list(c("intra","inter")),
                map_signif_level = T,
                test = t.test, test.args = c(paried = TRUE),
                size = .8,textsize = 10)+
    
    scale_y_continuous(limits = c(-0.4,max(idData$score)+0.2))+
    
    annotate("text", 1.5, min(c(between_score_z,within_score_z))-0.2, size = 10,
             label = paste0("95%CI \n[", CI95[1],", ", CI95[2],"]"))+
    
    annotate("text", 1.5, max(c(between_score_z,within_score_z))-0.1, size = 10,
             label = paste0("d = ",cohenD))+
    
    theme_classic()+
    theme(axis.title = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=20),
          legend.position = "none",#去掉图例
    )
  
  ggsave(filename = paste0(savedir, fileName,"vioPlot.tiff"), device = "tiff",
         width = 1600, height = 1600, units = "px")
}

#identification results
IDdf <- data.frame(taskPairs = "", rowID = 0,  colID = 0, comID = 0)#设置空列表用于储存识别率

for (i in 1:length(taskFiles)) {
  
  thisresults <- IDrateCalculator(taskFiles[i], fileName = taskPairs[i], savedir = picPath)
  
  IDdf[i,] <- thisresults
  
}

#————————heatmap of identification——————————

maintitle <- "HCP-YA LR"

ndim <- length(thistasks) 

comidmat <- matrix(nrow = ndim, ncol = ndim)

comidmat[lower.tri(comidmat)] <- as.numeric(IDdf[,4])

rownames(comidmat) <- thistasks
colnames(comidmat) <- thistasks

colp1 <- colorRampPalette(c("blue","yellow","red")) #设置热值图颜色

tiff(filename = paste0(picPath,maintitle," idmat.tiff"),width = 800,height = 800)

corrplot(comidmat,
         col =  colp1(200),
         col.lim = c(0,1),
         cl.cex = 2,
         type = "lower",
         addCoef.col = "black",
         method = "shade",
         #tl.pos = "n",
         tl.col = "black", tl.cex = 3,
         number.cex = 3,
         diag = FALSE)
dev.off()

