#identification rate after thresholding the explained variance for HCP-YA minimal data
library(data.table)
library(ggplot2)
library(reshape2)
library(corrplot)

#——————————load explained variance————————

varDir_HCP1200 <- "D:\\LWfiles\\liuwei\\minimalpre\\ICA\\HCP1200\\var\\LR"

varFiles_HCP1200 <- list.files(path = varDir_HCP1200,full.names = T)

thistasks <- c("RS1","RS2","EMO","GAM","LAN","MOT","REL","SOC","WM") 

for (i in 1:length(thistasks)) {
  
  thisvar <- read.csv(varFiles_HCP1200[i])
  
  assign(thistasks[i], thisvar)
}

#——————————load identification matrices——————————

ICA_id_mat_path <- "D:\\LWfiles\\liuwei\\minimalpre\\ICA\\HCP1200\\corr\\LR"

ICA_id_mat <- list.files(path = ICA_id_mat_path,full.names = TRUE)

FC_id_mat_path <- "D:\\LWfiles\\liuwei\\minimalpre\\FC\\HCP1200\\FCcorr\\LR"

FC_id_mat <- list.files(path = FC_id_mat_path,full.names = TRUE)

threshold <- 0.65

df_IDrate <- data.frame(pairs = substr(ICA_id_mat,nchar(ICA_id_mat)-34, nchar(ICA_id_mat)-4))

#——————load the subject NO.————————
noErrorSubs <-read.csv("D:\\LWfiles\\liuwei\\minimalpre\\ICA\\HCP1200\\NoErrorSubs.csv",header = T)

allsubs <- noErrorSubs$x

#——————按照阈值提取被试——————

i <- 1

thisICAIDrate <- c()

thisFCIDrate <- c()

nsub <- c()

for (j in 1:(length(varFiles_HCP1200)-1)) {
  
  for (k in (j+1):length(varFiles_HCP1200)) {
    
    thisICAcorr <- fread(ICA_id_mat[i], header = T)
    
    thisICAcorr <- as.matrix(thisICAcorr[,-1])
    
    thisFCcorr <- fread(FC_id_mat[i], header = T)
    
    thisFCcorr <- as.matrix(thisFCcorr[,-1])


    varfile1 <- get(thistasks[j])
    
    varfile2 <- get(thistasks[k])
    

    sub1 <- varfile1[,1]
    
    sub2 <- varfile2[,1]

    
    var1 <- varfile1[match(allsubs,sub1),]
    
    var2 <- varfile2[match(allsubs,sub2),]
    
    var1 <- var1[,2]
    
    var2 <- var2[,2]
    
    index_clean <- (var1>threshold) & (var2>threshold)#大于阈值的被试数
    
    nsub[i] <- sum(index_clean)
    
    thisICAcorr_clean <- thisICAcorr[index_clean,index_clean]
    
    thisFCcorr_clean <- thisFCcorr[index_clean,index_clean]
    
    #——————————MIC identification rate————————
    absthisCorrfile <- abs(thisICAcorr_clean)
    
    withinSubCorr <- diag(absthisCorrfile)

    rowmax <- apply(absthisCorrfile, 1, max)
    
    colmax <- apply(absthisCorrfile, 2, max)
    
    combindedID <- sum((withinSubCorr == colmax) & (withinSubCorr == rowmax))/length(withinSubCorr)
    
    thisICAIDrate[i] <- combindedID
    
    #——————————FC identification rate————————
    absthisCorrfile <- abs(thisFCcorr_clean)
    
    withinSubCorr <- diag(absthisCorrfile)
    
    rowmax <- apply(absthisCorrfile, 1, max)
    
    colmax <- apply(absthisCorrfile, 2, max)
    
    
    combindedID <- sum((withinSubCorr == colmax) & (withinSubCorr == rowmax))/length(withinSubCorr)
    
    thisFCIDrate[i] <- combindedID
    
    i <- i+1
  }
}

  
df_IDrate$ICAIDrate <- thisICAIDrate

df_IDrate$FCIDrate <- thisFCIDrate

df_IDrate$nsub <- nsub


#——————————————————————heatmap————————————————————


heatmat <- matrix(nrow = length(thistasks), ncol = length(thistasks))

rownames(heatmat) <- thistasks
colnames(heatmat) <- thistasks

heatmat[lower.tri(heatmat)] <- df_IDrate$ICAIDrate

heatmat <- t(heatmat)

heatmat[lower.tri(heatmat)] <- df_IDrate$FCIDrate

submat <- matrix(nrow = length(thistasks), ncol = length(thistasks))

rownames(submat) <- thistasks
colnames(submat) <- thistasks

submat[lower.tri(submat)] <- df_IDrate$nsub

submat <- t(submat)

submat[lower.tri(submat)] <- df_IDrate$nsub


colp1 <- colorRampPalette(c("blue","yellow","red")) 

picPath <- "D:\\LWfiles\\liuwei\\resultsMINIMAL\\pic\\MICvardis\\pic\\"

tiff(filename = paste0(picPath,"HCP1200thre",threshold," idmat.tiff"),width = 1200,height = 1200)

corrplot(heatmat,
         col =  colp1(200),
         col.lim = c(0,1),
         cl.cex = 1.4,
         #tl.pos = "n",
         type = "full",
         addCoef.col = "black",
         method = "shade",
         tl.col = "black", tl.cex = 2, 
         number.cex = 3,
         diag = FALSE)

dev.off()

#subject number

tiff(filename = paste0(picPath,"HCP1200thre",threshold," subNum.tiff"),width = 1200,height = 1200)

corrplot(submat,
         col =  NULL,
         is.corr = FALSE,
         cl.pos = "n",
         cl.cex = 1.4,
         #tl.pos = "n",
         type = "full",
         addCoef.col = "black",
         method = "shade",
         tl.col = "black", tl.cex = 2, 
         number.cex = 3,
         diag = FALSE)

dev.off()