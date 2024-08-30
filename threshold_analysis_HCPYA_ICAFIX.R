#slide threshold analysis
library(data.table)
library(ggplot2)
library(reshape2)

#——————————load the explained variance files————————

varDir_HCP1200 <- "D:\\LWfiles\\liuwei\\resultsICAFIX\\ICA\\HCP1200\\var\\LR"

varFiles_HCP1200 <- list.files(path = varDir_HCP1200,full.names = T)

var_R1LR <- read.csv(varFiles_HCP1200[1])

var_R2LR <- read.csv(varFiles_HCP1200[2])

var_R1LR_1IC <- var_R1LR[,2]

var_R2LR_1IC <- var_R2LR[,2]

#——————————load identification matrix——————————

id_mat <- c("D:\\LWfiles\\liuwei\\resultsICAFIX\\FC\\HCP1200\\FCcorr\\LR\\FCcorr_rfMRI_REST1_LR_rfMRI_REST2_LR.csv",
            "D:\\LWfiles\\liuwei\\resultsICAFIX\\ICA\\HCP1200\\corr\\LR\\ICAcor_rfMRI_REST1_LR_rfMRI_REST2_LR.csv")

threshold <- seq(0,0.75,0.05)

df_IDrate <- data.frame(thre = threshold)

for (j in 1:length(id_mat)) {
  
  thiscorr <- fread(id_mat[j], header = T)
  
  thiscorr <- as.matrix(thiscorr[,-1])
  
  thisIDrate <- c()
  
  nsub <- c()
  for (i in 1:length(threshold)) {
    
    index_clean <- (var_R1LR_1IC>threshold[i]) & (var_R2LR_1IC>threshold[i])
    
    thiscorr_clean <- thiscorr[index_clean,index_clean]
    
    #——————————identification rate————————
    absthisCorrfile <- abs(thiscorr_clean)
    
    withinSubCorr <- diag(absthisCorrfile)
    
    nsub[i] <- length(withinSubCorr)
    
    rowmax <- apply(absthisCorrfile, 1, max)
    
    colmax <- apply(absthisCorrfile, 2, max)
    
    
    combindedID <- sum((withinSubCorr == colmax) & (withinSubCorr == rowmax))/length(withinSubCorr)
    
    thisIDrate[i] <- combindedID
  }
  
  df_IDrate <- cbind(df_IDrate, thisIDrate)
}

colnames(df_IDrate) <- c("thre", "FC", "ICA")

df_IDrate$nsub <- nsub

df_IDrate_long <- melt(df_IDrate,id.vars = c("thre","nsub"))

colnames(df_IDrate_long) <- c("threshold", "n", "Methods","IDrate")

ggplot(df_IDrate_long,aes(x=threshold, y = IDrate, color = Methods))+
  geom_line(lwd = 1)+
  geom_point()+
  ylim(c(0.6,1.05))+
  geom_text(aes(label = round(IDrate,2)), vjust = -1, show.legend = FALSE,check_overlap = T)+
  scale_x_continuous(breaks = threshold,sec.axis = dup_axis(name = "n subjects",
                                         breaks = threshold,
                                         labels = as.character(nsub)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "D:\\LWfiles\\liuwei\\resultsICAFIX\\pic\\varAnalysis\\threshold.tiff",
       device = "tiff",dpi = 300, units = "px",
       width = 1600, height = 1600)
