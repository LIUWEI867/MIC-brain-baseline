#plot the histgram of the explained variance of MIC
library(ggplot2)

library(patchwork)
#——————————HCP-D——————————
varDir_HCPD <- "D:\\LWfiles\\liuwei\\resultsICAFIX\\ICA\\HCPD\\30\\var\\PA"

varFiles_HCPD <- list.files(path = varDir_HCPD,full.names = T)

HCPDtasks <- c("HCP-D RS1", "HCP-D RS2", "HCP-D CAR", "HCP-D EMO", "HCP-D GUE")

#——————————HCP-A——————————
varDir_HCPA <- "D:\\LWfiles\\liuwei\\resultsICAFIX\\ICA\\HCPA\\30\\var\\PA"

varFiles_HCPA <- list.files(path = varDir_HCPA,full.names = T)

HCPAtasks <- c("HCP-A RS1", "HCP-A RS2", "HCP-A CAR", "HCP-A FN", "HCP-A VM")

#——————————HCP-YA——————————
varDir_HCP1200 <- "D:\\LWfiles\\liuwei\\resultsICAFIX\\ICA\\HCP1200\\var\\LR"

varFiles_HCP1200 <- list.files(path = varDir_HCP1200,full.names = T)

HCP1200tasks <- c("HCP-YA RS1", "HCPA-YA RS2")

#——————————pic output path——————————

picdir <- "D:\\LWfiles\\liuwei\\resultsICAFIX\\pic\\varAnalysis\\"

#——————————plot——————————

alltasks <- c(HCPAtasks, HCPDtasks,HCP1200tasks)

allvarfiles <- c(varFiles_HCPA, varFiles_HCPD, varFiles_HCP1200)

picnames <- paste0("p", c(1:length(alltasks))) #绘图的名字

for (i in 1:length(allvarfiles)){
  
  thisVardf <- read.csv(allvarfiles[i])
  
  thisVar <- data.frame(vars = thisVardf[,2])
  
  VarDensity <- density(thisVar$vars)
  
  DesR <- range(VarDensity$y)
  
  thisp <- ggplot(thisVar, aes(vars))+
    geom_histogram(aes(y=after_stat(density)),fill = "red",alpha = 0.3)+
    geom_density(color = "skyblue", fill = "skyblue", alpha = 0.3)+
    coord_cartesian(xlim = c(0,1))+
    annotate("text", 0.25, DesR[2]*0.8, size = 3,
            label = paste0("range:\n [", round(min(thisVar$vars),2),", ",round(max(thisVar$vars),2),"]"))+
    theme_classic()+
    theme(axis.title = element_blank(),
          plot.title = element_text(size = 12))

  
  assign(picnames[i], thisp)
  
  ggsave(filename = paste0(picdir,picnames[i],".tiff"),
         device = "tiff",dpi = 300, units = "px",
         width = 600, height = 600)
  
}


  
