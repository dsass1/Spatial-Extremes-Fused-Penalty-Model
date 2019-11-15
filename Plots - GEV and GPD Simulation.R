#Plots for GEV and GPD simulation study
##########################################################################################
##########################################################################################
##########################################################################################
#directory with stationary GEV results
dir.gev <- "~/STAT - Research/Extremes Project/Final Paper/GEV_Final_Paper_Data/Results/"
#directory with Schlather results
dir.max <- "~/STAT - Research/Extremes Project/Final Paper/Schlather/Results/"
#directory with nonstationary GEV results
dir.gev2 = "~/STAT - Research/Extremes Project/Final Paper/NonStationary/7500 Iter 5000 Burn/Results/"
#directory with independent simulation results for Appendix
dir.ind <- "~/STAT - Research/Extremes Project/Final Paper/GEV_Final_Paper_Data/"
#directory with stationary GPD results
dir.gpd <- "~/STAT - Research/Extremes Project/Final Paper/GPD_Final_Paper_Data/Stationary/"
#directory with nostationary GPD results
dir.gpd2 <- "~/STAT - Research/Extremes Project/Final Paper/GPD_Final_Paper_Data/NonStationary/"

#set directory where plots should be saved
plot.dir = "~/STAT - Research/Extremes Project/Final Paper/BoxPlots/"


##GEV Results
#load results
files = list.files(dir.gev)
gev_data = read.csv(paste0(dir.gev, files[1]))
for(f in 2:length(files)) {
  gev_data = rbind(gev_data, read.csv(paste0(dir.gev, files[f])))
}

files = list.files(dir.max)
gev_max = read.csv(paste0(dir.max, files[1]))
for(f in 2:length(files)) {
  gev_max = rbind(gev_max, read.csv(paste0(dir.max, files[f])))
}

##################################################################################
#Table of averages for appendix
avgloc <- c(mean(gev_data$MSE.Loc.Spat), mean(gev_max$MSE.Loc.Max), mean(gev_data$MSE.Loc.Bayes), mean(gev_data$MSE.Loc.Ridge),mean(gev_data$MSE.Loc.Lasso))
avgscale <- c(mean(gev_data$MSE.Scale.Spat), mean(gev_max$MSE.Scale.Max), mean(gev_data$MSE.Scale.Bayes), mean(gev_data$MSE.Scale.Ridge),mean(gev_data$MSE.Scale.Lasso))
avgshape <- c(mean(gev_data$MSE.Shape.Spat), mean(gev_max$MSE.Shape.Max), mean(gev_data$MSE.Shape.Bayes), mean(gev_data$MSE.Shape.Ridge),mean(gev_data$MSE.Shape.Lasso))
avgRL10 <- c(mean(gev_data$MSE.RL10.Spat), mean(gev_max$MSE.RL10.Max), mean(gev_data$MSE.RL10.Bayes), mean(gev_data$MSE.RL10.Ridge),mean(gev_data$MSE.RL10.Lasso))
avgRL20 <- c(mean(gev_data$MSE.RL20.Spat), mean(gev_max$MSE.RL20.Max), mean(gev_data$MSE.RL20.Bayes), mean(gev_data$MSE.RL20.Ridge),mean(gev_data$MSE.RL20.Lasso))
avgRL50 <- c(mean(gev_data$MSE.RL50.Spat), mean(gev_max$MSE.RL50.Max),mean(gev_data$MSE.RL50.Bayes), mean(gev_data$MSE.RL50.Ridge),mean(gev_data$MSE.RL50.Lasso))
avgRL100 <- c(mean(gev_data$MSE.RL100.Spat), mean(gev_max$MSE.RL100.Max),mean(gev_data$MSE.RL100.Bayes), mean(gev_data$MSE.RL100.Ridge),mean(gev_data$MSE.RL100.Lasso))
time <- c(mean(gev_data$Time.Spat), mean(gev_max$Time.Max),mean(gev_data$Time.Bayes), mean(gev_data$Time.Ridge),mean(gev_data$Time.Lasso))
fail <- c(mean(gev_data$Fail.Spat), 5,0,mean(gev_data$Fail.Ridge), mean(gev_data$Fail.Lasso)) 
GEV_avgMSE <- matrix(c(avgloc,
                       avgscale,
                       avgshape,
                       avgRL10,
                       avgRL20,
                       avgRL50,
                       avgRL100,
                       time,
                       fail), ncol = 5, byrow = TRUE)
colnames(GEV_avgMSE)<-c("Spat","Schlather", "Bayes","Ridge","Lasso")
rownames(GEV_avgMSE)<-c("Loc","Scale","Shape","RL10","RL20","RL50","RL100","Time", "Fail")
GEV_MSE <- round(GEV_avgMSE,4)

GEV_RMSE <- sqrt(GEV_MSE)
GEV_RMSE
##################################################################################
#Final boxplots with ggplot
##################################################################################
root.loc <- c(sqrt(gev_data$MSE.Loc.Spat),sqrt(gev_data$MSE.Loc.Ridge),sqrt(gev_data$MSE.Loc.Lasso),sqrt(gev_data$MSE.Loc.Bayes), sqrt(gev_max$MSE.Loc.Max))
root.scale<- c(sqrt(gev_data$MSE.Scale.Spat),sqrt(gev_data$MSE.Scale.Ridge),sqrt(gev_data$MSE.Scale.Lasso), sqrt(gev_data$MSE.Scale.Bayes), sqrt(gev_max$MSE.Scale.Max))
root.shape <- c(sqrt(gev_data$MSE.Shape.Spat),sqrt(gev_data$MSE.Shape.Ridge),sqrt(gev_data$MSE.Shape.Lasso), sqrt(gev_data$MSE.Shape.Bayes), sqrt(gev_max$MSE.Shape.Max))
root.rl10 <-c(sqrt(gev_data$MSE.RL10.Spat),sqrt(gev_data$MSE.RL10.Ridge),sqrt(gev_data$MSE.RL10.Lasso), sqrt(gev_data$MSE.RL10.Bayes), sqrt(gev_max$MSE.RL10.Max))
root.rl20 <-c(sqrt(gev_data$MSE.RL20.Spat),sqrt(gev_data$MSE.RL20.Ridge),sqrt(gev_data$MSE.RL20.Lasso), sqrt(gev_data$MSE.RL20.Bayes), sqrt(gev_max$MSE.RL20.Max))
root.rl50 <- c(sqrt(gev_data$MSE.RL50.Spat),sqrt(gev_data$MSE.RL50.Ridge),sqrt(gev_data$MSE.RL50.Lasso), sqrt(gev_data$MSE.RL50.Bayes), sqrt(gev_max$MSE.RL50.Max))
root.rl100 <- c(sqrt(gev_data$MSE.RL100.Spat),sqrt(gev_data$MSE.RL100.Ridge),sqrt(gev_data$MSE.RL100.Lasso), sqrt(gev_data$MSE.RL100.Bayes), sqrt(gev_max$MSE.RL100.Max))
model <- c(rep("Spatial",length(gev_data$MSE.Loc.Spat)),rep("Ridge",length(gev_data$MSE.Loc.Ridge)),rep("Lasso",length(gev_data$MSE.Loc.Lasso)),rep("Bayes",length(gev_data$MSE.Loc.Bayes)),rep("Schlather",length(gev_max$MSE.Loc.Max)))
#gev_plot_data <- data.frame(RootLoc = root.loc, RootScale = root.scale, RootShape = root.shape, RootRL10 = root.rl10, RootRL20 = root.rl20, RootRL50=root.rl50, RootRL100=root.rl100, Model = model)

log.root.loc <- log(root.loc)
log.root.scale <- log(root.scale)
log.root.shape <- log(root.shape)
log.root.rl10 <- log(root.rl10)
log.root.rl20 <- log(root.rl20)
log.root.rl50 <- log(root.rl50)
log.root.rl100 <- log(root.rl100)
gev_plot_data <- data.frame(RootLoc = log.root.loc, RootScale = log.root.scale, RootShape = log.root.shape, RootRL10 = log.root.rl10, RootRL20 = log.root.rl20, RootRL50= log.root.rl50, RootRL100=log.root.rl100, Model = model)


#reorder data: Spatial 5, Schlather 4, Bayes 1, Ridge 3, Lasso 2
gev_plot_data$Model = factor(gev_plot_data$Model, levels = levels(gev_plot_data$Model)[c(5,4,1,3,2)])

ggplot(gev_plot_data, aes(x=Model, y=RootLoc))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="Location", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGEV_Loc.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootScale))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="Scale", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGEV_Scale.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootShape))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="Shape", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGEV_Shape.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootRL10))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="10-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGEV_RL10.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootRL20))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="20-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGEV_RL20.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootRL50))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="50-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGEV_RL50.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootRL100))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="100-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGEV_RL100.png"), width = 4, height = 6)

##########################################################################################
##########################################################################################
#GEV NONSTATIONARY
#nonstationary: 5000 burn in 7500 iterations
files = list.files(dir.gev2)
nonstat_data2 = read.csv(paste0(dir.gev2, files[1]))
for(f in 2:length(files)) {
  nonstat_data2 = rbind(nonstat_data2, read.csv(paste0(dir.gev2, files[f])))
}

avgloc <- c(mean(nonstat_data2$MSE.Loc.Spat), mean(nonstat_data2$MSE.Loc.Spat2),mean(nonstat_data2$MSE.Loc.Max), mean(nonstat_data2$MSE.Loc.Bayes), mean(nonstat_data2$MSE.Loc.Ridge),mean(nonstat_data2$MSE.Loc.Lasso))
avgscale <- c(mean(nonstat_data2$MSE.Scale.Spat), mean(nonstat_data2$MSE.Scale.Spat2),mean(nonstat_data2$MSE.Scale.Max),mean(nonstat_data2$MSE.Scale.Bayes),mean(nonstat_data2$MSE.Scale.Ridge),mean(nonstat_data2$MSE.Scale.Lasso))
avgshape <- c(mean(nonstat_data2$MSE.Shape.Spat), mean(nonstat_data2$MSE.Shape.Spat2),mean(nonstat_data2$MSE.Shape.Max),mean(nonstat_data2$MSE.Shape.Bayes),mean(nonstat_data2$MSE.Shape.Ridge),mean(nonstat_data2$MSE.Shape.Lasso))
avgRL10 <- c(mean(nonstat_data2$MSE.RL10.Spat), mean(nonstat_data2$MSE.RL10.Spat2),mean(nonstat_data2$MSE.RL10.Max),mean(nonstat_data2$MSE.RL10.Bayes),mean(nonstat_data2$MSE.RL10.Ridge),mean(nonstat_data2$MSE.RL10.Lasso))
avgRL20 <- c(mean(nonstat_data2$MSE.RL20.Spat), mean(nonstat_data2$MSE.RL20.Spat2),mean(nonstat_data2$MSE.RL20.Max),mean(nonstat_data2$MSE.RL20.Bayes),mean(nonstat_data2$MSE.RL20.Ridge),mean(nonstat_data2$MSE.RL20.Lasso))
avgRL50 <- c(mean(nonstat_data2$MSE.RL50.Spat), mean(nonstat_data2$MSE.RL50.Spat2),mean(nonstat_data2$MSE.RL50.Max),mean(nonstat_data2$MSE.RL50.Bayes),mean(nonstat_data2$MSE.RL50.Ridge),mean(nonstat_data2$MSE.RL50.Lasso))
avgRL100 <- c(mean(nonstat_data2$MSE.RL100.Spat), mean(nonstat_data2$MSE.RL100.Spat2),mean(nonstat_data2$MSE.RL100.Max),mean(nonstat_data2$MSE.RL100.Bayes),mean(nonstat_data2$MSE.RL100.Ridge),mean(nonstat_data2$MSE.RL100.Lasso))
time <-c(mean(nonstat_data2$Time.Spat),mean(nonstat_data2$Time.Spat),mean(nonstat_data2$Time.Max),mean(nonstat_data2$Time.Bayes),mean(nonstat_data2$Time.Ridge),mean(nonstat_data2$Time.Lasso))
fail <- c(mean(nonstat_data2$Fail.Spat),mean(nonstat_data2$Fail.Spat),0,0,mean(nonstat_data2$Fail.Ridge),mean(nonstat_data2$Fail.Lasso)) 


GEV_avgMSE <- matrix(c(avgloc,
                       avgscale,
                       avgshape,
                       avgRL10,
                       avgRL20,
                       avgRL50,
                       avgRL100,
                       time,
                       fail), ncol = 6, byrow = TRUE)
colnames(GEV_avgMSE)<-c("Spat","Spat (Ridge Sites)","Max", "Bayes", "Ridge","Lasso")
rownames(GEV_avgMSE)<-c("Loc","Scale","Shape","RL10","RL20","RL50","RL100","Time","Fail")
GEV_MSE <- round(GEV_avgMSE,4)

sqrt(GEV_MSE)


##################################################################################
#Final boxplots with ggplot
##################################################################################
root.loc <- c(sqrt(nonstat_data2$MSE.Loc.Spat),sqrt(nonstat_data2$MSE.Loc.Max),sqrt(nonstat_data2$MSE.Loc.Ridge),sqrt(nonstat_data2$MSE.Loc.Lasso),sqrt(nonstat_data2$MSE.Loc.Bayes))
root.scale<- c(sqrt(nonstat_data2$MSE.Scale.Spat),sqrt(nonstat_data2$MSE.Scale.Max),sqrt(nonstat_data2$MSE.Scale.Ridge),sqrt(nonstat_data2$MSE.Scale.Lasso), sqrt(nonstat_data2$MSE.Scale.Bayes))
root.shape <- c(sqrt(nonstat_data2$MSE.Shape.Spat),sqrt(nonstat_data2$MSE.Shape.Max),sqrt(nonstat_data2$MSE.Shape.Ridge),sqrt(nonstat_data2$MSE.Shape.Lasso), sqrt(nonstat_data2$MSE.Shape.Bayes))
root.rl10 <-c(sqrt(nonstat_data2$MSE.RL10.Spat),sqrt(nonstat_data2$MSE.RL10.Max),sqrt(nonstat_data2$MSE.RL10.Ridge),sqrt(nonstat_data2$MSE.RL10.Lasso), sqrt(nonstat_data2$MSE.RL10.Bayes))
root.rl20 <-c(sqrt(nonstat_data2$MSE.RL20.Spat),sqrt(nonstat_data2$MSE.RL20.Max),sqrt(nonstat_data2$MSE.RL20.Ridge),sqrt(nonstat_data2$MSE.RL20.Lasso), sqrt(nonstat_data2$MSE.RL20.Bayes))
root.rl50 <- c(sqrt(nonstat_data2$MSE.RL50.Spat),sqrt(nonstat_data2$MSE.RL50.Max),sqrt(nonstat_data2$MSE.RL50.Ridge),sqrt(nonstat_data2$MSE.RL50.Lasso), sqrt(nonstat_data2$MSE.RL50.Bayes))
root.rl100 <- c(sqrt(nonstat_data2$MSE.RL100.Spat),sqrt(nonstat_data2$MSE.RL100.Max),sqrt(nonstat_data2$MSE.RL100.Ridge),sqrt(nonstat_data2$MSE.RL100.Lasso), sqrt(nonstat_data2$MSE.RL100.Bayes))
model <- c(rep("Spatial",length(nonstat_data2$MSE.Loc.Spat)),rep("Schlather",length(nonstat_data2$MSE.Loc.Max)),rep("Ridge",length(nonstat_data2$MSE.Loc.Ridge)),rep("Lasso",length(nonstat_data2$MSE.Loc.Lasso)),rep("Bayes",length(nonstat_data2$MSE.Loc.Bayes)))
gev_plot_data <- data.frame(RootLoc = root.loc, RootScale = root.scale, RootShape = root.shape, RootRL10 = root.rl10, RootRL20 = root.rl20, RootRL50=root.rl50, RootRL100=root.rl100, Model = model)


log.root.loc <- log(root.loc)
log.root.scale <- log(root.scale)
log.root.shape <- log(root.shape)
log.root.rl10 <- log(root.rl10)
log.root.rl20 <- log(root.rl20)
log.root.rl50 <- log(root.rl50)
log.root.rl100 <- log(root.rl100)
gev_plot_data <- data.frame(RootLoc = log.root.loc, RootScale = log.root.scale, RootShape = log.root.shape, RootRL10 = log.root.rl10, RootRL20 = log.root.rl20, RootRL50= log.root.rl50, RootRL100=log.root.rl100, Model = model)


#reorder data: Spatial 4, Bayes 1, Ridge 3, Lasso 2
gev_plot_data$Model = factor(gev_plot_data$Model, levels = levels(gev_plot_data$Model)[c(5,4,1,3,2)])

ggplot(gev_plot_data, aes(x=Model, y=RootLoc))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="Location", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_NonStationary_Loc.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootScale))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="Scale", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_NonStationary_Scale.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootShape))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="Shape", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_NonStationary_Shape.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootRL10))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="10-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_NonStationary_RL10.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootRL20))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="20-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_NonStationary_RL20.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootRL50))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="50-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_NonStationary_RL50.png"), width = 4, height = 6)

ggplot(gev_plot_data, aes(x=Model, y=RootRL100))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","orange2","palegreen4","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="100-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_NonStationary_RL100.png"), width = 4, height = 6)

##########################################################################################
##########################################################################################
#GEV Independent Study: Appendix E
GEV_Ind <- read.csv(paste0(dir.ind,"GEV_Independent_Results.csv"))
avg_ridge <- c(mean(GEV_Ind$MSE.Loc.Ridge.Ind), mean(GEV_Ind$MSE.Scale.Ridge.Ind), mean(GEV_Ind$MSE.Shape.Ridge.Ind), 
               mean(GEV_Ind$MSE.RL20.Ridge.Ind), mean(GEV_Ind$MSE.RL50.Ridge.Ind), mean(GEV_Ind$MSE.RL100.Ridge.Ind))
avg_lasso <- c(mean(GEV_Ind$MSE.Loc.Lasso.Ind), mean(GEV_Ind$MSE.Scale.Lasso.Ind), mean(GEV_Ind$MSE.Shape.Lasso.Ind), 
               mean(GEV_Ind$MSE.RL20.Lasso.Ind), mean(GEV_Ind$MSE.RL50.Lasso.Ind), mean(GEV_Ind$MSE.RL100.Lasso.Ind))

sqrt_ridge <- sqrt(avg_ridge)
sqrt_lasso <- sqrt(avg_lasso)

#percent change

ridge_chg <- c( (round(sqrt_ridge[1],3)/round(GEV_RMSE[1,4],3)-1)*100,
                (round(sqrt_ridge[2],3)/round(GEV_RMSE[2,4],3)-1)*100,
                (round(sqrt_ridge[3],4)/round(GEV_RMSE[3,4],4)-1)*100,
                (round(sqrt_ridge[4],3)/round(GEV_RMSE[5,4],3)-1)*100,
                (round(sqrt_ridge[5],3)/round(GEV_RMSE[6,4],3)-1)*100,
                (round(sqrt_ridge[6],3)/round(GEV_RMSE[7,4],3)-1)*100 )

spatial_chg_r <- c( (round(sqrt_ridge[1],3)/round(GEV_RMSE[1,1],3)-1)*100,
                (round(sqrt_ridge[2],3)/round(GEV_RMSE[2,1],3)-1)*100,
                (round(sqrt_ridge[3],4)/round(GEV_RMSE[3,1],4)-1)*100,
                (round(sqrt_ridge[4],3)/round(GEV_RMSE[5,1],3)-1)*100,
                (round(sqrt_ridge[5],3)/round(GEV_RMSE[6,1],3)-1)*100,
                (round(sqrt_ridge[6],3)/round(GEV_RMSE[7,1],3)-1)*100 )

lasso_chg <- c( (round(sqrt_lasso[1],3)/round(GEV_RMSE[1,5],3)-1)*100,
                (round(sqrt_lasso[2],3)/round(GEV_RMSE[2,5],3)-1)*100,
                (round(sqrt_lasso[3],4)/round(GEV_RMSE[3,5],4)-1)*100,
                (round(sqrt_lasso[4],3)/round(GEV_RMSE[5,5],3)-1)*100,
                (round(sqrt_lasso[5],3)/round(GEV_RMSE[6,5],3)-1)*100,
                (round(sqrt_lasso[6],3)/round(GEV_RMSE[7,5],3)-1)*100 )

spatial_chg_l <- c( (round(sqrt_lasso[1],3)/round(GEV_RMSE[1,1],3)-1)*100,
                    (round(sqrt_lasso[2],3)/round(GEV_RMSE[2,1],3)-1)*100,
                    (round(sqrt_lasso[3],4)/round(GEV_RMSE[3,1],4)-1)*100,
                    (round(sqrt_lasso[4],3)/round(GEV_RMSE[5,1],3)-1)*100,
                    (round(sqrt_lasso[5],3)/round(GEV_RMSE[6,1],3)-1)*100,
                    (round(sqrt_lasso[6],3)/round(GEV_RMSE[7,1],3)-1)*100 )

Table7 <- round(rbind(spatial_chg_r, ridge_chg, spatial_chg_l,lasso_chg),2)
colnames(Table7)<-c("Loc","Scale","Shape","RL20","RL50","RL100")
Table7
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##GPD Results
gpd_data = read.csv(paste0(dir.gpd,"GPD_Stationary_Results_1234.csv"))

avgscale <- c(mean(gpd_data$MSE.Scale.Spat),mean(gpd_data$MSE.Scale.Ridge),mean(gpd_data$MSE.Scale.Lasso))
avgshape <- c(mean(gpd_data$MSE.Shape.Spat),mean(gpd_data$MSE.Shape.Ridge),mean(gpd_data$MSE.Shape.Lasso))
avgRL20 <- c(mean(gpd_data$MSE.RL20.Spat),mean(gpd_data$MSE.RL20.Ridge),mean(gpd_data$MSE.RL20.Lasso))
avgRL50 <- c(mean(gpd_data$MSE.RL50.Spat),mean(gpd_data$MSE.RL50.Ridge),mean(gpd_data$MSE.RL50.Lasso))
avgRL100 <- c(mean(gpd_data$MSE.RL100.Spat),mean(gpd_data$MSE.RL100.Ridge),mean(gpd_data$MSE.RL100.Lasso))
time <- c(mean(gpd_data$Time.Spat),mean(gpd_data$Time.Ridge),mean(gpd_data$Time.Lasso))
fail <- c(mean(gpd_data$Fail.Spat),mean(gpd_data$Fail.Ridge),mean(gpd_data$Fail.Lasso)) 

GPD_avgMSE <- matrix(c(avgscale,
                       avgshape,
                       avgRL20,
                       avgRL50,
                       avgRL100,
                       time,
                       fail), ncol = 3, byrow = TRUE)
colnames(GPD_avgMSE)<-c("Spat" ,"Ridge","Lasso")
rownames(GPD_avgMSE)<-c("Scale","Shape","RL20","RL50","RL100","Time","Fail")
GPD_MSE <- round(GPD_avgMSE,4)

sqrt(GPD_MSE)

##################################################################################
#Final boxplots with ggplot
##################################################################################
root.scale<- c(sqrt(gpd_data$MSE.Scale.Spat),sqrt(gpd_data$MSE.Scale.Ridge),sqrt(gpd_data$MSE.Scale.Lasso))
root.shape <- c(sqrt(gpd_data$MSE.Shape.Spat),sqrt(gpd_data$MSE.Shape.Ridge),sqrt(gpd_data$MSE.Shape.Lasso))
root.rl20 <-c(sqrt(gpd_data$MSE.RL20.Spat),sqrt(gpd_data$MSE.RL20.Ridge),sqrt(gpd_data$MSE.RL20.Lasso))
root.rl50 <- c(sqrt(gpd_data$MSE.RL50.Spat),sqrt(gpd_data$MSE.RL50.Ridge),sqrt(gpd_data$MSE.RL50.Lasso))
root.rl100 <- c(sqrt(gpd_data$MSE.RL100.Spat),sqrt(gpd_data$MSE.RL100.Ridge),sqrt(gpd_data$MSE.RL100.Lasso))
model <- c(rep("Spatial",length(gpd_data$MSE.Scale.Spat)),rep("Ridge",length(gpd_data$MSE.Scale.Ridge)),rep("Lasso",length(gpd_data$MSE.Scale.Lasso)))
gpd_plot_data <- data.frame(RootScale = root.scale, RootShape = root.shape, RootRL20 = root.rl20, RootRL50=root.rl50, RootRL100=root.rl100, Model = model)

log.root.scale <- log(root.scale)
log.root.shape <- log(root.shape)
log.root.rl20 <- log(root.rl20)
log.root.rl50 <- log(root.rl50)
log.root.rl100 <- log(root.rl100)
gpd_plot_data <- data.frame(RootScale = log.root.scale, RootShape = log.root.shape, RootRL20 = log.root.rl20, RootRL50= log.root.rl50, RootRL100= log.root.rl100, Model = model)


#reorder data: Spatial 3, Ridge 2, Lasso 1
gpd_plot_data$Model = factor(gpd_plot_data$Model, levels = levels(gpd_plot_data$Model)[c(3,2,1)])

ggplot(gpd_plot_data, aes(x=Model, y=RootScale))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="Scale", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGPD_Scale.png"), width = 4, height = 6)

ggplot(gpd_plot_data, aes(x=Model, y=RootShape))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="Shape", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGPD_Shape.png"), width = 4, height = 6)

ggplot(gpd_plot_data, aes(x=Model, y=RootRL20))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="20-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGPD_RL20.png"), width = 4, height = 6)

ggplot(gpd_plot_data, aes(x=Model, y=RootRL50))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="50-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGPD_RL50.png"), width = 4, height = 6)


ggplot(gpd_plot_data, aes(x=Model, y=RootRL100))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="100-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGPD_RL100.png"), width = 4, height = 6)

########################################
########################################
#Nonstationary data

gpd_data = read.csv(paste0(dir.gpd2,"GPD_Nonstationary_Results_1234.csv"))

avgscale <- c(mean(gpd_data$MSE.Scale.Spat),mean(gpd_data$MSE.Scale.Ridge),mean(gpd_data$MSE.Scale.Lasso))
avgshape <- c(mean(gpd_data$MSE.Shape.Spat),mean(gpd_data$MSE.Shape.Ridge),mean(gpd_data$MSE.Shape.Lasso))
avgRL20 <- c(mean(gpd_data$MSE.RL20.Spat),mean(gpd_data$MSE.RL20.Ridge),mean(gpd_data$MSE.RL20.Lasso))
avgRL50 <- c(mean(gpd_data$MSE.RL50.Spat),mean(gpd_data$MSE.RL50.Ridge),mean(gpd_data$MSE.RL50.Lasso))
avgRL100 <- c(mean(gpd_data$MSE.RL100.Spat),mean(gpd_data$MSE.RL100.Ridge),mean(gpd_data$MSE.RL100.Lasso))

time <- c(mean(gpd_data$Time.Spat),mean(gpd_data$Time.Ridge),mean(gpd_data$Time.Lasso))
fail <- c(mean(gpd_data$Fail.Spat),mean(gpd_data$Fail.Ridge),mean(gpd_data$Fail.Lasso)) 

GPD_avgMSE <- matrix(c(avgscale,
                       avgshape,
                       avgRL20,
                       avgRL50,
                       avgRL100,
                       time,
                       fail), ncol = 3, byrow = TRUE)
colnames(GPD_avgMSE)<-c("Spat" ,"Ridge","Lasso")
rownames(GPD_avgMSE)<-c("Scale","Shape","RL20","RL50","RL100","Time","Fail")
GPD_MSE <- round(GPD_avgMSE,4)

GPD_MSE
sqrt(GPD_MSE)

##################################################################################
#Final boxplots with ggplot
##################################################################################
root.scale<- c(sqrt(gpd_data$MSE.Scale.Spat),sqrt(gpd_data$MSE.Scale.Ridge),sqrt(gpd_data$MSE.Scale.Lasso))
root.shape <- c(sqrt(gpd_data$MSE.Shape.Spat),sqrt(gpd_data$MSE.Shape.Ridge),sqrt(gpd_data$MSE.Shape.Lasso))
root.rl20 <-c(sqrt(gpd_data$MSE.RL20.Spat),sqrt(gpd_data$MSE.RL20.Ridge),sqrt(gpd_data$MSE.RL20.Lasso))
root.rl50 <- c(sqrt(gpd_data$MSE.RL50.Spat),sqrt(gpd_data$MSE.RL50.Ridge),sqrt(gpd_data$MSE.RL50.Lasso))
root.rl100 <- c(sqrt(gpd_data$MSE.RL100.Spat),sqrt(gpd_data$MSE.RL100.Ridge),sqrt(gpd_data$MSE.RL100.Lasso))
model <- c(rep("Spatial",length(gpd_data$MSE.Scale.Spat)),rep("Ridge",length(gpd_data$MSE.Scale.Ridge)),rep("Lasso",length(gpd_data$MSE.Scale.Lasso)))
#gpd_plot_data <- data.frame(RootScale = root.scale, RootShape = root.shape, RootRL20 = root.rl20, RootRL50=root.rl50, Model = model)

log.root.scale <- log(root.scale)
log.root.shape <- log(root.shape)
log.root.rl20 <- log(root.rl20)
log.root.rl50 <- log(root.rl50)
log.root.rl100 <- log(root.rl100)
gpd_plot_data <- data.frame(RootScale = log.root.scale, RootShape = log.root.shape, RootRL20 = log.root.rl20, RootRL50= log.root.rl50, RootRL100= log.root.rl100, Model = model)

#reorder data: Spatial 3, Ridge 2, Lasso 1
gpd_plot_data$Model = factor(gpd_plot_data$Model, levels = levels(gpd_plot_data$Model)[c(3,2,1)])

ggplot(gpd_plot_data, aes(x=Model, y=RootScale))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="Scale", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGPD_Nonstationary_Scale.png"), width = 4, height = 6)

ggplot(gpd_plot_data, aes(x=Model, y=RootShape))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="Shape", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGPD_Nonstationary_Shape.png"), width = 4, height = 6)

ggplot(gpd_plot_data, aes(x=Model, y=RootRL20))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="20-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGPD_Nonstationary_RL20.png"), width = 4, height = 6)

ggplot(gpd_plot_data, aes(x=Model, y=RootRL50))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="50-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGPD_Nonstationary_RL50.png"), width = 4, height = 6)

ggplot(gpd_plot_data, aes(x=Model, y=RootRL100))+
  theme_classic()+
  geom_boxplot(col = c("tomato2","steelblue3","orchid4"))+
  theme(plot.title = element_text(hjust=0.5, size = 20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20, angle =45, hjust = 1))+
  theme(axis.text.y = element_text(size=15))+
  labs(title="100-year Return Level", y="Logarithmic RMSE", x= NULL)
ggsave(paste0(plot.dir,"BoxPlot_SimGPD_Nonstationary_RL100.png"), width = 4, height = 6)



