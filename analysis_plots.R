# Analysis of expert's evaluation throught multivariate methods
# H. Achicanoy
# CIAT, 2014

# ========================================================= #
# Load packages & source code
# ========================================================= #

library(FactoMineR)
library(ggplot2)
library(reshape)
library(ForImp)
library(aspect)
library(GGally)
library(plyr)

# ========================================================= #
# Analyzing results
# ========================================================= #

plot_dir = paste(wk_dir,"/Results/Plots",sep="")
if(!file.exists(plot_dir)){dir.create(plot_dir)}

results = results[-1]
lapply(results, class)
id_res = unlist(lapply(results, function(x){class(x)[1]}))
nm_mfa = names(id_res)[which(id_res=="MFA")]

# ========================================================= #
# Cumulative percentage of variance for each crop
# MFA analysis
# ========================================================= #

pc_mfa = lapply(nm_mfa, function(x){results[[x]]$eig[,3]})
names(pc_mfa) = nm_mfa

pc_mfa = ldply(pc_mfa, cbind)
colnames(pc_mfa) = c("Crop","Percentage")
pc_mfa$Dimensions = c(seq(1,sum(pc_mfa$Crop=="avena")),
                      seq(1,sum(pc_mfa$Crop=="beet")),
                      seq(1,sum(pc_mfa$Crop=="cajanus")),
                      seq(1,sum(pc_mfa$Crop=="helianthus")),
                      seq(1,sum(pc_mfa$Crop=="wheat")))

require(ggplot2)
p <- ggplot(pc_mfa, aes(x=Dimensions, y=Percentage, colour=Crop)) + geom_line(size=1) + ylim(0,100)
p <- p + theme_bw()
p <- p + ggtitle("Cumulative variance by crop, MFA analysis")
p <- p + scale_colour_discrete(labels=c("Avena", "Beet", "Cajanus", "Helianthus", "Wheat"))
p <- p + theme(plot.title=element_text(face="bold")) + geom_vline(xintercept=2, colour="red", linetype = "longdash")
p <- p + xlab("Number of dimensions") + ylab("Cumulative percentage of variance (%)")
p <- p + scale_x_continuous(breaks=1:14)
p

ggsave(filename=paste(plot_dir,"/cumulative_variance_MFA.png",sep=""),plot=p,width=20,height=17,units="cm")

# ========================================================= #
# Cumulative percentage of variance for each crop
# MFA and PCA analysis (1)
# ========================================================= #

nm_all = names(id_res)[which(id_res=="MFA"|id_res=="PCA")]
pc_all = lapply(nm_all, function(x){results[[x]]$eig[,3]})
names(pc_all) = nm_all

pc_all = ldply(pc_all, cbind)
colnames(pc_all) = c("Crop","Percentage")
pc_all$Dimensions = c(seq(1,sum(pc_all$Crop=="avena")),
                      seq(1,sum(pc_all$Crop=="beet")),
                      seq(1,sum(pc_all$Crop=="cajanus")),
                      seq(1,sum(pc_all$Crop=="daucus")),
                      seq(1,sum(pc_all$Crop=="dioscorea")),
                      seq(1,sum(pc_all$Crop=="helianthus")),
                      seq(1,sum(pc_all$Crop=="lettuce")),
                      seq(1,sum(pc_all$Crop=="malus")),
                      seq(1,sum(pc_all$Crop=="quinoa")),
                      seq(1,sum(pc_all$Crop=="rice")),
                      seq(1,sum(pc_all$Crop=="safflower")),
                      seq(1,sum(pc_all$Crop=="wheat")))

require(ggplot2)
p <- ggplot(pc_all, aes(x=Dimensions, y=Percentage, colour=Crop)) + geom_line(size=1) + ylim(0,100)
p <- p + theme_bw()
p <- p + ggtitle("Cumulative variance by crop, MFA and PCA analysis")
p <- p + scale_colour_discrete(labels=c("Avena", "Beet", "Cajanus", "Daucus", "Dioscorea", "Helianthus", "Lettuce", "Malus", "Quinoa", "Rice", "Safflower", "Wheat"))
p <- p + theme(plot.title=element_text(face="bold")) + geom_vline(xintercept=2, colour="red", linetype = "longdash")
p <- p + xlab("Number of dimensions") + ylab("Cumulative percentage of variance (%)")
p <- p + scale_x_continuous(breaks=1:14)
p

ggsave(filename=paste(plot_dir,"/cumulative_variance_MFA_PCA.png",sep=""),plot=p,width=20,height=17,units="cm")

# ========================================================= #
# Cumulative percentage of variance for each crop
# MFA and PCA analysis (2)
# ========================================================= #

id_df = data.frame(Crop=names(id_res),Analysis=id_res)
rownames(id_df) = 1:nrow(id_df)
id_df = id_df[which(id_df$Analysis=="MFA"|id_df$Analysis=="PCA"),]
rownames(id_df) = 1:nrow(id_df)
pc_all$Analysis = 0

for(i in 1:nrow(id_df)){
  pc_all$Analysis[pc_all$Crop==id_df$Crop[i]] <- as.character(paste(id_df$Analysis[i]))
}
rm(i); rm(id_df)

require(ggplot2)
p <- ggplot(pc_all, aes(x=Dimensions, y=Percentage, shape=Crop, colour=Analysis)) + geom_line(size=1) + ylim(0,100)
p <- p + theme_bw()
p <- p + ggtitle("Cumulative variance by crop, MFA and PCA analysis")
p <- p + theme(plot.title=element_text(face="bold")) + geom_vline(xintercept=2, colour="red", linetype = "longdash")
p <- p + xlab("Number of dimensions") + ylab("Cumulative percentage of variance (%)")
p <- p + scale_x_continuous(breaks=1:14)
p

ggsave(filename=paste(plot_dir,"/cumulative_variance_MFA_PCA2.png",sep=""),plot=p,width=20,height=17,units="cm")

# ========================================================= #
# Analysis of components for crop
# ========================================================= #

eval_crop <- function(analysis, crop){
  
  crop    <- crop
  a_class <- class(analysis)[1]
  
  if(a_class=="MFA"){
    
    plot(analysis, choix="ind",
         title=paste("Individual factor map for ",crop,sep=""))
    #plot(analysis, choix="ind", partial="all")
    plot(analysis, choix="var", habillage="group",
         title=paste("Correlation circle of variables for ",crop,sep=""))
    
  } else {
    if(a_class=="PCA"){
      
      plot(analysis, choix="ind",
           title=paste("Individual factor map for ",crop,sep=""))
      plot(analysis, choix="var",
           title=paste("Correlation circle of variables for ",crop,sep=""))
      
    } else {
      cat("No results available\n")
    }
  }
  return("Done!")
}
mapply(results, names(results), FUN=function(x,y){eval_crop(x,y)})

# eval_crop for cajanus
# eval_crop(results[["cajanus"]],"cajanus")

# ========================================================= #
# Evaluation index
# ========================================================= #

index_cal  <- function(analysis){
  
  analysis <- analysis
  a_class  <- class(analysis)[1]
  
  if(a_class=="MFA"){
    
    w <- analysis$eig[,1] # Extracción de valores propios
    q <- w/sum(w)         # Porcentaje de varianza explicada
    i <- q[1]*analysis$quanti.var$coord[,1] + q[2]*analysis$quanti.var$coord[,2] # Ponderaciones de variables
    l <- i/sum(i)         # Estandarización
    Ind0 <- scale(analysis$call$XTDC)%*%l   # Ponderación individuos
    Ind1 <- (exp(Ind0))/(1+(exp(Ind0)))*100 # Índice final
    
  } else {
    
    if(a_class=="PCA"){
      
      w <- analysis$eig[,1] # Extracción de valores propios
      q <- w/sum(w)         # Porcentaje de varianza explicada
      i <- q[1]*analysis$var$coord[,1] + q[2]*analysis$var$coord[,2] # Ponderaciones de variables
      l <- i/sum(i)         # Estandarización
      Ind0 <- scale(analysis$call$X)%*%l      # Ponderación individuos
      Ind1 <- (exp(Ind0))/(1+(exp(Ind0)))*100 # Índice final
      
    } else {
      
      cat("The crop tal hasn't appropiated data for index construction \n")
      
    }
  
  }
  
}
index_res  <- lapply(results, function(x){index_cal(x)})
index_res1 <- unlist(index_res)

### Evaluation index

summary(index_res1)

# Histogram
colscale <- colorRampPalette(c("red","yellow","forestgreen"), space="rgb")(100)
hist(index_res1,prob=T,xlab="Evaluation index",col="gray90",
     main="Evaluation index distribution",border="gray80")
box()

# Density plot
dens <- density(index_res1)
plot(dens, type="n", xlim=c(0,100),
     xlab="Evaluation index",
     main="Evaluation index density distribution")
polygon(dens$x, dens$y, col="gray90", border="gray80")
box()

# Taxa ordered
plot(sort(index_res1),
     xlab="Ordered taxa",
     ylab="Evaluation index",
     main="Ordered taxa vs evaluation index",
     type="p", pch=21, cex=2,
     col="gray80", bg="gray90")
abline(h=50,col=2,lty=2)
# 107 taxones con valoraciones superiores a 50 puntos en el índice
# ~60% de los taxones evaluados

# Taxa ordered gradient colors
colscale <- colorRampPalette(c("red","yellow","forestgreen"), space="rgb")(100)
colscale <- data.frame(index=1:100,colors=colscale)

index_sort <- sort(index_res1)
auxiliar   <- round(index_sort)
col.scale  <- unlist(lapply(auxiliar,function(x){ colscale$colors[match(x,table=colscale$index)] }))
col.scale  <- as.character(col.scale)

plot(index_sort,
     xlab="Ordered taxa",
     ylab="Evaluation index",
     main="Ordered taxa vs evaluation index",
     type="p", pch=21, cex=2,
     col=col.scale,bg=col.scale) # "gray90"
abline(h=50,col=2,lty=2)
# 107 taxones con valoraciones superiores a 50 puntos en el índice
# ~60% de los taxones evaluados

# Helianthus case

helianthus <- as.data.frame(index_res[["helianthus"]])
helianthus <- cbind(helianthus,rownames(helianthus))
rownames(helianthus) <- 1:nrow(helianthus)
colnames(helianthus) <- c("index","taxon")
helianthus <- helianthus[order(helianthus$index),]
rownames(helianthus) <- 1:nrow(helianthus)

# Color scale
colscale <- colorRampPalette(c("red","yellow","forestgreen"), space="rgb")(100)
colscale <- data.frame(index=1:100,colors=colscale)

auxiliar   <- round(helianthus$index)
col.scale  <- unlist(lapply(auxiliar,function(x){ colscale$colors[match(x,table=colscale$index)] }))
col.scale  <- as.character(col.scale)

# Evaluation index for helianthus
par(mar=c(5, 18, 4, 2) + 0.1)
barplot(helianthus$index, horiz=T, names.arg=helianthus$taxon,
        las=1, col=col.scale, border=NA, xlim=c(0,100),
        xlab="Evaluation index", main="Evaluation index for helianthus taxa")
box()

# Lettuce case

lettuce <- as.data.frame(index_res[["lettuce"]])
lettuce <- cbind(lettuce,rownames(lettuce))
rownames(lettuce) <- 1:nrow(lettuce)
colnames(lettuce) <- c("index","taxon")
lettuce <- lettuce[order(lettuce$index),]
rownames(lettuce) <- 1:nrow(lettuce)

# Color scale
colscale <- colorRampPalette(c("red","yellow","forestgreen"), space="rgb")(100)
colscale <- data.frame(index=1:100,colors=colscale)

auxiliar   <- round(lettuce$index)
col.scale  <- unlist(lapply(auxiliar,function(x){ colscale$colors[match(x,table=colscale$index)] }))
col.scale  <- as.character(col.scale)

# Evaluation index for lettuce
par(mar=c(5, 12, 4, 2) + 0.1)
barplot(lettuce$index, horiz=T, names.arg=lettuce$taxon,
        las=1, col=col.scale, border=NA, xlim=c(0,100),
        xlab="Evaluation index", main="Evaluation index for lettuce taxa")
box()

# Cajanus case

cajanus <- as.data.frame(index_res[["cajanus"]])
cajanus <- cbind(cajanus,rownames(cajanus))
rownames(cajanus) <- 1:nrow(cajanus)
colnames(cajanus) <- c("index","taxon")
cajanus <- cajanus[order(cajanus$index),]
rownames(cajanus) <- 1:nrow(cajanus)

# Color scale
colscale <- colorRampPalette(c("red","yellow","forestgreen"), space="rgb")(100)
colscale <- data.frame(index=1:100,colors=colscale)

auxiliar   <- round(cajanus$index)
col.scale  <- unlist(lapply(auxiliar,function(x){ colscale$colors[match(x,table=colscale$index)] }))
col.scale  <- as.character(col.scale)

cajanus$Label <- gsub(pattern="ajanus_",replacement=". ",x=as.character(cajanus$taxon))

# Evaluation index for cajanus
par(mar=c(5, 10, 4, 2) + 0.1)
x <- barplot(cajanus$index, horiz=T, names.arg=cajanus$Label)
barplot(cajanus$index, horiz=T, names.arg=cajanus$Label,
        las=1, col=col.scale, border=NA, xlim=c(0,100),font.axis=1,
        xlab="Evaluation index", main="Evaluation index for cajanus taxa",yaxt="n")
box()
axis(2,x[,1],as.character(cajanus$Label),font=3,las=2)

par(mar=c(5, 10, 4, 2) + 0.1)
plot(cajanus$index,1:length(cajanus$index),xlim=c(0,100),ylim=c(0,16),
     xlab="Evaluation index", main="Evaluation index for cajanus taxa",
     ylab="",yaxt="n",col=col.scale,pch=20,cex=2)
axis(2,1:15,as.character(cajanus$Label),font=3,las=2)



# Sería interesante explorar la correlación entre la diferencia relativa
# de cada uno de los taxones y su correspondiente índice de evaluación
# según los resultados individuales

# ========================================================= #



















# Final priority score and expert priority score in the same plot
# with ggplot2

dir <- "F:/CIAT/CWR_project/Papers"

library(ggplot2)

fps_cajanus <- read.csv(paste(dir,"/cajanus/priorities_cajanus.csv",sep=""),header=T)
fps_cajanus


# Vertical way
p <- ggplot(fps_cajanus, aes(x=FPS, y=reorder(TAXON,FPS),colour="black")) + theme_bw()
p <- p + xlim(0,10) + geom_blank()
p <- p + annotate("rect", xmin=0, xmax=3, ymin=0, ymax=17, alpha=.3, fill="red")
p <- p + annotate("rect", xmin=3, xmax=5, ymin=0, ymax=17, alpha=.5, fill="orange")
p <- p + annotate("rect", xmin=5, xmax=7.5, ymin=0, ymax=17, alpha=.5, fill="yellow")
p <- p + annotate("rect", xmin=7.5, xmax=10, ymin=0, ymax=17, alpha=.5, fill="forestgreen")
p <- p + xlab(label="Priority score") + ylab("")
p <- p + scale_x_continuous(breaks=0:10)
p <- p + geom_point(stat="identity",size=7)
p <- p + geom_point(data=fps_cajanus, aes(x=SRS, y=reorder(TAXON,FPS), colour="blue3"),size=3.5)
p <- p + geom_point(data=fps_cajanus, aes(x=GRS, y=reorder(TAXON,FPS), colour="firebrick1"),size=3.5)
p <- p + geom_point(data=fps_cajanus, aes(x=ERS, y=reorder(TAXON,FPS), colour="forestgreen"),size=3.5)
p <- p + theme(panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank(),
               panel.grid.major.y = element_line(colour="grey60", linetype="dashed"))
p <- p + theme(axis.text.x=element_text(size=15),
               axis.text.y=element_text(face="italic",size=15),
               axis.title.x=element_text(face="bold",size=15),
               axis.title.y=element_text(face="bold",size=15))
cols <- c("black","blue3","firebrick1","forestgreen")
p <- p + scale_colour_manual(name="",labels=c("FPS","SRS","GRS","ERS"),values=cols, guide=guide_legend(fill=NULL))
p <- p + theme(legend.position="top")
p

# Cargar información de expertos

scores_data[scores_data$Crop_code=="Cajanus",]
dim(scores_data[scores_data$Crop_code=="Cajanus",])










