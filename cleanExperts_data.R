
# Cleaning experts data
# H. Achicanoy
# CIAT, 2014

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Load packages & source code
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# source("F:/CIAT/CWR_project/Experts_evaluation/Scripts/downExperts_data.R")

# Scores data
# download_data(data_kind="scores",
#               final_dir="F:/CIAT/CWR_project/Experts_evaluation",
#               email="harold22010@gmail.com",
#               password="**********")

# Maps data
# download_data(data_kind="maps",
#               final_dir="F:/CIAT/CWR_project/Experts_evaluation",
#               email="harold22010@gmail.com",
#               password="**********")

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Read data
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

wk_dir = "F:/CIAT/CWR_project/Experts_evaluation"
path   = paste(wk_dir,"/Downloaded_data_(2014-04-28)",sep="")

# Data by taxa
scores_data = read.csv(paste(path,"/scoresEval.csv",sep=""),header=T)
maps_data   = read.csv(paste(path,"/mapsEval.csv",sep=""),header=T)

# Data by crop
scores_text = read.csv(paste(path,"/scores_textEval.csv",sep=""),header=T)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Cleaning data: step 1
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Necesary information
info_scores = scores_data[,c("Expert_id","Expert_nm","Crop_code")]
info_scores = unique(info_scores)
rownames(info_scores) = 1:nrow(info_scores)

# Merge the information
data_merge = merge(info_scores, maps_data, by=c("Expert_nm","Crop_code"))
rm(info_scores)

data_merge$Expert_id.y = NULL
names(data_merge)[match("Expert_id.x",names(data_merge))] <- "Expert_id"
data_merge = data_merge[,names(maps_data)]

# Include experts missing
expList = as.character(maps_data$Expert_nm)
expList = unique(expList)

exp_miss = setdiff(expList, data_merge$Expert_nm)
dat_miss = subset(maps_data,maps_data$Expert_nm==exp_miss)
rm(expList); rm(dat_miss)

data_merge = rbind(data_merge,dat_miss)
rownames(data_merge) = 1:nrow(data_merge)

crop_check = as.character(unique(data_merge$Crop_code[data_merge$Expert_nm==exp_miss]))

if(crop_check %in% scores_data$Crop_code){
  id = as.character(unique(scores_data$Expert_id[scores_data$Crop_code==crop_check]))
} else {
  cat("The crop has an unique expert \n")
}

if(!is.null(id)){
  data_merge$Expert_id[data_merge$Crop_code==crop_check & data_merge$Expert_nm==exp_miss] <- "Expert_2"
}

maps_data = data_merge[order(data_merge$Crop_code),]
rownames(maps_data) = 1:nrow(maps_data)
rm(data_merge)

rm(crop_check)
rm(exp_miss)
rm(id)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Cleaning data: step 2
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

dataList <- list(maps_data,scores_data,scores_text)

clean_data <- function(data, ...){
  
  data = as.matrix(data)
  data[data=="N/A"] <- NA
  data = as.data.frame(data)
  
  return(data)
  
}

dataList <- lapply(dataList, clean_data)
names(dataList) <- c("maps_data","scores_data","scores_text")

maps_data   <- dataList[["maps_data"]]
scores_data <- dataList[["scores_data"]]
scores_text <- dataList[["scores_text"]]

rm(dataList)

path <- paste(path,"/clean",sep="")
if(!file.exists(path)){dir.create(path)}

write.csv(maps_data, paste(path,"/maps_corrected.csv",sep=""), row.names=FALSE)
write.csv(scores_data, paste(path,"/scores_corrected.csv",sep=""), row.names=FALSE)
write.csv(scores_text, paste(path,"/sc_text_corrected.csv",sep=""), row.names=FALSE)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Cleaning data: step 3
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

wk_dir = "F:/CIAT/CWR_project/Experts_evaluation"
path   = paste(wk_dir,"/Downloaded_data_(2014-04-28)/clean",sep="")

maps_data = read.csv(paste(path,"/maps_corrected.csv",sep=""),header=T)
scores_data = read.csv(paste(path,"/scores_corrected.csv",sep=""),header=T)
scores_text = read.csv(paste(path,"/sc_text_corrected.csv",sep=""),header=T)
fps_data    = read.csv(paste(wk_dir,"/allPriorities2014-04-21.csv",sep=""),header=T)

### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
### Delete words: var. and subsp.
### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###

### ========================================== ###
# Maps data: var.
tNames = unique(maps_data$Taxon)
tNames = as.vector(tNames)
fout   = tNames[grep("var._",tNames,fixed=T)]
maps_data$Taxon = as.vector(maps_data$Taxon)
for(f in fout) {
  true_name  = strsplit(paste(f),"var._")[[1]]
  true_name3 = paste(true_name[1],true_name[2],sep="")
  maps_data$Taxon[which(maps_data$Taxon==f)] = array(true_name3,dim=length(which(maps_data$Taxon==f)))
}; rm(f)

# Maps data: subsp.
tNames = unique(maps_data$Taxon)
tNames = as.vector(tNames)
fout   = tNames[grep("subsp.",tNames,fixed=T)]
maps_data$Taxon = as.vector(maps_data$Taxon)
for(f in fout) {
  true_name  = strsplit(paste(f),"subsp._")[[1]]
  true_name3 = paste(true_name[1],true_name[2],sep="")
  maps_data$Taxon[which(maps_data$Taxon==f)] = array(true_name3,dim=length(which(maps_data$Taxon==f)))
}; rm(f)

### ========================================== ###
# Scores data: var.
tNames = unique(scores_data$Taxon)
tNames = as.vector(tNames)
fout   = tNames[grep("var._",tNames,fixed=T)]
scores_data$Taxon = as.vector(scores_data$Taxon)
for(f in fout) {
  true_name  = strsplit(paste(f),"var._")[[1]]
  true_name3 = paste(true_name[1],true_name[2],sep="")
  scores_data$Taxon[which(scores_data$Taxon==f)] = array(true_name3,dim=length(which(scores_data$Taxon==f)))
}; rm(f)

# Scores data: subsp.
tNames = unique(scores_data$Taxon)
tNames = as.vector(tNames)
fout   = tNames[grep("subsp.",tNames,fixed=T)]
scores_data$Taxon = as.vector(scores_data$Taxon)
for(f in fout) {
  true_name  = strsplit(paste(f),"subsp._")[[1]]
  true_name3 = paste(true_name[1],true_name[2],sep="")
  scores_data$Taxon[which(scores_data$Taxon==f)] = array(true_name3,dim=length(which(scores_data$Taxon==f)))
}; rm(f)

### ========================================== ###
# FPS data: var.
tNames = unique(fps_data$TAXON)
tNames = as.vector(tNames)
fout   = tNames[grep("var._",tNames,fixed=T)]
fps_data$TAXON = as.vector(fps_data$TAXON)
for(f in fout) {
  true_name  = strsplit(paste(f),"var._")[[1]]
  true_name3 = paste(true_name[1],true_name[2],sep="")
  fps_data$TAXON[which(fps_data$TAXON==f)] = array(true_name3,dim=length(which(fps_data$TAXON==f)))
}; rm(f)

# FPS data: subsp.
tNames = unique(fps_data$TAXON)
tNames = as.vector(tNames)
fout   = tNames[grep("subsp.",tNames,fixed=T)]
fps_data$TAXON = as.vector(fps_data$TAXON)
for(f in fout) {
  true_name  = strsplit(paste(f),"subsp._")[[1]]
  true_name3 = paste(true_name[1],true_name[2],sep="")
  fps_data$TAXON[which(fps_data$TAXON==f)] = array(true_name3,dim=length(which(fps_data$TAXON==f)))
}; rm(f)

### ========================================== ###
# Write data
write.csv(maps_data,"F:/CIAT/CWR_project/Experts_evaluation/Downloaded_data_(2014-04-28)/clean/maps_corrected.csv",row.names=F)
write.csv(scores_data,"F:/CIAT/CWR_project/Experts_evaluation/Downloaded_data_(2014-04-28)/clean/scores_corrected.csv",row.names=F)
write.csv(fps_data,"F:/CIAT/CWR_project/Experts_evaluation/Downloaded_data_(2014-04-28)/clean/fps_corrected.csv",row.names=F)

### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
