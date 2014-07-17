# Analysis of expert's evaluation throught multivariate methods
# H. Achicanoy
# CIAT, 2014

# ---------------------------------------------------------------------------- #
# Load packages & source code
# ---------------------------------------------------------------------------- #

g <- gc()
rm(list=ls())

library(FactoMineR)
library(ggplot2)
library(reshape)
library(ForImp)
library(aspect)
library(GGally)
library(plyr)

# source("F:/CIAT/CWR_project/Experts_evaluation/Scripts/cleanExperts_data.R")


# ---------------------------------------------------------------------------- #
# Define paths
# ---------------------------------------------------------------------------- #

wk_dir = "F:/CIAT/CWR_project/Experts_evaluation"
path   = paste(wk_dir,"/Downloaded_data_(2014-04-28)/clean",sep="")


# ---------------------------------------------------------------------------- #
# Read data
# ---------------------------------------------------------------------------- #

maps_data   = read.csv(paste(path,"/maps_corrected.csv",sep=""),header=T)
scores_data = read.csv(paste(path,"/scores_corrected.csv",sep=""),header=T)
scores_text = read.csv(paste(path,"/sc_text_corrected.csv",sep=""),header=T)
fps_data    = read.csv(paste(wk_dir,"/allpriorities_2014-04-21.csv",sep=""),header=T)
fps_modf    = read.csv(paste(wk_dir,"/modpriorities_2014-05-15.csv",sep=""),header=T)


# ---------------------------------------------------------------------------- #
# Merge data: SCORES and MAPS data
# ---------------------------------------------------------------------------- #

data_merge = merge(maps_data,scores_data,by=c("Expert_id","Expert_nm","Crop_code","Taxon","Priority_crop"))
data_merge$Position.x = NULL
names(data_merge)[match("Position.y",names(data_merge))] <- "Position"
names.col  = c("Expert_id","Expert_nm","Priority_crop","Crop_code","Taxon","Position","Comparable","Contextual","Evaluation","Occ_data","SDM_map","Gap_map")
data_merge = data_merge[,names.col]; rm(names.col)
data_merge = data_merge[order(data_merge$Crop_code),]
rownames(data_merge) = 1:nrow(data_merge)


# ---------------------------------------------------------------------------- #
# FPS data
# ---------------------------------------------------------------------------- #

fps_data = fps_data[,c("CROP_CODE","TAXON","SRS","GRS","ERS","FPS","FPCAT")]
#ggpairs(fps_data[,c("SRS","GRS","ERS","FPS")],alpha=0.4)
#ggpairs(fps_modf[,c("SRS","GRS","ERS","FPS")],alpha=0.4)


# ---------------------------------------------------------------------------- #
# Function to run Multiple Factor Analysis for expert's evaluation data
# ---------------------------------------------------------------------------- #

# List of crops
cropList = sort(unique(as.character(data_merge$Crop_code)))

# Function inputs
# 1. Data set
# 2. List of crops

mfa_fun <- function(data, list){
  
  cropList   <- list
  data_merge <- data
  results    <- list(0)
  
  for(i in 1:length(cropList)){
    
    cat("\n =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n")
    cat(" Processing:",tolower(cropList[i]))
    cat("\n =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n")
    cat("\n")
    
    # Organize data horizontal way
    options(warn=-1)
    x = subset(data_merge, data_merge$Crop_code==cropList[i])
    x$Position = NULL
    x$Expert_nm = NULL
    x$Crop_code = NULL
    x$Priority_crop = NULL
    x = reshape(x, timevar="Expert_id", idvar="Taxon", direction="wide")
    
    # Identify Final Priority Score for each taxon
    if(cropList[i] %in% c("beet","dioscorea","groundnut","lettuce","quinoa","safflower","soy","spinach","tomato")){
      taxList = as.character(x$Taxon)
      taxID   = lapply(taxList, function(x){grep(pattern=x,fps_modf$TAXON,fixed=T,value=F)})
      names(taxID) = taxList
      id_null = unlist(lapply(taxID, length))
      if(any(id_null==0)){
        id_null = which(unlist(lapply(taxID, length))==0)
        fps = unlist(lapply(taxID, function(x){min(fps_modf$FPS[x],na.rm=T)}))
        fps[id_null] = 0
      } else {
        fps = unlist(lapply(taxID, function(x){min(fps_modf$FPS[x],na.rm=T)}))
      }
      rm(taxList); rm(taxID); rm(id_null)
    } else {
      taxList = as.character(x$Taxon)
      taxID   = lapply(taxList, function(x){grep(pattern=x,fps_data$TAXON,fixed=T,value=F)})
      names(taxID) = taxList
      id_null = unlist(lapply(taxID, length))
      if(any(id_null==0)){
        id_null = which(unlist(lapply(taxID, length))==0)
        fps = unlist(lapply(taxID, function(x){min(fps_data$FPS[x],na.rm=T)}))
        fps[id_null] = 0
      } else {
        fps = unlist(lapply(taxID, function(x){min(fps_data$FPS[x],na.rm=T)}))
      }
      rm(taxList); rm(taxID); rm(id_null)
    }
    
    # Add taxon name to rows of data frame
    rownames(x) = x$Taxon
    x$Taxon = NULL
    
    # Identify number of experts that responded the surveys
    n_exp = names(x)
    n_exp = unlist(lapply(n_exp, function(x){strsplit(x, split="_")}))
    n_chr = unlist(lapply(n_exp, function(x){nchar(x)}))
    n_exp = as.numeric(n_exp[which(n_chr==1)]); rm(n_chr)
    n_exp = unique(n_exp)
    
    cat(" Calculate relative difference \n")
    for(ex in n_exp){
      label = paste("x$Comparable.Expert_",ex," = 1-abs((fps-x$Comparable.Expert_",ex,")/10)",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Comparable.Expert_",ex," = as.numeric(x$Comparable.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Contextual.Expert_",ex," = 1-abs((fps-x$Contextual.Expert_",ex,")/10)",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Contextual.Expert_",ex," = as.numeric(x$Contextual.Expert_",ex,")",sep="")
      eval(parse(text=label))
      label = paste("x$Evaluation.Expert_",ex," = as.character(x$Evaluation.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Occ_data.Expert_",ex," = as.character(x$Occ_data.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$SDM_map.Expert_",ex," = as.character(x$SDM_map.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Gap_map.Expert_",ex," = as.character(x$Gap_map.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
    }; rm(ex)
    
    # Exclude for analysis datasets with less than 6 taxa
    if(nrow(x) > 6){
      
      cat(" Filtering data \n")
      
      filter_fun <- function(X){
        
        na_var = apply(X, 2, function(x){sum(is.na(x))})
        mtch   = names(na_var[which(na_var!=nrow(X))]); rm(na_var)
        X      = X[,mtch]; rm(mtch)
        
        lapply(names(X), function(x){
          label  = paste("sd(as.numeric(as.factor(X[,","'",x,"'","])),na.rm=T)",sep="")
          sd_val = eval(parse(text=label)); rm(label)
          if(sd_val == 0){
            cat(" ",x,"variable hasn't associated variability. Thus eliminiating. \n")
            label = paste("X$",x," = NULL",sep="")
            eval(parse(text=label)); rm(label)
          } else{
            cat(" ",x,"variable is ok for the analysis \n")
          }
        })
        
        return(X)
        
      }
      x = filter_fun(x)
      
      cat(" Verify if exists missing data \n")
      
      imp_fun = function(X, ...){
        
        if(any(is.na(X))){
          
          cat(" The dataset has missing data: imputing. \n")
          
          # Step 1: convert to numeric results only for categorical data
          X = as.matrix(X)
          X = gsub(pattern="[Ss][Tt][Rr][Oo][Nn][Gg][Ll][Yy] [Dd][Ii][Ss][Aa][Gg][Rr][Ee][Ee]", replacement=1, X, perl=T)
          X = gsub(pattern="[Ss][Tt][Rr][Oo][Nn][Gg][Ll][Yy] [Aa][Gg][Rr][Ee][Ee]", replacement=5, X, perl=T)
          X = gsub(pattern="[Dd][Ii][Ss][Aa][Gg][Rr][Ee][Ee]", replacement=2, X, perl=T)
          X = gsub(pattern="[Nn][Ee][Uu][Tt][Rr][Aa][Ll]", replacement=3, X, perl=T)
          X = gsub(pattern="[Aa][Gg][Rr][Ee][Ee]", replacement=4, X, perl=T)
          
          y = matrix(as.numeric(unlist(X)),nrow=nrow(X))
          rownames(y) = rownames(X)
          colnames(y) = colnames(X)
          
          X = matrix(as.numeric(unlist(X)),nrow=nrow(y))
          rownames(X) = rownames(y)
          colnames(X) = colnames(y)
          
          # Step 2: data imputation only for categorical data
          
          cat_data = c("Evaluation.","Occ_data.","SDM_map.","Gap_map.")
          mtch = sort(unlist(lapply(cat_data,function(x){grep(pattern=x,colnames(y))})))
          y = y[,mtch]; rm(mtch)
          
          y = modeimp(y)
          X[,colnames(y)] = y[,colnames(y)]
          
          X[which(X=="1")] = "Strongly Disagree"
          X[which(X=="2")] = "Disagree"
          X[which(X=="3")] = "Neutral"
          X[which(X=="4")] = "Agree"
          X[which(X=="5")] = "Strongly Agree"
          
          X = as.data.frame(X)
          
          for(exp in n_exp){
            
            if(paste("Comparable.Expert_",exp,sep="") %in% colnames(X)){
              label = paste("X$Comparable.Expert_",exp," = as.character(X$Comparable.Expert_",exp,")",sep="")
              eval(parse(text=label)); rm(label)
              label = paste("X$Comparable.Expert_",exp," = as.numeric(X$Comparable.Expert_",exp,")",sep="")
              eval(parse(text=label)); rm(label)
            }
            
            if(paste("Contextual.Expert_",exp,sep="") %in% colnames(X)){
              label = paste("X$Contextual.Expert_",exp," = as.character(X$Contextual.Expert_",exp,")",sep="")
              eval(parse(text=label))
              label = paste("X$Contextual.Expert_",exp," = as.numeric(X$Contextual.Expert_",exp,")",sep="")
              eval(parse(text=label))
            }
            
          }; rm(exp); rm(label)
          rm(y)
          
          return(X)
          
        } else {
          
          X = as.matrix(X)
          X = gsub(pattern="[Ss][Tt][Rr][Oo][Nn][Gg][Ll][Yy] [Dd][Ii][Ss][Aa][Gg][Rr][Ee][Ee]", replacement="1", X, useBytes=TRUE)
          X = gsub(pattern="[Ss][Tt][Rr][Oo][Nn][Gg][Ll][Yy] [Aa][Gg][Rr][Ee][Ee]", replacement="5", X, useBytes=TRUE)
          X = gsub(pattern="[Dd][Ii][Ss][Aa][Gg][Rr][Ee][Ee]", replacement="2", X, useBytes=TRUE)
          X = gsub(pattern="[Nn][Ee][Uu][Tt][Rr][Aa][Ll]", replacement="3", X, useBytes=TRUE)
          X = gsub(pattern="[Aa][Gg][Rr][Ee][Ee]", replacement="4", X, useBytes=TRUE)
          X[which(X=="1")] = "Strongly Disagree"
          X[which(X=="2")] = "Disagree"
          X[which(X=="3")] = "Neutral"
          X[which(X=="4")] = "Agree"
          X[which(X=="5")] = "Strongly Agree"
          X = as.data.frame(X)
          
          for(exp in n_exp){
            
            if(paste("Comparable.Expert_",exp,sep="") %in% colnames(X)){
              label = paste("X$Comparable.Expert_",exp," = as.character(X$Comparable.Expert_",exp,")",sep="")
              eval(parse(text=label)); rm(label)
              label = paste("X$Comparable.Expert_",exp," = as.numeric(X$Comparable.Expert_",exp,")",sep="")
              eval(parse(text=label)); rm(label)
            }
            
            if(paste("Contextual.Expert_",exp,sep="") %in% colnames(X)){
              label = paste("X$Contextual.Expert_",exp," = as.character(X$Contextual.Expert_",exp,")",sep="")
              eval(parse(text=label))
              label = paste("X$Contextual.Expert_",exp," = as.numeric(X$Contextual.Expert_",exp,")",sep="")
              eval(parse(text=label))
            }
            
          }; rm(exp); rm(label)
          
          cat(" The dataset hasn't missing data ... \n")
          cat(" Continue with the process. \n")
          return(X)
          
        }
        
      }
      x = imp_fun(x)
      
      cat(" Formating data \n")
      mlevels = c("Strongly Disagree","Disagree","Neutral","Agree","Strongly Agree")
      format_data <- function(X, ...){
        
        # ID column
        col_names = names(X)
        id = unlist(lapply(X[col_names], function(x){!is.numeric(x)}))
        id = id[which(id==TRUE)]
        id = names(id)
        rm(col_names)
        
        # Define categories for each variable
        for(j in 1:length(id)){
          
          var = X[,id[j]]
          var = unique(as.character(var))
          
          if(length(var) > 1){
            var = var[na.omit(match(mlevels, var))]
            X[,id[j]] = factor(X[,id[j]], levels=var, ordered=TRUE)
          } else {
            X[,id[j]] = NULL
          }
          
        }; rm(j)
        
        return(X)
        
      }
      x = format_data(x)
      rm(mlevels)
      
      cat(" Applying optimal scaling of categorical variables \n")
      require(aspect)
      for(ex in n_exp){
        n_var = names(x)
        cat_var = c(paste("Evaluation.Expert_",ex,sep=""),
                    paste("Occ_data.Expert_",ex,sep=""),
                    paste("SDM_map.Expert_",ex,sep=""),
                    paste("Gap_map.Expert_",ex,sep=""))
        n_var = unlist(lapply(cat_var,function(x){grep(n_var,pattern=x,fixed=T)}))
        if(length(n_var) > 1){
          optScal = corAspect(x[,n_var], aspect="aspectEigen", level=c(rep("ordinal",length(n_var))))
          x[,n_var] = optScal$scoremat
          follow = TRUE
        } else {
          follow = FALSE
        }
      }; rm(ex); rm(optScal); rm(n_var)
      
      if(follow){
        n_var = names(x)
        szgrp <- unlist(lapply(n_exp, function(x){length(grep(n_var,pattern=x))}))
        
        require(FactoMineR)
        if(length(n_exp) > 1){
          cat(" Applying Factor Multiple Analysis \n")
          results[[tolower(cropList[i])]] = MFA(x, group=szgrp, type=rep("s",length(n_exp)),
                                                ncp=5, name.group=paste("Expert ",n_exp,sep=""),
                                                axes=c(1,2), graph=TRUE)
          cat(" Done. \n")
        } else {
          cat(" Applying Principal Component Analysis \n")
          results[[tolower(cropList[i])]] = PCA(x, scale.unit=T, ncp=5, axes=c(1,2), graph=TRUE)
          cat(" Done. \n")
        }
      } else {
        cat("The data are not suitable for analysis \n")
        cat("Return data matrix. \n")
        results[[tolower(cropList[i])]] = x
      }
      
    } else {
      cat("The number of taxa is too small to run the analysis \n")
      cat("Return data matrix. \n")
      results[[tolower(cropList[i])]] = x
    }
    
    cat("\n")
    
  }
  
  return(results)
  
}
results <- mfa_fun(data_merge,cropList)
results <- results[-1]


# ---------------------------------------------------------------------------- #
# Construction of index
# ---------------------------------------------------------------------------- #

# Cálculo del índice de evaluación
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

# Extracción de datos utilizados para la construcción del índice
extract_fun <- function(data, list){
  
  cropList   <- list
  data_merge <- data
  data_info  <- list(0)
  
  for(i in 1:length(cropList)){
    
    cat("\n =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n")
    cat(" Processing:",tolower(cropList[i]))
    cat("\n =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n")
    cat("\n")
    
    # Organize data horizontal way
    options(warn=-1)
    x = subset(data_merge, data_merge$Crop_code==cropList[i])
    x$Position = NULL
    x$Expert_nm = NULL
    x$Crop_code = NULL
    x$Priority_crop = NULL
    x = reshape(x, timevar="Expert_id", idvar="Taxon", direction="wide")
    
    # Identify Final Priority Score for each taxon
    if(cropList[i] %in% c("beet","dioscorea","groundnut","lettuce","quinoa","safflower","soy","spinach","tomato")){
      taxList = as.character(x$Taxon)
      taxID   = lapply(taxList, function(x){grep(pattern=x,fps_modf$TAXON,fixed=T,value=F)})
      names(taxID) = taxList
      id_null = unlist(lapply(taxID, length))
      if(any(id_null==0)){
        id_null = which(unlist(lapply(taxID, length))==0)
        fps = unlist(lapply(taxID, function(x){min(fps_modf$FPS[x],na.rm=T)}))
        fps[id_null] = 0
      } else {
        fps = unlist(lapply(taxID, function(x){min(fps_modf$FPS[x],na.rm=T)}))
      }
      rm(taxList); rm(taxID); rm(id_null)
    } else {
      taxList = as.character(x$Taxon)
      taxID   = lapply(taxList, function(x){grep(pattern=x,fps_data$TAXON,fixed=T,value=F)})
      names(taxID) = taxList
      id_null = unlist(lapply(taxID, length))
      if(any(id_null==0)){
        id_null = which(unlist(lapply(taxID, length))==0)
        fps = unlist(lapply(taxID, function(x){min(fps_data$FPS[x],na.rm=T)}))
        fps[id_null] = 0
      } else {
        fps = unlist(lapply(taxID, function(x){min(fps_data$FPS[x],na.rm=T)}))
      }
      rm(taxList); rm(taxID); rm(id_null)
    }
    
    # Add taxon name to rows of data frame
    rownames(x) = x$Taxon
    x$Taxon = NULL
    
    # Identify number of experts that responded the surveys
    n_exp = names(x)
    n_exp = unlist(lapply(n_exp, function(x){strsplit(x, split="_")}))
    n_chr = unlist(lapply(n_exp, function(x){nchar(x)}))
    n_exp = as.numeric(n_exp[which(n_chr==1)]); rm(n_chr)
    n_exp = unique(n_exp)
    
    # Add FPS to matrix data
    x$FPS <- fps
    
    cat(" Calculate relative difference \n")
    for(ex in n_exp){
      label = paste("x$Comparable.Trans.Expert_",ex," = 1-abs((fps-x$Comparable.Expert_",ex,")/10)",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Comparable.Trans.Expert_",ex," = as.numeric(x$Comparable.Trans.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Contextual.Trans.Expert_",ex," = 1-abs((fps-x$Contextual.Expert_",ex,")/10)",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Contextual.Trans.Expert_",ex," = as.numeric(x$Contextual.Trans.Expert_",ex,")",sep="")
      eval(parse(text=label))
      label = paste("x$Evaluation.Expert_",ex," = as.character(x$Evaluation.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Occ_data.Expert_",ex," = as.character(x$Occ_data.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$SDM_map.Expert_",ex," = as.character(x$SDM_map.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Gap_map.Expert_",ex," = as.character(x$Gap_map.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
    }; rm(ex)
    
    # Exclude for analysis datasets with less than 6 taxa
    if(nrow(x) > 6){
      
      cat(" Filtering data \n")
      
      filter_fun <- function(X){
        
        na_var = apply(X, 2, function(x){sum(is.na(x))})
        mtch   = names(na_var[which(na_var!=nrow(X))]); rm(na_var)
        X      = X[,mtch]; rm(mtch)
        
        lapply(names(X), function(x){
          label  = paste("sd(as.numeric(as.factor(X[,","'",x,"'","])),na.rm=T)",sep="")
          sd_val = eval(parse(text=label)); rm(label)
          if(sd_val == 0){
            cat(" ",x,"variable hasn't associated variability. Thus eliminiating. \n")
            label = paste("X$",x," = NULL",sep="")
            eval(parse(text=label)); rm(label)
          } else{
            cat(" ",x,"variable is ok for the analysis \n")
          }
        })
        
        return(X)
        
      }
      x = filter_fun(x)
      
      # No es necesario
      #cat(" Verify if exists missing data, then impute \n")
      
      cat(" Formating data \n")
      mlevels = c("Strongly Disagree","Disagree","Neutral","Agree","Strongly Agree")
      format_data <- function(X, ...){
        
        # ID column
        col_names = names(X)
        id = unlist(lapply(X[col_names], function(x){!is.numeric(x)}))
        id = id[which(id==TRUE)]
        id = names(id)
        rm(col_names)
        
        # Define categories for each variable
        for(j in 1:length(id)){
          
          var = X[,id[j]]
          var = unique(as.character(var))
          
          if(length(var) > 1){
            var = var[na.omit(match(mlevels, var))]
            X[,id[j]] = factor(X[,id[j]], levels=var, ordered=TRUE)
          } else {
            X[,id[j]] = NULL
          }
          
        }; rm(j)
        
        return(X)
        
      }
      x = format_data(x)
      rm(mlevels)
      
      cat(" Applying optimal scaling of categorical variables \n")
      for(ex in n_exp){
        n_var = names(x)
        cat_var = c(paste("Evaluation.Expert_",ex,sep=""),
                    paste("Occ_data.Expert_",ex,sep=""),
                    paste("SDM_map.Expert_",ex,sep=""),
                    paste("Gap_map.Expert_",ex,sep=""))
        n_var = unlist(lapply(cat_var,function(x){grep(n_var,pattern=x,fixed=T)}))
        if(length(n_var) > 1){
          follow = TRUE
        } else {
          follow = FALSE
        }
      }; rm(ex)
      
      if(follow){
        n_var = names(x)
        szgrp <- unlist(lapply(n_exp, function(x){length(grep(n_var,pattern=x))}))
        
        if(length(n_exp) > 1){
          cat(" More of one expert. \n")
          data_info[[tolower(cropList[i])]] = x
          cat(" Done. \n")
        } else {
          cat(" One expert. \n")
          data_info[[tolower(cropList[i])]] = x
          cat(" Done. \n")
        }
      } else {
        cat("The data are not suitable for analysis \n")
        cat("Return data matrix. \n")
        data_info[[tolower(cropList[i])]] = x
      }
      
    } else {
      cat("The number of taxa is too small to run the analysis \n")
      cat("Return data matrix. \n")
      data_info[[tolower(cropList[i])]] = x
    }
    
    cat("\n")
    
  }
  
  return(data_info)
  
}
data_info <- extract_fun(data_merge,cropList)
data_info <- data_info[-1]

# Unir la información
final_data <- mapply(index_res, data_info, FUN=function(x,y){y$Index <- x; return(y)})
mapply(final_data, names(final_data), FUN=function(x,y){write.csv(x,file=paste("F:/CIAT/CWR_project/Experts_evaluation/Final_results/",y,".csv",sep=""),row.names=T)})

# Read data organized
crops <- list.files("F:/CIAT/CWR_project/Experts_evaluation/Final_results")
all_data <- lapply(crops, function(x){read.csv(paste("F:/CIAT/CWR_project/Experts_evaluation/Final_results/",x,sep=""),header=T,row.names=1)})
all_data <- lapply(all_data, function(x){x <- x[,c("FPS","Index")]})

all_data <- Reduce(function(...) merge(..., all=T), all_data)
plot(all_data, xlab="Final Priority Score", ylab="Evaluation Index")
abline(lm(Index~FPS,data=all_data))

cor(all_data, use="complete.obs")

# ---------------------------------------------------------------------------- #
# Merge all information by crop
# ---------------------------------------------------------------------------- #










