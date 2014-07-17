# mfa_fun original
# H. Achicanoy
# CIAT, 2014

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
      label = paste("x$Comparable.Expert_",ex," = (x$Comparable.Expert_",ex,"-fps)/10",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Comparable.Expert_",ex," = as.numeric(x$Comparable.Expert_",ex,")",sep="")
      eval(parse(text=label)); rm(label)
      label = paste("x$Contextual.Expert_",ex," = (x$Contextual.Expert_",ex,"-fps)/10",sep="")
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
