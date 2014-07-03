
# Download experts data
# H. Achicanoy
# CIAT, 2014

# ---------------------------------------------------------------------------- #
### INPUTS
# data_kind: "scores" or "maps"
# final_dir: path for save the results
# email: e-mail direction for google (Full!)
# password: password of your e-mail direction
# ---------------------------------------------------------------------------- #
### EXAMPLES
# download_data(data_kind="scores",
#               final_dir="C:/Users/haachicanoy/Documents/Experts_evaluation",
#               email="harold22010@gmail.com",
#               password="**********")
# download_data(data_kind="maps",
#               final_dir="C:/Users/haachicanoy/Documents/Experts_evaluation",
#               email="harold22010@gmail.com",
#               password="**********")
# ---------------------------------------------------------------------------- #

download_data = function(data_kind, final_dir, email, password, ...){
  
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  # Load packages & source code
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  
  require(XML)
  require(RCurl)
  require(RGoogleDocs)
  
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  # Read data online
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  
  connection   = getGoogleDocsConnection(getGoogleAuth(paste0(email),paste0(password),service="wise"))
  spreadsheets = getDocs(connection, ssl.verifypeer=FALSE)
  sheets       = names(spreadsheets)
  
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  # Classify data
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  
  # Priority crops
  pList <- c("Avena","Cajanus","Cicer","Daucus","Eggplant",
             "Eleusine","Helianthus","Hordeum","Ipomoea",
             "Lathyrus","Lens","Malus","Medicago","Musa",
             "Pennisetum","Phaseolus","Pisum","Potato",
             "Rice","Secale","Sorghum","Vicia","Vigna",
             "Wheat")
  
  # List of priority taxa (Utilizar CWR_GA_Data)
  taxaPrior = getWorksheets(spreadsheets[["PtaxaofPcrops_all"]],connection)
  comp      = taxaPrior$"names_codes_new455"
  taxaPrior = sheetAsMatrix(comp,header=T,as.data.frame=T,trim=T)
  rm(comp)
  
  # Non priority crops
  npList <- c("allium","asparagus","beet","brassica","breadfruit",
              "cacao","capsicum","cassava","citrus","cocoyam",
              "cotton","cucumis","dioscorea","grape","groundnut",
              "lettuce","maize","mango","millet_panicum",
              "millet_setaria","papaya","pear","pineapple",
              "prunus","quinoa","safflower","soy","spinach",
              "squash","strawberry","sugar_cane","tomato","vigna",
              "watermelon")
  
  # List of non priority taxa (Utilizar NonP_PrioritiesTables)
  taxaNonPrior = getWorksheets(spreadsheets[["CWR_non-priority_crops_processing"]],connection)
  comp         = taxaNonPrior$"non_priority_taxa-2013-09-13"
  taxaNonPrior = sheetAsMatrix(comp,header=T,as.data.frame=T,trim=T)
  taxaNonPrior = taxaNonPrior[,1:12]
  rm(comp)
  
  # Publications crops
  publs <- c("Helianthus(Paper)")
  
  ### DOWNLOAD SCORES DATA ###
  if(data_kind == "scores"){
    
    g_shts = sheets[grep(pattern="CWR_Expert_Evaluation_", sheets)]
    g_shts = unique(g_shts); g_shts = sort(g_shts)
    
    # =-=-= Priority crops
    gp_shts = g_shts[unlist(sapply(pList,function(x){grep(pattern=paste("_",x," ",sep=""),g_shts,fixed=T)}))] #"\\", "$"
    
    # =-=-= Non priority
    gnp_shts = g_shts[unlist(sapply(npList,function(x){grep(pattern=paste("_",x," ",sep=""),g_shts,fixed=T)}))] #"\\"
    
    # =-=-= Publications crops
    pub_shts = g_shts[unlist(sapply(publs,function(x){grep(pattern=paste("_",x," ",sep=""),g_shts,fixed=T)}))] #"\\", "$"
    
    # Organize information
    g_crop = unlist(strsplit(g_shts, split=" "))
    g_crop = g_crop[seq(1,length(g_crop),by=2)]
    g_crop = unlist(strsplit(g_crop, split="CWR_Expert_Evaluation_",fixed=T))
    g_crop = g_crop[seq(2,length(g_crop),by=2)]
    
    require(gdata)
    
    count = 1
    for(i in 1:length(g_shts)){
      
      cat("\n")
      
      cat("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n")
      cat("Procesing information for",g_shts[i],"\n")
      cat("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n")
      
      cat("\n")
      
      # Crear conexión con el archivo en gdocs
      con1 = getWorksheets(spreadsheets[[g_shts[i]]],connection)
      comp = con1$"Form Responses"
      con2 = try(sheetAsMatrix(comp,header=T,as.data.frame=T,trim=T),silent=T)
      
      # Revisar si hay datos
      if(is.data.frame(con2)){
        
        # Revisar si los datos digitados por cada experto son validos
        aux = 1
        for(nr in 1:nrow(con2)){
          
          info = as.vector(t(con2[nr,])) # Extraer la información de cada experto
          responses = c("High Priority: 0","No Need for Further Collection: 10","Strongly Disagree","Disagree","Neutral","Agree","Strongly Agree",paste(1:9,sep=""))
          
          if(!any(info %in% responses)){
            
            cat("The information typed in the row",nr,"isn't valid \n")
            cat("Storing information \n")
            
            if(aux == 1){
              rows = nr
            } else {
              rows = c(rows,nr)
              rows = as.numeric(rows)
            }
            
            aux = aux + 1
            
          } else {
            cat("The information typed in the row",nr,"is valid \n")
          }
          
        }
        
        cat("\n")
        cat("\n")
        rm(aux); rm(nr)
        rm(info); rm(responses)
        
        # Si existe información invalida, eliminar los registros
        if(exists("rows")){
          cat("Delete information \n")
          cat("\n")
          cat("\n")
          con2 = con2[-rows,]
          rm(rows)
        }
        
        rownames(con2) = 1:nrow(con2)
        
        # Etiquetar la información de cada experto
        for(j in 1:nrow(con2)){
          con2$Expert[j] = paste("Expert_",j,sep="")
        }; rm(j)
        
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        ### TEXTUAL DATA
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        
        data_crop = data.frame(Expert_id       = con2$Expert,
                               Expert_nm       = con2$Name,
                               Position        = con2$"Position title",
                               Institution     = con2$Institution,
                               e_mail          = con2$Email,
                               Notes           = con2$Notes,
                               Prof_background = if(is.null(con2$"Professional Background")){NA} else{con2$"Professional Background"},
                               Additional_taxa = if(is.null(con2$Additional_Taxa)){NA} else{con2$Additional_Taxa},
                               Eval_gap_scores = if(is.null(con2$"Evaluation of Gap Analysis Scores")){NA} else{con2$"Evaluation of Gap Analysis Scores"},
                               Comments_gap_analysis = if(is.null(con2$"Comments on Gap Analysis Scores")){NA} else{con2$"Comments on Gap Analysis Scores"})
        
        data_crop$Priority_crop = if(g_crop[i] %in% pList){1}else{0}
        data_crop$Crop_code = g_crop[i]
        
        if(count == 1){
          text_data = data_crop
        } else {
          text_data = rbind(text_data,data_crop)
        }
        
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        ### QUANTITATIVE DATA
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        
        # List of taxa by crop
        if(g_crop[i] %in% pList){ # PRIORITY CROPS
          if(g_crop[i] == "Phaseolus"){ # Phaseolus case
            taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Bean")]),
                             as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Lima_bean")])))
          } else {
            if(g_crop[i] == "Vicia"){   # Vicia case
              taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Faba_bean")]),
                               as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Vetch")])))
            } else {
              if(g_crop[i] == "Wheat"){ # Wheat case
                taxList = sort(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Triticum")]))
              } else {
                if(g_crop[i] == "Vigna"){ # Vigna case
                  taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Bambara")]),
                                   as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Cowpea")])))
                } else {
                  taxList = sort(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower(g_crop[i])]))
                }
              }
            }
          }
        } else {
          if(g_crop[i] %in% npList){ # NON PRIORITY CROPS
            if(g_crop[i] == "dioscorea"){ # Dioscorea case
              taxList = sort(c(as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("yam_lagos")]),
                               as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("yam_water")]),
                               as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("yam_whiteguinea")]),
                               "Dioscorea_cayennensis_subsp._rotundata"))
              taxList = unique(taxList)
            } else {
              if(g_crop[i] == "soy"){
                taxList = sort(as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("soybean")]))
              } else {
                taxList = sort(as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower(g_crop[i])]))
              }
            }
          } else {
            if(g_crop[i] %in% publs){ # PUBLICATION CROPS
              taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Helianthus")]),"Helianthus_pauciflorus_subsp._pauciflorus","Helianthus_pauciflorus_subsp._subrhomboideus","Helianthus_winteri"))
            }
          }
        }
        
        taxList = as.vector(taxList)
        cols    = as.vector(colnames(con2)) # Names of columns in data set
        
        # Intersect names of taxa
        mtch = unlist(sapply(taxList, function(x){grep(pattern=paste("[",x,"]",sep=""), cols, fixed=TRUE)}))
        
        if(is.null(dim(mtch))){
          mtch = unlist(sapply(names(mtch),FUN=function(x){substr(x, 1, nchar(x)-1)}))
          mtch = sort(unique(mtch))
        } else {
          mtch = colnames(mtch)
        }
        
        #--- Comparable data field
        comp = matchcols(con2, with=paste("Comparable Expert Priority Score [",mtch,"]",sep=""), method="or", fixed=T)
        names(comp) = NULL; comparable = con2[,comp]; rm(comp)
        
        comparable = t(comparable); colnames(comparable) = con2$Expert
        comparable[which(comparable=="High Priority: 0")] <- "0"
        comparable[which(comparable=="No Need for Further Collection: 10")] <- "10"
        options(warn=-1)
        class(comparable) = "numeric"
        
        comparable = as.data.frame(comparable)
        comparable = stack(comparable, select=colnames(comparable))
        colnames(comparable) = c("Comparable","Expert_id")
        
        #--- Contextual field
        cont = matchcols(con2, with=paste("Contextual Expert Priority Score [",mtch,"]",sep=""), method="or", fixed=T)
        names(cont) = NULL; contextual = con2[,cont]; rm(cont)
        
        contextual = t(contextual); colnames(contextual) = con2$Expert
        contextual[which(contextual=="High Priority: 0")] <- "0"
        contextual[which(contextual=="No Need for Further Collection: 10")] <- "10"
        options(warn=-1)
        class(contextual) = "numeric"
        
        contextual = as.data.frame(contextual)
        contextual = stack(contextual, select=colnames(contextual))
        colnames(contextual) = c("Contextual","Expert_id")
        
        #--- Evaluation field
        eval = matchcols(con2, with=paste("Evaluation of Gap Analysis Scores [",mtch,"]",sep=""), method="or", fixed=T)
        names(eval) = NULL; evaluation = con2[,eval]; rm(eval)
        
        evaluation$Expert = con2$Expert
        evaluation$time = 1
        
        evaluation = reshape(evaluation, direction="long", v.names="Evaluation", idvar="Expert", varying=1:length(mtch))
        evaluation = evaluation[order(evaluation$Expert),]
        rownames(evaluation) = 1:nrow(evaluation)
        evaluation$time = NULL
        
        #--- Position field
        ne = 1
        for(j in 1:nrow(con2)){
          if(ne == 1){
            position = rep(con2$"Position title"[j],length(mtch))
            position = as.character(position)
          } else {
            position = c(position, as.character(rep(con2$"Position title"[j],length(mtch))))
          }
          ne = ne + 1
        }; rm(j); rm(ne)
        
        #--- Name field
        ne = 1
        for(j in 1:nrow(con2)){
          if(ne == 1){
            name = rep(con2$Name[j],length(mtch))
            name = as.character(name)
          } else {
            name = c(name, as.character(rep(con2$Name[j],length(mtch))))
          }
          ne = ne + 1
        }; rm(j); rm(ne)
        
        data_taxon = data.frame(Expert_id     = comparable$Expert_id,
                                Expert_nm     = name,
                                Position      = position,
                                Comparable    = comparable$Comparable,
                                Contextual    = contextual$Contextual,
                                Evaluation    = evaluation$Evaluation,
                                Priority_crop = if(g_crop[i] %in% pList){1}else{0},
                                Crop_code     = g_crop[i],
                                Taxon         = rep(mtch,nrow(con2)))
        
        rm(list=c("name","position","comparable","contextual","evaluation"))
        
        if(count == 1){
          quan_data = data_taxon
        } else {
          quan_data = rbind(quan_data,data_taxon)
        }
        
        rm(cols); rm(taxList)
        
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        
        count = count + 1
        
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        
      } else {
        cat("Awaiting for experts evaluation \n")
        cat("\n")
        cat("\n")
        rm(con2)
      }
      
    }
    
    rm(con1); rm(con2); rm(count)
    rm(data_crop); rm(data_taxon)
    rm(i); rm(mtch)
    
    sumData <- data.frame(Crop_code = tolower(unique(quan_data$Crop_code)),
                          nTaxa     = unlist(sapply(as.character(unique(quan_data$Crop_code)),FUN=function(x){length(unique(quan_data$Taxon[quan_data$Crop_code==x]))})),
                          nExperts  = unlist(sapply(as.character(unique(quan_data$Crop_code)),FUN=function(x){length(unique(quan_data$Expert_id[quan_data$Crop_code==x]))})),
                          Priority  = unlist(sapply(as.character(unique(quan_data$Crop_code)),FUN=function(x){unique(quan_data$Priority_crop[quan_data$Crop_code==x])})))
    rownames(sumData) <- 1:nrow(sumData)
    
    # Print results
    
    if(!exists(final_dir)){
      dir.create(final_dir)
    }
    
    write.csv(quan_data, paste0(final_dir,"/scoresEval.csv"), row.names=F)
    write.csv(text_data, paste0(final_dir,"/scores_textEval.csv"), row.names=F)
    write.csv(sumData,   paste0(final_dir,"/summ_scoresEval.csv"), row.names=F)
    
  }
  
  ### DOWNLOAD MAPS DATA ###
  if(data_kind == "maps"){
    
    m_shts   = sheets[grep(pattern="CWR_Map_Evaluation_", sheets)]
    m_shts   = sort(m_shts)
    
    # =-=-= Priority crops
    mp_shts  = m_shts[unlist(sapply(pList,function(x){grep(pattern=paste("_",x," ",sep=""),m_shts,fixed=T)}))]
    
    # =-=-= Non priority
    mnp_shts = m_shts[unlist(sapply(npList,function(x){grep(pattern=paste("_",x," ",sep=""),m_shts,fixed=T)}))]
    
    # =-=-= Publications crops
    mpb_shts = m_shts[unlist(sapply(publs,function(x){grep(pattern=paste("_",x," ",sep=""),m_shts,fixed=T)}))]
    
    g_cropM = unlist(strsplit(m_shts, split=" "))
    g_cropM = g_cropM[seq(1,length(g_cropM),by=2)]
    g_cropM = unlist(strsplit(g_cropM, split="CWR_Map_Evaluation_",fixed=T))
    g_cropM = g_cropM[seq(2,length(g_cropM),by=2)]
    
    # Solo descargar datos para taxones
    
    require(gdata)
    
    count = 1
    for(i in 1:length(m_shts)){
      
      cat("\n")
      
      cat("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n")
      cat("Procesing information for",m_shts[i],"\n")
      cat("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n")
      
      cat("\n")
      
      # Crear conexión con el archivo en gdocs
      con1 = getWorksheets(spreadsheets[[m_shts[i]]],connection)
      comp = con1$"Form Responses"
      con2 = try(sheetAsMatrix(comp,header=T,as.data.frame=T,trim=T),silent=T)
      
      # Revisar si hay datos
      if(is.data.frame(con2)){
        
        # Revisar si los datos digitados por cada experto son validos
        aux = 1
        for(nr in 1:nrow(con2)){
          
          info = as.vector(t(con2[nr,])) # Extraer la información de cada experto
          responses = c("Strongly Disagree","Disagree","Neutral","Agree","Strongly Agree") # paste(1:5,sep="")
          
          if(!any(info %in% responses)){
            
            cat("The information typed in the row",nr,"isn't valid \n")
            cat("Storing information \n")
            
            if(aux == 1){
              rows = nr
            } else {
              rows = c(rows,nr)
              rows = as.numeric(rows)
            }
            
            aux = aux + 1
            
          } else {
            cat("The information typed in the row",nr,"is valid \n")
          }
          
        }
        
        cat("\n")
        cat("\n")
        rm(aux); rm(nr)
        rm(info); rm(responses)
        
        # Si existe información invalida, eliminar los registros
        if(exists("rows")){
          cat("Delete information \n")
          cat("\n")
          cat("\n")
          con2 = con2[-rows,]
          rm(rows)
        }
        
        rownames(con2) = 1:nrow(con2)
        
        # Etiquetar la información de cada experto
        for(j in 1:nrow(con2)){
          con2$Expert[j] = paste("Expert_",j,sep="")
        }; rm(j)
        
        c_names = colnames(con2)
        if(!(length(unique(c_names))==length(c_names))){
          con2 = con2[,!is.na(con2[1,])]
        }; rm(c_names)
        
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        ### TEXTUAL DATA
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        
        #data_crop = data.frame(Expert_id         = con2$Expert,
        #                       Expert_nm         = con2$Name,
        #                       Position          = con2$"Position Title",
        #                       Institution       = con2$Institution,
        #                       e_mail            = con2$Email,
        #                       Regions_expertise = if(is.null(con2$"Regions of Expertise")){NA} else{con2$"Regions of Expertise"},
        #                       Prof_background   = if(is.null(con2$"Professional Background")){NA} else{con2$"Professional Background"},
        #                       Comm_Occ_data     = if(is.null(con2$"Comments on Occurrence Data")){NA} else{con2$"Comments on Occurrence Data"},
        #                       Comm_SDM_maps     = if(is.null(con2$"Comments on Potential Distribution Models")){NA} else{con2$"Comments on Potential Distribution Models"},
        #                       Comm_Gap_maps     = if(is.null(con2$"Comments on Collecting Priorities Maps")){NA} else{con2$"Comments on Collecting Priorities Maps"},
        #                       Eval_Taxon_rich   = if(is.null(con2$"Evaluation of Crop Gene Pool Level Maps [Taxon Richness]")){NA} else{con2$"Evaluation of Crop Gene Pool Level Maps [Taxon Richness]"},
        #                       Eval_Coll_hot     = if(is.null(con2$"Evaluation of Crop Gene Pool Level Maps [Collecting Hotspots]")){NA} else{con2$"Evaluation of Crop Gene Pool Level Maps [Collecting Hotspots]"},
        #                       Comm_Taxon_rich   = if(is.null(con2$"Comments on Gene Pool Taxon Richness Map")){NA} else{con2$"Comments on Gene Pool Taxon Richness Map"},
        #                       Comm_Coll_hot     = if(is.null(con2$"Comments on Gene Pool Collecting Hotspots Map")){NA} else{con2$"Comments on Gene Pool Collecting Hotspots Map"},
        #                       Eval_interface    = if(is.null(con2$"Map Interface")){NA} else{con2$"Map Interface"},
        #                       Comm_interface    = if(is.null(con2$"Comments on Map Interface")){NA} else{con2$"Comments on Map Interface"})
        
        #data_crop$Priority_crop = if(g_cropM[i] %in% pList){1}else{0}
        #data_crop$Crop_code = g_cropM[i]
        
        #if(count == 1){
        #  text_data = data_crop
        #} else {
        #  text_data = rbind(text_data,data_crop)
        #}
        
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        ### QUALITATIVE DATA
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        
        ### Data by taxon
        # Priority_crop, Crop_code, Taxon, Expert_id
        # Expert_nm, Position, Occ_data, SDM_map, Gap_map
        
        # List of taxa by crop
        if(g_cropM[i] %in% pList){ # PRIORITY CROPS
          if(g_cropM[i] == "Phaseolus"){ # Phaseolus case
            taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Bean")]),
                             as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Lima_bean")])))
          } else {
            if(g_cropM[i] == "Vicia"){   # Vicia case
              taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Faba_bean")]),
                               as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Vetch")])))
            } else {
              if(g_cropM[i] == "Wheat"){ # Wheat case
                taxList = sort(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Triticum")]))
              } else {
                if(g_cropM[i] == "Vigna"){ # Vigna case
                  taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Bambara")]),
                                   as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Cowpea")])))
                } else {
                  if(g_cropM[i] == "Avena"){
                    taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Avena")]),
                                     "Avena_sterilis_subsp._ludoviciana"))
                  } else {
                    taxList = sort(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower(g_cropM[i])]))
                  }
                }
              }
            }
          }
        } else {
          if(g_cropM[i] %in% npList){ # NON PRIORITY CROPS
            if(g_cropM[i] == "dioscorea"){ # Dioscorea case
              taxList = sort(c(as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("yam_lagos")]),
                               as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("yam_water")]),
                               as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("yam_whiteguinea")]),
                               "Dioscorea_cayennensis_subsp._rotundata"))
              taxList = unique(taxList)
            } else {
              if(g_cropM[i] == "soy"){
                taxList = sort(as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("soybean")]))
              } else {
                taxList = sort(as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower(g_cropM[i])]))
              }
            }
          } else {
            if(g_cropM[i] %in% publs){ # PUBLICATION CROPS
              taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Helianthus")]),"Helianthus_pauciflorus_subsp._pauciflorus","Helianthus_pauciflorus_subsp._subrhomboideus","Helianthus_winteri"))
            }
          }
        }
        
        taxList = as.vector(taxList)
        cols    = as.vector(colnames(con2)) # Names of columns in data set
        
        # Intersect names of taxa
        mtch = unlist(sapply(taxList, function(x){grep(pattern=paste("[",x,"]",sep=""), cols, fixed=TRUE)}))
        
        if(is.null(dim(mtch))){
          
          mtch = names(mtch)
          options(warn=-1)
          
          for(j in 1:length(mtch)){
            if(!is.na(as.numeric(substr(mtch[j],nchar(mtch[j]),nchar(mtch[j]))))){
              mtch[j] = substr(mtch[j], 1, nchar(mtch[j])-1)
            }
          }; rm(j)
          
          mtch = sort(unique(mtch))
          
        } else {
          mtch = colnames(mtch)
        }
        
        #--- Occurrence data
        occ = matchcols(con2, with=paste("Evaluation of Occurrence Data [",mtch,"]",sep=""), method="or", fixed=T)
        names(occ) = NULL
        
        # check = unlist(lapply(occ,function(x){identical(x,character(0))}))
        # checkList = data.frame(Taxon=mtch,Check=check); rm(check)
        
        occ = unlist(occ)
        occ_data = con2[,occ]
        
        if(is.factor(occ_data)){
          occ_data = as.data.frame(occ_data)
          colnames(occ_data) = occ
        }
        
        for(j in 1:length(mtch)){
          if(!paste("Evaluation of Occurrence Data [",mtch[j],"]",sep="") %in% colnames(occ_data)){
            eq = paste("occ_data$'Evaluation of Occurrence Data [",mtch[j],"]' <- NA",sep="")
            eval(parse(text=eq))
          }
        }; rm(j); rm(eq)
        
        col_names = sort(names(occ_data))
        occ_data = occ_data[,col_names]
        
        if(is.factor(occ_data)){
          occ_data = as.data.frame(occ_data)
          colnames(occ_data) = col_names
        }
        
        rm(occ); rm(col_names)
        
        occ_data$Expert = con2$Expert
        occ_data$time = 1
        
        occ_data = reshape(occ_data, direction="long", v.names="Occ_data", idvar="Expert", varying=1:length(mtch))
        occ_data = occ_data[order(occ_data$Expert),]
        rownames(occ_data) = 1:nrow(occ_data)
        occ_data$time = NULL
        
        #--- Potential distribution model
        sdm = matchcols(con2, with=paste("Evaluation of Potential Distribution Models [",mtch,"]",sep=""), method="or", fixed=T)
        names(sdm) = NULL
        
        sdm = unlist(sdm)
        sdm_map = con2[,sdm]
        
        if(is.factor(sdm_map)){
          sdm_map = as.data.frame(sdm_map)
          colnames(sdm_map) = sdm
        }
        
        for(j in 1:length(mtch)){
          if(!paste("Evaluation of Potential Distribution Models [",mtch[j],"]",sep="") %in% colnames(sdm_map)){
            eq = paste("sdm_map$'Evaluation of Potential Distribution Models [",mtch[j],"]' <- NA",sep="")
            eval(parse(text=eq))
          }
        }; rm(j); rm(eq)
        
        col_names = sort(names(sdm_map))
        sdm_map = sdm_map[,col_names]
        
        if(is.factor(sdm_map)){
          sdm_map = as.data.frame(sdm_map)
          colnames(sdm_map) = col_names
        }
        
        rm(sdm); rm(col_names)
        
        sdm_map$Expert = con2$Expert
        sdm_map$time = 1
        
        sdm_map = reshape(sdm_map, direction="long", v.names="SDM_map", idvar="Expert", varying=1:length(mtch))
        sdm_map = sdm_map[order(sdm_map$Expert),]
        rownames(sdm_map) = 1:nrow(sdm_map)
        sdm_map$time = NULL
        
        #--- Collecting Prioroties maps
        gap = matchcols(con2, with=paste("Evaluation of Collecting Priorities Maps [",mtch,"]",sep=""), method="or", fixed=T)
        names(gap) = NULL
        
        gap = unlist(gap)
        gap_map = con2[,gap]
        
        if(is.factor(gap_map)){
          gap_map = as.data.frame(gap_map)
          colnames(gap_map) = gap
        }
        
        for(j in 1:length(mtch)){
          if(!paste("Evaluation of Collecting Priorities Maps [",mtch[j],"]",sep="") %in% colnames(gap_map)){
            eq = paste("gap_map$'Evaluation of Collecting Priorities Maps [",mtch[j],"]' <- NA",sep="")
            eval(parse(text=eq))
          }
        }; rm(j); rm(eq)
        
        col_names = sort(names(gap_map))
        gap_map = gap_map[,col_names]
        
        if(is.factor(gap_map)){
          gap_map = as.data.frame(gap_map)
          colnames(gap_map) = col_names
        }
        
        rm(gap); rm(col_names)
        
        gap_map $Expert = con2$Expert
        gap_map $time = 1
        
        gap_map  = reshape(gap_map , direction="long", v.names="Gap_map", idvar="Expert", varying=1:length(mtch))
        gap_map  = gap_map [order(gap_map $Expert),]
        rownames(gap_map ) = 1:nrow(gap_map )
        gap_map $time = NULL
        
        #--- Position
        ne = 1
        for(j in 1:nrow(con2)){
          if(ne == 1){
            position = rep(con2$"Position Title"[j],length(mtch))
            position = as.character(position)
          } else {
            position = c(position, as.character(rep(con2$"Position Title"[j],length(mtch))))
          }
          ne = ne + 1
        }; rm(j); rm(ne)
        
        #--- Name
        ne = 1
        for(j in 1:nrow(con2)){
          if(ne == 1){
            name = rep(con2$Name[j],length(mtch))
            name = as.character(name)
          } else {
            name = c(name, as.character(rep(con2$Name[j],length(mtch))))
          }
          ne = ne + 1
        }; rm(j); rm(ne)
        
        # Priority_crop, Crop_code, Taxon, Expert_id
        # Expert_nm, Position, Occ_data, SDM_map, Gap_map
        
        data_taxon = data.frame(Expert_id     = occ_data$Expert,
                                Expert_nm     = name,
                                Position      = position,
                                Occ_data      = occ_data$Occ_data,
                                SDM_map       = sdm_map$SDM_map,
                                Gap_map       = gap_map$Gap_map,
                                Priority_crop = if(g_cropM[i] %in% pList){1}else{0},
                                Crop_code     = g_cropM[i],
                                Taxon         = rep(mtch,nrow(con2)))
        
        rm(list=c("name","position","occ_data","sdm_map","gap_map"))
        
        if(count == 1){
          quan_data = data_taxon
        } else {
          quan_data = rbind(quan_data,data_taxon)
        }
        
        rm(cols); rm(taxList)
        
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        
        count = count + 1
        
        ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
        
      } else {
        cat("Awaiting for experts evaluation \n")
        cat("\n")
        cat("\n")
        rm(con2)
      }
      
    }
    
    rm(con1); rm(con2); rm(count)
    rm(data_crop); rm(data_taxon)
    rm(i); rm(mtch)
    
    sumData <- data.frame(Crop_code = tolower(unique(quan_data$Crop_code)),
                          nTaxa     = unlist(sapply(as.character(unique(quan_data$Crop_code)),FUN=function(x){length(unique(quan_data$Taxon[quan_data$Crop_code==x]))})),
                          nExperts  = unlist(sapply(as.character(unique(quan_data$Crop_code)),FUN=function(x){length(unique(quan_data$Expert_id[quan_data$Crop_code==x]))})),
                          Priority  = unlist(sapply(as.character(unique(quan_data$Crop_code)),FUN=function(x){unique(quan_data$Priority_crop[quan_data$Crop_code==x])})))
    rownames(sumData) <- 1:nrow(sumData)
    
    # Print results
    
    write.csv(quan_data, paste0(final_dir,"/mapsEval.csv"), row.names=F)
    write.csv(sumData,   paste0(final_dir,"/summ_mapsEval.csv"), row.names=F)
    
  }
  
  return(cat("Done!"))
  
}

