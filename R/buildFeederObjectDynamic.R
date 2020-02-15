#
#  This file is part of the CNO software
#
#  Copyright (c) 2018 - RWTH Aachen - JRC COMBINE
#
#  File author(s): E.Gjerga (enio.gjerga@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$

# This function estimates the possible mechanisms of interactions to be added to the PKN
# from a database of interactions for improving the fitting cost.

# Inputs:
# Mandatory:  A cnolist object containing the data (cnolist)
#             A model to optimize (model)
#             An indices object (i.e. as returned from computeMSE.R function) indicating poor fits
# Optional:   A path length parameter for the maximal path length of additional interactions (pathLength = 3 by default)
#             A database of interactions
#             A 2 value vector (err1, err2) defining the error model of the data as sd^2 = err1^2 + (err2*data)^2, default to c(0.1, 0)
#             A parameter that determine the threshold of significancy of the effect of stimuli and inhibitors, default to 2

buildFeederObjectDynamic <- function(model = model, cnolist = cnolist, 
                                     indices = indices, database = NULL, 
                                     DDN = TRUE, pathLength = 3, k = 2, 
                                     measErr = c(0.1, 0), timePoint = NA){
  
  ##
  # Intial checks
  if(is.null(database) && DDN==FALSE){
    stop("You should either allow data-driven inference or make use of a 
         database of interactions or both in order to integrate the new links 
         in the PKN!")
  }
  
  ##
  # Identifying all the erroneous measurements
  indMeas = unique(names(indices$indices))
  
  ##
  # Identifying the interactions from the FEED algorithm
  # BTable <- makeBTables(CNOlist=cnolist, k=k, measErr=measErr)
  BTable <- makeBTables(CNOlist = cnolist, k = k, measErr = measErr, 
                        timePoint = timePoint)
  # for(ii in 1:length(BTable$tables)){
  #   currMeas = names(BTable$tables)[ii]
  #   if(!(currMeas%in%indMeas)){
  #     for(jj in 1:nrow(BTable$tables[[ii]])){
  #       for(kk in 1:ncol(BTable$tables[[ii]])){
  #         
  #         BTable$tables[[ii]][jj, kk] = 0
  #         BTable$NotMatStim[[ii]][jj, kk] = 0
  #         BTable$NotMatInhib[[ii]][jj, kk] = 0
  #         
  #       }
  #     }
  #   }
  # }
  
  #doing integration with FEEDER
  modelIntegr <- mapBTables2model(BTable=BTable,model=model,allInter=TRUE)
  integrSIF <- model2sif(modelIntegr)
  
  #now retaining only those interactions involving poorly fitted measurements
  reactions = modelIntegr$reacID[modelIntegr$indexIntegr]
  reacTable = matrix(data = , nrow = length(reactions), ncol = 2)
  for(ii in 1:length(reactions)){
    reacTable[ii, 1] = gsub(pattern = "!", replacement = "", 
                            x = strsplit(x = reactions[ii], split = "=", 
                                         fixed = TRUE)[[1]][1])
    reacTable[ii, 2] = gsub(pattern = "!", replacement = "", 
                            x = strsplit(x = reactions[ii], split = "=", 
                                         fixed = TRUE)[[1]][2])
  }
  
  gg <- graph_from_data_frame(d = as.data.frame(integrSIF[, c(1, 3)]),
                              directed = TRUE)
  adj <- get.adjacency(graph = gg)
  idx2keep <- c()
  for(ii in 1:length(indices[[1]])){
    currMeas <- names(indices[[1]])[ii]
    currCues <- colnames(cnolist@cues)[which(cnolist@cues[indices[[1]][[ii]][2], ]==1)]
      spList <- list()
      for(jj in length(currCues)){
        sP <- all_simple_paths(graph = gg, 
                               from = which(rownames(adj)==currCues[jj]), 
                               to = which(rownames(adj)==currMeas))
        if(length(sP)>0){
          for(kk in 1:length(sP)){
            if(length(sP[[kk]])>1){
              for(ll in 1:(length(sP[[kk]])-1)){
                ss <- rownames(adj)[sP[[kk]]][ll]
                tt <- rownames(adj)[sP[[kk]]][ll+1]
                idx1 <- which(reacTable[, 1]==ss)
                idx2 <- which(reacTable[, 2]==tt)
                idx <- intersect(x = idx1, y = idx2)
                if(length(idx)>0){idx2keep <- c(idx2keep, idx)}
              }
            }
          }
        }
      }
  }
  
  idx2keep <- unique(idx2keep)
  if(length(idx2keep)>0){
    if(length(idx2keep)==1){
      reacTable <- t(as.matrix(reacTable[idx2keep, ]))
    } else {
      reacTable <- reacTable[idx2keep, ]
    }
  }
  
  if(!is.null(database)){
    
    ##
    # Transform path lengs from based on interactions to based on nodes
    if((pathLength < 1) || !is.numeric(pathLength)){
      stop("Path length should be numeric and equal or greater to 1")
    } else {
      pathLength = pathLength + 1
    }
    
    ##
    # Compacting cnolist
    if(class(cnolist)=="CNOlist"){
      cnolist = compatCNOlist(object = cnolist)
    }
    
    ##
    # Identifying the possible links
    modelSIF <- model2sif(model = model)
    if(length(indices$indices) > 0){
      indices <- indices$indices
    } else {
      indices <- NULL
    }
    
    df <- as.data.frame(x = database[, c(1, 3)])
    gg <- graph_from_data_frame(d = df, directed = TRUE)
    adj <- get.adjacency(gg)
    
    # all shortest paths connecting the measurements with the perturbed cues in the indices list
    sP_all <- list()
    if(!is.null(pathLength) && !is.null(indices)){
      
      for(ii in 1:length(indices)){
        
        measurement <- cnolist$namesSignals[indices[[ii]][1]]
        
        inhSet <- cnolist$namesCues[which(cnolist$valueCues[indices[[ii]][2], ]==1)]
        
        if(length(inhSet) > 0){
          
          for(jj in 1:length(inhSet)){
            
            if((inhSet[jj]%in%rownames(adj)) && (measurement%in%rownames(adj))){
              
              sP <- get.all.shortest.paths(graph = gg, from = which(rownames(adj)==inhSet[jj]), to = which(rownames(adj)==measurement))
              
              if(length(sP[[1]]) > 0){
                
                if(length(sP[[1]][[1]]) <= pathLength){
                  
                  for(kk in 1:length(sP[[1]])){
                    
                    sP_all[[length(sP_all)+1]] <- sP[[1]][[kk]]
                    
                  }
                  
                }
                
              }
              
            }
            
          }
          
        }
        
      }
      
    }
    
    if(length(sP_all)==0){
      stop("There cannot be found interactions from the database for the used settings. Please either change the settings or set database = NULL")
    }
    
    ##
    # now creating the feedInteractions list containing the signed interactions between the species
    feedInteractions <- list()
    for(ii in 1:length(sP_all)){
      
      feederMatrix <- matrix(data = , nrow = 1, ncol = 3)
      if(length(sP_all[[ii]]) > 1){
        
        for(jj in 1:(length(sP_all[[ii]])-1)){
          
          ss <- rownames(adj)[sP_all[[ii]][jj]]
          tt <- rownames(adj)[sP_all[[ii]][jj+1]]
          
          idx <- intersect(x = which(modelSIF[, 1]==ss), y = which(modelSIF[, 3]==tt))
          
          if(length(idx)==0){
            
            cc <- t(as.matrix(c(ss, database[intersect(which(database[, 1]==ss), which(database[, 3]==tt))[1], 2], tt)))
            
            feederMatrix <- unique(rbind(feederMatrix, cc))
            
          } else {
            
            cc <- t(as.matrix(c(ss, modelSIF[intersect(which(modelSIF[, 1]==ss), which(modelSIF[, 3]==tt))[1], 2], tt)))
            
            feederMatrix <- unique(rbind(feederMatrix, cc))
            
          }
          
        }
        
        if(nrow(feederMatrix)==2){
          
          feedInteractions[[length(feedInteractions)+1]] <- cc
          
        } else {
          
          feedInteractions[[length(feedInteractions)+1]] <- feederMatrix[-1, ]
          
        }
        
      }
      
    }
    
    ## 
    # removing interactions which involve self-activating interactions
    if(length(feedInteractions)<=0){
      
      stop("No links to add for these settings..")
      
    } else {
      
      feedInteractions <- unique(feedInteractions)
      idx2rem <- c()
      for(ii in 1:length(feedInteractions)){
        if(class(feedInteractions[[ii]])=="matrix"){
          if(((feedInteractions[[ii]][1, 1]%in%model$namesSpecies)==FALSE) || (feedInteractions[[ii]][1, 1]==feedInteractions[[ii]][1, 3])){
            idx2rem <- c(idx2rem, ii)
          }
        } else {
          idx2rem <- c(idx2rem, ii)
        }
      }
      
      if(length(idx2rem)>0){
        feedInteractions <- feedInteractions[-idx2rem]
      }
      
      ##
      # removing those feedInteraction cases which are already present in the PKN
      idx2rem <- c()
      for(ii in 1:length(feedInteractions)){
        
        cnt = 0
        for(jj in 1:nrow(feedInteractions[[ii]])){
          
          idx1 = which(modelSIF[, 1]==feedInteractions[[ii]][jj, 1])
          idx2 = which(modelSIF[, 3]==feedInteractions[[ii]][jj, 3])
          idx <- intersect(x = idx1, y = idx2)
          
          if(length(idx)>0){
            cnt <- cnt + 1
          }
          
        }
        
        if(cnt==nrow(feedInteractions[[ii]])){
          idx2rem <- c(idx2rem, ii)
          
        }
        
      }
      
      if(length(idx2rem)>0){
        feedInteractions <- feedInteractions[-idx2rem]
      }
      
      ##
      # removing mechanisms which might have incoming interactions in the stimuli's
      idx2rem = c()
      for(ii in 1:length(feedInteractions)){
        
        vv = unique(as.vector(feedInteractions[[ii]][, c(1, 3)]))[2:length(unique(as.vector(feedInteractions[[ii]][, c(1, 3)])))]
        
        int = intersect(x = vv, y = cnolist$namesStimuli)
        
        if(length(int) > 0){
          idx2rem = c(idx2rem, ii)
        }
        
      }
      
      if(length(idx2rem)>0){
        feedInteractions <- feedInteractions[-idx2rem]
      }
      
      ##
      # removing those feedInteractions cases which do not involve any cues
      idx2rem <- c()
      for(ii in 1:length(feedInteractions)){
        feedSpecies <- unique(c(feedInteractions[[ii]][, 1], feedInteractions[[ii]][, 3]))
        if(length(intersect(x = feedSpecies, y = cnolist$namesCues))==0){
          idx2rem <- c(idx2rem, ii)
        }
      }
      
      if(length(idx2rem)>0){
        feedInteractions <- feedInteractions[-idx2rem]
      }
      
      ##
      # keeping only those interactions which involve most of the current present species in the PKN
      idx2keep <- c()
      dM <- matrix(data = , nrow = 1, ncol = 3)
      for(ii in 1:length(feedInteractions)){
        
        curr <- feedInteractions[[ii]]
        ss <- curr[1, 1]
        tt <- curr[nrow(curr), 3]
        currSpecies <- unique(c(curr[, 1], curr[, 3]))
        
        dM <- rbind(dM, t(as.matrix(c(ss, tt, length(intersect(x = currSpecies, y = model$namesSpecies)))))) 
        
      }
      
      dM <- dM[-1, ]
      
      if(is.null(nrow(dM))){
        dM = t(as.matrix(dM))
      }
      
      for(ii in 1:nrow(dM)){
        
        idx1 <- which(dM[, 1]==dM[ii, 1])
        idx2 <- which(dM[, 2]==dM[ii, 2])
        idx <- intersect(x = idx1, idx2)
        
        if(dM[ii, 3]==max(dM[idx, 3])){
          
          idx2keep <- c(idx2keep, ii)
          
        }
        
      }
      
      feedInteractions <- feedInteractions[idx2keep]
      
      for(ii in 1:length(feedInteractions)){
        
        feedInteractions[[ii]][, 2] = as.character(as.numeric(feedInteractions[[ii]][, 2]))
        
      }
      
      ##
      # Now adding the links from the FEED method
      if(nrow(reacTable)>0){
        
        sifFEED = matrix(data = , nrow = nrow(reacTable), ncol = 3)
        sifFEED[, 1] = reacTable[, 1]
        sifFEED[, 3] = reacTable[, 2]
        sgn = c()
        
        for(ii in 1:nrow(reacTable)){
          
          idx1 = which(database[, 1]==reacTable[ii, 1])
          idx2 = which(database[, 3]==reacTable[ii, 2])
          
          idx = intersect(x = idx1, y = idx2)
          
          if(length(idx)==1){
            
            sgn = c(sgn, as.numeric(database[idx, 2]))
            
          } else {
            
            sgn = c(sgn, 1)
            
          }
          
        }
        
        sifFEED[, 2] = as.character(sgn)
        
      }
      
      if(DDN){
        
        feedInteractions[[length(feedInteractions)+1]] = sifFEED
        
      }
      
      if(length(feedInteractions)>0){
        
        ##
        # saving the object
        names(feedInteractions) <- paste0("Feed - ", 1:length(feedInteractions))
        
        object <- list()
        object[[length(object)+1]] <- modelSIF
        object[[length(object)+1]] <- feedInteractions
        
        names(object) <- c("Original PKN", "Feed mechanisms")
        
        return(object)
        
      } else {
        
        stop("No links to add for these settings..")
        
      }
      
    }
    
  } else {
    
    feedInteractions = list()
    if(nrow(reacTable)>0){
      
      sifFEED = matrix(data = , nrow = nrow(reacTable), ncol = 3)
      sifFEED[, 1] = reacTable[, 1]
      sifFEED[, 3] = reacTable[, 2]
      sgn = c()
      
      for(ii in 1:nrow(reacTable)){
        
        idx1 = which(database[, 1]==reacTable[ii, 1])
        idx2 = which(database[, 3]==reacTable[ii, 2])
        
        idx = intersect(x = idx1, y = idx2)
        
        if(length(idx)==1){
          
          sgn = c(sgn, as.numeric(database[idx, 2]))
          
        } else {
          
          sgn = c(sgn, 1)
          
        }
        
      }
      
      sifFEED[, 2] = as.character(sgn)
      
    }
    
    if(DDN){
      
      feedInteractions[[length(feedInteractions)+1]] = sifFEED
      
      if(length(feedInteractions)>0){
        
        ##
        # saving the object
        names(feedInteractions) <- paste0("Feed - ", 1:length(feedInteractions))
        
        object <- list()
        object[[length(object)+1]] <- model2sif(model = model)
        object[[length(object)+1]] <- feedInteractions
        
        names(object) <- c("Original PKN", "Feed mechanisms")
        
        return(object)
        
      } else {
        
        stop("No links to add for these settings..")
        
      }
      
    }
    
  }
  
}
