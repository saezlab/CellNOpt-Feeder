#
#  This file is part of the CNO software
#
#  Copyright (c) Prof. Julio Saez-Rodriguez group - JRC COMBINE, RWTH Aachen
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
## Enio Gjerga - 17/01/2017

preprocessingWeighted<-function(data=NULL, model, cutNONC=TRUE, compression=TRUE,
                                expansion=TRUE, ignoreList=NA, maxInputsPerGate=2,verbose=FALSE, weights, weightsID){
  
  # why not doing this check here ? Does not cost too much
  if (is.null(data)!=TRUE){
    checkSignals(CNOlist=data,model=model)
  }
  
  pknmodel = model
  
  if(is.null(weights) && is.null(weightsID)){
    print("\nWeights not given. Performing a simple preprocessing.\n\n")
    return(preprocessing(data, model, cutNONC, compression = compression, expansion = expansion, ignoreList = ignoreList,
                         maxInputsPerGate = maxInputsPerGate, verbose = verbose))
  }
  
  if((length(weights)==length(weightsID)) == FALSE){
    stop("\nWeights and weightsID should have the same array length")
  }
  
  if(class(weightsID) == "character"){
    idTemp <- as.numeric(matrix(0, nrow = 1, ncol = length(weightsID)))
    for(i in 1:length(weightsID)){
      idTemp[i] = which(pknmodel$reacID==weightsID[i])
    }
    
    weightsID <- idTemp
  }
  
  if(length(weights)==length(pknmodel$reacID) && length(weightsID)==length(pknmodel$reacID)){
    weights1 <- weights
    weightsID1 <- weightsID
  }
  
  if(length(weights) < length(pknmodel$reacID)){
    #weights1 <- as.numeric(matrix(0.5, 1, length(pknmodel$reacID)))
    weights1 <- rep(1, length(pknmodel$reacID))
    for(i in 1:length(weightsID)){
      weights1[weightsID[i]] <- weights[i]
    }
  }
  
  ww1 <- as.data.frame(weights1)
  rownames(ww1) <- pknmodel$reacID
  weightsID1 <- c(1:length(pknmodel$reacID))
  
  # a copy of the model
  cutModel <- model
  
  if (cutNONC==TRUE && is.null(data)!=TRUE){
    # Find the indices, in the model, of the species that are inh/stim/sign
    indices<-indexFinder(CNOlist=data, model=model,    verbose=verbose)
    
    # Find the indices of the non-osb/non-contr
    temp_indices <- findNONC(model=model, indexes=indices,verbose=verbose)
    # Cut the nonc off the model
    cutModel <-cutNONC(model=model, NONCindexes=temp_indices)
  }
  
  if (compression == TRUE && is.null(data)!=TRUE){
    # Recompute the indices
    temp_indices<-indexFinder(CNOlist=data, model=cutModel)
    
    # Compress the model
    cutModel<-compressModel(model=cutModel,indexes=temp_indices)
    
    #################
    # set the weights after compression
    
    diffCut <- setdiff(cutModel$reacID, pknmodel$reacID)
    diffPKN <- setdiff(pknmodel$reacID, cutModel$reacID)
    int <- intersect(pknmodel$reacID, cutModel$reacID)
    
    
    weights2 <- matrix(, 1, length(cutModel$reacID))
    weightsID2 <- c(1:length(cutModel$reacID))
    for(i in 1:length(int)){
      for(j  in 1:length(cutModel$reacID)){
        if(int[i] == cutModel$reacID[j]){
          weights2[1, j] <- weights1[which(pknmodel$reacID == int[i])]
        }
      }
    }
    weights2 <- as.numeric(weights2)
    ww2 <- as.data.frame(weights2)
    rownames(ww2) <- cutModel$reacID
    
    ###
    
    SIF1 <- model2sif(pknmodel)
    g1 <- sif2graph(sif = SIF1)
    pknGraph <- graph_from_graphnel(g1)
    
    SIF2 <- model2sif(cutModel)
    g2 <- sif2graph(sif = SIF2)
    cutGraph <- graph_from_graphnel(g2)
    
    for(i in 1:length(diffCut)){
      kk <- which(colnames(cutModel$interMat) == diffCut[i])
      source <- row.names(cutModel$interMat)[which(cutModel$interMat[, kk]==-1)]
      target <- row.names(cutModel$interMat)[which(cutModel$interMat[, kk]==1)]
      
      idSpkn <- which(row.names(get.adjacency(pknGraph)) == source)
      idTpkn <- which(row.names(get.adjacency(pknGraph)) == target)
      
      idScut <- which(row.names(get.adjacency(cutGraph)) == source)
      idTcut <- which(row.names(get.adjacency(cutGraph)) == target)
      
      
      allPKN <- all_simple_paths(pknGraph, from = idSpkn, to = idTpkn)
      allCUT <- all_simple_paths(cutGraph, from = idScut, to = idTcut)
      #diffPath <- setdiff(allPKN, allCUT)
      simplePaths <- list()
      for(j in 1:length(allPKN)){
        pp <- as.numeric(unlist(allPKN[[j]]))
        flag = 0
        for(v in 1:length(allCUT)){
          cc <- as.numeric(unlist(allCUT[[v]]))
          if(identical(rownames(pknmodel$interMat)[pp], rownames(cutModel$interMat)[cc])){
            flag <- flag + 1
          }
        }
        if(flag == 0){
          simplePaths[[length(simplePaths)+1]] <- allPKN[[j]]
        }
      }
      
      mins <- list()
      for(j in 1:length(simplePaths)){
        mins[[length(mins)+1]] <- c(1)
        sp <- unlist(simplePaths[[j]])
        for(tt in 1:(length(sp)-1)){
          for(zz in (tt+1):length(sp)){
            eID <- intersect(grep(rownames(pknmodel$interMat)[unlist(sp)[tt]], pknmodel$reacID), 
                             grep(rownames(pknmodel$interMat)[unlist(sp)[zz]], pknmodel$reacID))
            
            mins[[j]] <- c(mins[[j]], weights1[eID])
          }
          #mins <- c(mins, weights1[eID])
        }
      }
      #weights2[kk] <- min(mins)
      if(length(simplePaths)==1){
        weights2[kk] <- prod(mins[[1]], na.rm = TRUE)
      }
      else{
        ss = 1
        for(mm in 1:length(mins)){
          ss <- ss*(1-prod(mins[[mm]], na.rm = TRUE)) 
        }
        weights2[kk] <- 1-ss
      }
    }
    
    ww2 <- as.data.frame(weights2)
    rownames(ww2) <- cutModel$reacID
    
    ww1 <- as.data.frame(weights1)
    rownames(ww1) <- pknmodel$reacID
    
    if(verbose == TRUE){
      print("Weights values before compression:")
      print(ww1)
      print("Weight values after compression:")
      print(ww2)
    }
    
    www <- weights2
    wwwID <- as.numeric(1:length(cutModel$reacID))
    
    procModel <- cutModel
  }
  
  ######### 
  # Expand the gates and set the weights after expansion
  if (expansion == TRUE){
    expModel <- expandGates(model=cutModel, ignoreList=ignoreList,maxInputsPerGate=maxInputsPerGate)
    
    int <- intersect(cutModel$reacID, expModel$reacID)
    
    weights3 <- matrix(, 1, length(expModel$reacID))
    weightsID3 <- c(1:length(expModel$reacID))
    for(i in 1:length(int)){
      for(j  in 1:length(expModel$reacID)){
        if(int[i] == expModel$reacID[j]){
          weights3[1, j] <- weights2[which(cutModel$reacID == int[i])]
        }
      }
    }
    
    addAnd <- c()
    andReac <- expModel$reacID[is.na(weights3[1, ])]
    for(i in 1:length(andReac)){
      andID <- which(expModel$reacID==andReac[i])
      sources <- rownames(expModel$interMat)[which(expModel$interMat[, andID] == -1)]
      target <- rownames(expModel$interMat)[which(expModel$interMat[, andID] == 1)]
      
      andMins <- c()
      for(j in 1:length(sources)){
        andMins <- c(andMins, weights2[intersect(grep(sources[j], cutModel$reacID), grep(target, cutModel$reacID))])
      }
      addAnd <- c(addAnd, prod(andMins, na.rm = TRUE))
    }
    
    weights3 <- as.numeric(weights3)
    weights3[is.na(weights3)] <- addAnd
    
    #weights3[is.na(weights3)] <- 0.5
    #weights3 <- as.numeric(weights3)
    ww3 <- as.data.frame(weights3)
    rownames(ww3) <- expModel$reacID
    
    if(verbose == TRUE){
      print("Weight values after expansion:")
      print(ww3)
    }
    
    www <- weights3
    wwwID <- as.numeric(1:length(expModel$reacID))
    
    procModel <- expModel
  }
  
  if(compression==FALSE && expansion==FALSE){
    
    procModel = cutModel
    www = weights
    
  }
  
  wwwID = procModel$reacID
  
  sif = model2sif(model = procModel)
  idx2rem = which(duplicated(x = sif[, c(1, 3)]))
  if(length(idx2rem) > 0){
    sif = sif[-idx2rem, ]
    www = www[-idx2rem]
    wwwID = wwwID[-idx2rem]
    
  }
  write.table(x = sif, file = "temporary_sif.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
  procModel <- readSIF(sifFile = "temporary_sif.txt")
  file.remove("temporary_sif.txt")
  
  
  # since version 1.3.28 return only model, indices are recomputed in other
  # functions
  
  res = list()
  res[[length(res)+1]] = procModel
  res[[length(res)+1]] = www
  res[[length(res)+1]] = wwwID
  
  names(res) = c("model", "weights", "weightsID")
  
  return(res)
}