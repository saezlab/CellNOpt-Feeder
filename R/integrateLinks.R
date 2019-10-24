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

# This function integrates the new links inferred via the FEED method or from the
# database to the original PKN

# Inputs:
# Mandatory:  A cnolist object containing the data (cnolist)
#             A model to optimize (model)
#             A feeder object as returned by buildFeederObjectDynamic.R function

integrateLinks <- function(feederObject = feederObject, cnolist = cnolist, database = NULL){
  
  # optimizing all the added interactions together
  object = feederObject
  initial_sif <- object$`Original PKN`
  sif <- object$`Original PKN`
  
  write.table(x = sif, file = "temporary_sif.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
  model <- readSIF(sifFile = "temporary_sif.txt")
  file.remove("temporary_sif.txt")
  
  for(ii in 1:length(object$`Feed mechanisms`)){
    
    sif <- unique(rbind(sif, object$`Feed mechanisms`[[ii]]))
    
  }
  
  write.table(x = sif, file = "temporary_sif.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
  currModel <- readSIF(sifFile = "temporary_sif.txt")
  file.remove("temporary_sif.txt")
  
  caseCNO <- cnolist
  
  if(is.null(database) || ncol(database)<=3){
    
    currModel <- preprocessing(data = caseCNO, model = currModel, compression = FALSE, expansion = FALSE)
    curr_sif = model2sif(model = currModel)
    idx2rem = which(duplicated(x = curr_sif[, c(1, 3)]))
    if(length(idx2rem) > 0){
      curr_sif = curr_sif[-idx2rem, ]
    }
    write.table(x = curr_sif, file = "temporary_sif.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
    currModel <- readSIF(sifFile = "temporary_sif.txt")
    file.remove("temporary_sif.txt")
    
    reacDiff <- setdiff(currModel$reacID, model$reacID)
    speciesDiff <- setdiff(currModel$namesSpecies, model$namesSpecies)
    
    reacDiffIdx = which(currModel$reacID%in%reacDiff)
    speciesDiffIdx = which(currModel$namesSpecies%in%speciesDiff)
    
    returnModel = list()
    
    returnModel[[length(returnModel)+1]] = currModel
    returnModel[[length(returnModel)+1]] = reacDiffIdx
    returnModel[[length(returnModel)+1]] = speciesDiffIdx
    returnModel[[length(returnModel)+1]] = rep(0, length(currModel$reacID))
    
    names(returnModel) = c("model", "integLinksIdx", "integSpeciesIdx", "databaseWeight")
    
    return(returnModel)
    
  } else {
    
    database[, 2] = as.character(as.numeric(database[, 2]))
    database[, 4] = as.character(as.numeric(database[, 4]))
    
    reacDiff <- setdiff(currModel$reacID, model$reacID)
    speciesDiff <- setdiff(currModel$namesSpecies, model$namesSpecies)
    
    weights = rep(0, length(currModel$reacID))
    weightsID = currModel$reacID
    
    for(ii in 1:length(reacDiff)){
      
      if(grepl(pattern = "!", x = reacDiff[ii], fixed = TRUE)){
        
        currReac = gsub(pattern = "!", replacement = "", x = reacDiff[ii], fixed = TRUE)
        ss = strsplit(x = currReac, split = "=", fixed = TRUE)[[1]][1]
        tt = strsplit(x = currReac, split = "=", fixed = TRUE)[[1]][2]
        sgn = "-1"
        
        idx = intersect(x = intersect(x = which(database[, 1]==ss), y = which(database[, 3]==tt)), y = which(database[, 2]==sgn))[1]
        
        if(length(idx) > 0){
          weights[which(weightsID==reacDiff[ii])] = 1 - as.numeric(database[idx, 4])
        } else {
          weights[which(weightsID==reacDiff[ii])] = Inf
        }
        
      } else {
        
        currReac = reacDiff[ii]
        ss = strsplit(x = currReac, split = "=", fixed = TRUE)[[1]][1]
        tt = strsplit(x = currReac, split = "=", fixed = TRUE)[[1]][2]
        sgn = "1"
        
        idx = intersect(x = intersect(x = which(database[, 1]==ss), y = which(database[, 3]==tt)), y = which(database[, 2]==sgn))
        
        if(length(idx)>0){
          weights[which(weightsID==reacDiff[ii])] = 1 - as.numeric(database[idx, 4])
        } else {
          weights[which(weightsID==reacDiff[ii])] = Inf
        }
        
      }
      
    }
    
    currModel = preprocessingWeighted(data = caseCNO, model = currModel, compression = FALSE, expansion = FALSE, weights = weights, weightsID = weightsID)
    
    returnModel = list()
    temp = currModel$model
    returnModel[[length(returnModel)+1]] = temp
    returnModel[[length(returnModel)+1]] = which(temp$reacID%in%setdiff(temp$reacID, model$reacID))
    returnModel[[length(returnModel)+1]] = which(temp$namesSpecies%in%setdiff(temp$namesSpecies, model$namesSpecies))
    returnModel[[length(returnModel)+1]] = currModel$weights
    
    names(returnModel) = c("model", "integLinksIdx", "integSpeciesIdx", "databaseWeight")
    
    return(returnModel)
    
  }
  
}
