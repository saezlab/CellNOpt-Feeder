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

# This function identifies poorly fitted measurements for specific experimental conditions.
# It returns a list of possible indices and mse's pointing to possible connections to be added
# during the feeding process -- MSE method

# Inputs:
# Mandatory:  A cnolist object containing the data (cnolist)
#             A model to optimize (model)
#             A simulation data object for a specific set of dynamic parameters as returned by the getLBodeSimFunction.R function (simData)
# Optional:   A thrreshold parameter for minimal misfit to be considered (mseThresh = 0.05 by default)

computeMSE <- function(cnolist = cnolist, model = model, mseThresh = 0.05, simData = simData){
  
  ##
  # Compacting cnolist
  if(class(cnolist)=="CNOlist"){
    cnolist = compatCNOlist(object = cnolist)
  }
  
  ##
  # Listing all the data-points for each measurement across each condition
  cnoSplines <- list()
  
  for(i in 1:ncol(cnolist$valueSignals[[1]])){
    
    temp <- list()
    
    for(j in 1:nrow(cnolist$valueSignals[[1]])){
      
      vals <- c()
      
      for(k in 1:length(cnolist$timeSignals)){
        
        vals <- c(vals, cnolist$valueSignals[[k]][j, i])
        
      }
      
      temp[[length(temp)+1]] <- vals
      
    }
    
    cnoSplines[[length(cnoSplines)+1]] <- temp
    
  }
  
  names(cnoSplines) <- cnolist$namesSignals
  
  ##
  # Listing all the simulated values for each measurement across each condition
  sim_data = simData
  
  simSplines <- list()
  
  for(i in 1:ncol(sim_data[[1]])){
    
    temp <- list()
    
    for(j in 1:nrow(sim_data[[1]])){
      
      cc <- sim_data[[1]][j, i]
      for(k in 2:length(cnoSplines[[1]][[1]])){
        
        cc <- c(cc, sim_data[[k]][j, i])
        
      }
      
      temp[[length(temp)+1]] <- cc
      
    }
    
    simSplines[[length(simSplines)+1]] <- temp
    
  }
  
  names(simSplines) <- cnolist$namesSignals
  
  ##
  # computing mse between simulations and data for each measurement at each condition
  mse <- matrix(data = , nrow = nrow(cnolist$valueSignals[[1]]), ncol = ncol(cnolist$valueSignals[[1]]))
  for(j in 1:length(simSplines)){
    
    for(i in 1:length(simSplines[[j]])){
      
      cc <- (cnoSplines[[j]][[i]][1]-simSplines[[j]][[i]][1])^2
      
      for(k in 2:length(simSplines[[j]][[i]])){
        
        cc <- c(cc, (cnoSplines[[j]][[i]][k]-simSplines[[j]][[i]][k])^2)
        
      }
      
      ss <- mean(cc)
      
      mse[i, j] <- ss
      
    }
    
  }
  
  colnames(mse) <- cnolist$namesSignals
  
  ##
  # identifying the list of indices indicating at which experiment a measurement is poorly fitted in comparison to the specified threshold
  indices <- list()
  for(i in 1:nrow(mse)){
    
    for(j in 1:ncol(mse)){
      
      if((mse[i, j] > mseThresh) && !is.na(mse[i, j])){
        
        indices[[length(indices)+1]] <- c(j, i, mse[i, j])
        
      }
      
    }
    
  }
  
  idx <- indices
  
  ##
  # returing a list containing the indices and the mse matrix
  if(length(idx) > 0){
    
    idxNames = c()
    for(ii in 1:length(idx)){
      
      idxNames = c(idxNames, names(idx[[ii]])[3])
      names(idx[[ii]]) = NULL
      
    }
    
    names(idx) = idxNames
    
    indices <- list()
    indices[[length(indices)+1]] <- idx
    indices[[length(indices)+1]] <- mse
    
    names(indices) <- c("indices", "mse")
    
    return(indices)
    
  } else {
    
    print("No measurement error falls within the specified error threshold specified")
    
    indices <- list()
    indices[[length(indices)+1]] <- idx
    indices[[length(indices)+1]] <- mse
    
    names(indices) <- c("indices", "mse")
    
    return(indices)
    
  }
  
}