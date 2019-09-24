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
# during the feeding process -- MI method

# Inputs:
# Mandatory:  A cnolist object containing the data (cnolist)
#             A model to optimize (model)
#             A simulation data object for a specific set of dynamic parameters as returned by the getLBodeSimFunction.R function (simData)
#
# Optional:   A percentage of measurements identified as misfits based on the MI value (percMI = 0.1 or 10% by default).
#             A number of equally spced points to split the time-course data and observation (nSplines = 100 by default)
#             A binning number (nBins = 10 by default)
#             A spline interpolation method (method = "fmm" by default. For more check the spline() function)

computeMI <- function(cnolist = cnolist, model = model, simData = simData, percMI = 0.1, nSplines = 100, nBins = 10, spliningMethod = "fmm"){
  
  ##
  # Compacting cnolist
  if(class(cnolist)=="CNOlist"){
    cnolist = compatCNOlist(object = cnolist)
  }
  
  cnoSplines <- list()
  
  for(i in 1:ncol(cnolist$valueSignals[[1]])){
    
    temp <- list()
    
    for(j in 1:nrow(cnolist$valueSignals[[1]])){
      
      vals <- c()
      
      for(k in 1:length(cnolist$timeSignals)){
        
        vals <- c(vals, cnolist$valueSignals[[k]][j, i])
        
      }
      
      temp[[length(temp)+1]] <- spline(x = vals, y = NULL, n = nSplines + 1, method = spliningMethod, xmin = min(vals), xmax = max(vals))$x
      
    }
    
    cnoSplines[[length(cnoSplines)+1]] <- temp
    
  }
  
  names(cnoSplines) <- cnolist$namesSignals
  
  ##
  # simData=plotLBodesimData(cnolist, model, opt_pars, timeSignals=seq(0,cnolist$timeSignals[length(cnolist$timeSignals)],cnolist$timeSignals[length(cnolist$timeSignals)]/100))
  
  simSplines <- list()
  for(i in 1:length(cnolist$namesSignals)){
    
    # idx <- which(model$namesSpecies==cnolist$namesSignals[i])
    
    if(length(i) > 0){
      
      temp <- list()
      
      for(j in 1:nrow(simData[[1]])){
        
        vals <- c()
        
        for(k in 1:length(simData)){
          
          vals <- c(vals, simData[[k]][j, i])
          
        }
        
        # temp[[length(temp)+1]] <- vals
        temp[[length(temp)+1]] <- spline(x = vals, y = NULL, n = nSplines + 1, method = spliningMethod, xmin = min(vals), xmax = max(vals))$x
        
      }
      
      simSplines[[length(simSplines)+1]] <- temp
      
    }
    
  }
  
  ##
  cnoBin <- list()
  simBin <- list()
  
  for(i in 1:length(cnoSplines)){
    
    cnoTemp <- list()
    simTemp <- list()
    
    for(j in 1:length(cnoSplines[[i]])){
      
      cnoTemp[[length(cnoTemp)+1]] <- floor(nBins*cnoSplines[[i]][[j]])
      simTemp[[length(simTemp)+1]] <- floor(nBins*simSplines[[i]][[j]])
      
    }
    
    cnoBin[[length(cnoBin)+1]] <- cnoTemp
    simBin[[length(simBin)+1]] <- simTemp
    
  }
  
  ##
  matMI <- matrix(, nrow = length(cnoSplines), ncol = length(cnoSplines[[1]]))
  
  for(i in 1:nrow(matMI)){
    
    for(j in 1:ncol(matMI)){
      
      matMI[i, j] <- mutinformation(X = cnoBin[[i]][[j]], Y = simBin[[i]][[j]])
      
    }
    
  }
  
  matMI <- t(matMI)
  
  colnames(matMI) <- names(cnoSplines)
  
  idxZero <- which(matMI==0, arr.ind = TRUE)
  
  for(i in 1:nrow(idxZero)){
    
    if(all(cnoBin[[idxZero[i, 2]]][[idxZero[i, 1]]]==simBin[[idxZero[i, 2]]][[idxZero[i, 1]]])){
      
      matMI[idxZero[i, 1], idxZero[i, 2]] <- Inf
      
    }
    
  }
  
  ##
  misFitIdx <- which(x = matMI<=sort(x = matMI, decreasing = TRUE)[round(percMI*nrow(matMI)*ncol(matMI))], arr.ind = TRUE)
  miThresh <- max(matMI[misFitIdx])
  
  indices <- list()
  for(i in 1:nrow(matMI)){
    
    for(j in 1:ncol(matMI)){
      
      if((matMI[i, j] <= miThresh) && !is.na(matMI[i, j])){
        
        indices[[length(indices)+1]] <- c(round(j), round(i), matMI[i, j])
        names(indices)[length(indices)] <- names(matMI[i, j])
        names(indices[[length(indices)]]) <- NULL
        
      }
      
    }
    
  }
  
  ##
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
  
  idx <- indices
  
  indices <- list()
  indices[[length(indices)+1]] <- idx
  indices[[length(indices)+1]] <- matMI
  indices[[length(indices)+1]] <- mse
  
  names(indices) <- c("indices", "MI matrix", "mse")
  
  return(indices)
  
}