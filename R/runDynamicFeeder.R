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

# This function evaluates the effects of possible feeder mechanisms which could 
# potentially be added to the PKN. A case study is model is build by binding 
# feeder mechanisms with a subset of the original PKN linking all the cues set 
# which have an influence over a specific measurement over which we are doing 
# the evaluation and their corresponding directly linked measurements. A 
# mechanism is only accepted if the overall fit of that measurement is improved
# by considering the feeder mechanism which we add and evaluate and which at the
# same time does not worsen the fit of the other measurements included in the 
# case model within a specific tolerance parameter.

# Inputs:
# Mandatory:  A cnolist object containing the data (cnolist)
#             A model to optimize (model)
#             A feeder object as returned by buildFeederObjectDynamic.R function

runDynamicFeeder <- function(cnolist = cnolist, 
                             integratedModel = integratedModel, 
                             ode_parameters = ode_parameters,
                             penFactor_k = 100, penFactor_tau = 1, 
                             penFactorPIN_k = 10,
                             paramsSSm = defaultParametersSSm()){
  
  ##
  model = integratedModel$model
  reacDiff = model$reacID[integratedModel$integLinksIdx]
  speciesDiff = model$namesSpecies[integratedModel$integSpeciesIdx]
  
  ##
  lambda_tau = rep(paramsSSm$lambda_tau, length(ode_parameters$index_tau)) 
  for(ii in 1:length(ode_parameters$index_tau)){
    if(ode_parameters$parNames[ode_parameters$index_tau[ii]]%in%
       paste0("tau_", speciesDiff)){
      lambda_tau[ii] = penFactor_tau
    }
  }
  paramsSSm$lambda_tau = lambda_tau
  
  ##
  lambda_k = rep(paramsSSm$lambda_k, length(ode_parameters$index_k))
  ll = strsplit(x = reacDiff, split = "=", fixed = TRUE)
  ss = unlist(lapply(ll, function(ll) ll[[1]]))
  ss = gsub(pattern = "!", replacement = "", x = ss, fixed = TRUE)
  tt = unlist(lapply(ll, function(ll) ll[[2]]))
  reactions = paste0(ss, "_k_", tt)
  for(ii in 1:length(ode_parameters$index_k)){
    if(ode_parameters$parNames[ode_parameters$index_k[ii]]%in%reactions){
      ss = strsplit(x = ode_parameters$parNames[ode_parameters$index_k[ii]], 
                    split = "_k_", fixed = TRUE)[[1]][1]
      tt = strsplit(x = ode_parameters$parNames[ode_parameters$index_k[ii]], 
                    split = "_k_", fixed = TRUE)[[1]][2]
      idx = which(model$reacID%in%c(paste0(ss, "=", tt), 
                                    paste0("!", ss, "=", tt)))[1]
      if(is.infinite(integratedModel$databaseWeight[idx])){
        lambda_k[ii] = penFactor_k
      } else {
        lambda_k[ii] = penFactorPIN_k*(1+integratedModel$databaseWeight[idx])
      }
    }
  }
  paramsSSm$lambda_k = lambda_k
  
  #
  opt_pars = parEstimationLBodeWeighted(cnolist,model, method="essm", 
                                        ode_parameters=ode_parameters, 
                                        paramsSSm=paramsSSm)
  
  res = list()
  res[[1]] <- opt_pars
  res[[2]] <- integratedModel
  res[[3]] <- cnolist
  
  names(res) <- c("Parameters", "Integrated-Model", "CNOList")
  
  return(res)
  
}
