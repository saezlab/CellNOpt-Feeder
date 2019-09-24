map2cys <- function(model = model, cnolist = cnolist, opt_pars = opt_pars){
  
  ##
  # Some initial checks
  if(!(all(colnames(cnolist@signals[[1]])%in%model$namesSpecies) && all(colnames(cnolist@cues%in%model$namesSpecies)))){
    stop("model and cnolist objects do not match !!")
  }
  
  modelSIF = model2sif(model = model)
  if(!(all(opt_pars$parNames[opt_pars$index_k]%in%paste0(modelSIF[, 1], "_k_", modelSIF[, 3])) && 
       all(opt_pars$parNames[opt_pars$index_n]%in%paste0(modelSIF[, 1], "_n_", modelSIF[, 3])) &&
       all(opt_pars$parNames[opt_pars$index_tau]%in%paste0("tau_", model$namesSpecies)))){
    stop("model and opt_pars objects do not match !!")
  }
  
  ##
  # Assigning edge attributes
  edgeMatrix = matrix(data = , nrow = length(model$reacID), ncol = 5)
  colnames(edgeMatrix) = c("Source", "Sign", "Target", "k_param", "n_param")
  edgeMatrix[, 1] = modelSIF[, 1]
  edgeMatrix[, 2] = modelSIF[, 2]
  edgeMatrix[, 3] = modelSIF[, 3]
  
  for(ii in 1:nrow(edgeMatrix)){
    edgeMatrix[ii, 4] = opt_pars$parValues[which(opt_pars$parNames==paste0(edgeMatrix[ii, 1], "_k_", edgeMatrix[ii, 3]))]
    edgeMatrix[ii, 5] = opt_pars$parValues[which(opt_pars$parNames==paste0(edgeMatrix[ii, 1], "_n_", edgeMatrix[ii, 3]))]
  }
  
  ##
  # Assigning node attributes
  nodeMatrix = matrix(data = , nrow = length(model$namesSpecies), ncol = 3)
  colnames(nodeMatrix) = c("Node", "Status", "tau_param")
  nodeMatrix[, 1] = model$namesSpecies
  nodeMatrix[, 2] = "Unmeasured"
  nodeMatrix[, 3] = "1"
  
  for(ii in 1:length(model$namesSpecies)){
    
    if(model$namesSpecies[ii]%in%colnames(cnolist@inhibitors)){
      nodeMatrix[ii, 2] = "Inhibited"
    }
    
    if(model$namesSpecies[ii]%in%colnames(cnolist@stimuli)){
      nodeMatrix[ii, 2] = "Stimulated"
    }
    
    if(model$namesSpecies[ii]%in%colnames(cnolist@signals[[1]])){
      if(model$namesSpecies[ii]%in%colnames(cnolist@inhibitors)){
        nodeMatrix[ii, 2] = "Inh&Meas"
      } else {
        nodeMatrix[ii, 2] = "Measured"
      }
    }
    
    if(paste0("tau_", model$namesSpecies[ii])%in%opt_pars$parNames){
      nodeMatrix[ii, 3] = as.character(opt_pars$parValues[which(opt_pars$parNames==paste0("tau_", nodeMatrix[ii, 1]))])
    }
    
  }
  
  ##
  # Returning object
  
  attributes = list()
  attributes[[length(attributes) + 1]] = edgeMatrix
  attributes[[length(attributes) + 1]] = nodeMatrix
  
  names(attributes) = c("Edge Attributes", "Node Attributes")
  
  return(attributes)
  
}