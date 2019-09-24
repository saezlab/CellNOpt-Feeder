aicCNO <- function(model = model, cnolist = cnolist, opt_pars = opt_pars, transfer_function = 4){
  
  simData = plotLBodeFitness(cnolist = cnolist, model = model, ode_parameters = opt_pars, transfer_function = transfer_function)
  dev.off()
  
  netSize = length(model$reacID)
  
  rssInitial = 0
  nInitial = 0
  for(ii in 1:length(cnolist@signals)){
    
    for(jj in 1:nrow(cnolist@signals[[ii]])){
      
      for(kk in 1:ncol(cnolist@signals[[ii]])){
        
        rssInitial <- rssInitial + (cnolist@signals[[ii]][jj, kk] - simData[[ii]][jj, kk])^2
        nInitial = nInitial+1
        
      }
      
    }
    
  }
  # kInitial <- length(opt_pars$index_k) + length(opt_pars$index_tau)
  kInitial = length(c(opt_pars$index_k, opt_pars$index_tau)) - length(which(opt_pars$parValues[opt_pars$index_k]==0)) - length(which(opt_pars$parValues[opt_pars$index_tau]==0))
  
  aic = 2*kInitial + nInitial*log(x = rssInitial/nInitial)
  names(aic) = NULL
  aicc = aic + (2*kInitial^2 + 2*kInitial)/(nInitial - kInitial - 1)
  names(aicc) = NULL
  bic = nInitial*log(rssInitial/nInitial) + kInitial*log(nInitial)
  names(bic) = NULL
  
  res = list()
  res[[length(res)+1]] = aic
  res[[length(res)+1]] = aicc
  res[[length(res)+1]] = bic
  
  names(res) = c("AIC", "AICc", "BIC")
  
  return(res)
  
}