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

# This function evaluates the effects of possible feeder mechanisms which could potentially
# bee added to the PKN. A case study is model is build by binding feeder mechanisms with a
# subset of the original PKN linking all the cues set which have an influence over a specific
# measurement over which we are doing the evaluation and their corresponding directly linked
# measurements. A mechanism is only accepted if the overall fit of that measurement is improved
# by considering the feeder mechanism which we add and evaluate and which at the same time does
# not worsen the fit of the other measurements included in the case model within a specific 
# tolerance parameter.

# Inputs:
# Mandatory:  A cnolist object containing the data (cnolist)
#             A model to optimize (model)
#             A feeder object as returned by buildFeederObjectDynamic.R function

runDynamicFeeder <- function(cnolist = cnolist, integratedModel = integratedModel, ode_parameters = ode_parameters,
                             penFactor_k = 100, penFactor_tau = 1, penFactorPIN_k = 10,
                             paramsSSm = defaultParametersSSm()){
  
  ##
  # getLBodeContObjFunctionWeighted
  getLBodeContObjFunctionWeighted<-
    function
  (
    cnolist,                model,                    ode_parameters,
    indices=NULL,           time=1,                   verbose=0,
    transfer_function=3,    reltol=1e-4,            atol=1e-3,
    maxStepSize=Inf,        maxNumSteps=100000,        maxErrTestsFails=50,
    nan_fac=1, lambda_tau=0, lambda_k=0, bootstrap=F,
    SSpenalty_fac=0, SScontrolPenalty_fac=0, boot_seed=sample(1:10000,1)
  )
    {
    
    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
    adjMatrix=incidence2Adjacency(model);
    
    if(is.null(indices))indices=indexFinder(cnolist,model,verbose=FALSE);
    
    sim_function=getLBodeSimFunction(cnolist,model,adjMatrix1=adjMatrix,
                                     indices1=indices, odeParameters1=ode_parameters$parValues, time1=time,verbose1=verbose,
                                     transfer_function1=transfer_function,reltol1=reltol,atol1=atol,maxStepSize1=maxStepSize,
                                     maxNumSteps1=maxNumSteps,maxErrTestsFails1=maxErrTestsFails)
    
    logic_ode_continuous_objective_function<-function
    (
      x,                            cnolist1=cnolist,        model1=model,
      adjMatrix1=adjMatrix,        indices1=indices,        ode_parameters1=ode_parameters,
      sim_function1=sim_function,        nan_fac1=nan_fac
    )
    {
      
      
      ode_parameters1$parValues[ode_parameters1$index_opt_pars]=x;
      sim=sim_function1(cnolist1,model1,ode_parameters1$parValues);
      temp=sim;
      sim=list();
      for(i in 1:length(cnolist$timeSignals)){
        sim[[i]]=temp[[i]];
      }
      
      sim<-as.vector(unlist(lapply(sim,function(x) x[,indices1$signals])));
      measured_values<-as.vector(unlist(lapply(cnolist1$valueSignals,function(x)x)));
      
      # do bootstrap if required
      if (bootstrap==T){
        res_boot<-(sim-measured_values)^2
        set.seed(boot_seed)
        res_boot<-sample(res_boot, length(res_boot), replace=T)
        error<-sum(res_boot, na.rm = T)
        # print(sample(1:5, 5, replace=T))
      }else{
        error<-sum((sim-measured_values)^2, na.rm = T)
      }
      
      NaNs_sim=which(is.na(sim));
      
      
      NApenalty<-length(NaNs_sim)*nan_fac1
      SSpenalty<-SSpenalty_fac*sum(temp[[length(cnolist$timeSignals)+1]]^2)
      SScontrolPenalty<-SScontrolPenalty_fac*sum(temp[[length(cnolist$timeSignals)+2]][1,]^2)
      absVal = abs(ode_parameters1$parValues[ode_parameters1$index_k])
      L1reg_k=0
      for(ll in 1:length(lambda_k)){L1reg_k=L1reg_k+lambda_k[ll]*absVal[ll]}
      absVal = abs(ode_parameters1$parValues[ode_parameters1$index_tau])
      L1reg_tau=0
      for(ll in 1:length(lambda_tau)){L1reg_tau=L1reg_tau+lambda_tau[ll]*absVal[ll]}
      L1reg<-L1reg_k+L1reg_tau
      
      res=error+NApenalty+SSpenalty+SScontrolPenalty+L1reg;
      
      
      if (verbose==T){
        cat("res =", res, "\n",
            "NA penalty = ", NApenalty, "\n",
            "error =", error, "\n",
            "steady state penalty ( fac =", SSpenalty_fac, ") = ", SSpenalty, "\n",
            "control steady state penalty ( fac =", SScontrolPenalty_fac, ") = ",  SScontrolPenalty, "\n",
            "L1-reg: lambda =", lambda_tau, "- P penalty (only tau)", L1reg_tau, "\n",
            "L1-reg: lambda =", lambda_k, "- P penalty (only k)", L1reg_k, "\n",
            "L1-reg: P penalty (total)", L1reg, "\n",
            "--------\n")
      }
      
      if(is.nan(res) || is.na(res))res=1e10;
      
      return(res);
    }
    return(logic_ode_continuous_objective_function);
  }
  
  ##
  # parEstimationLBodeSSmWeighted
  parEstimationLBodeSSmWeighted <-function
  (
    cnolist,				model,					ode_parameters=NULL,
    indices=NULL,			maxeval=Inf,			maxtime=100,			
    ndiverse=NULL,			dim_refset=NULL, 		local_solver=NULL,      
    time=1,					verbose=0, 				transfer_function=3,	
    reltol=1e-4,			atol=1e-3,				maxStepSize=Inf,		
    maxNumSteps=100000,		maxErrTestsFails=50,	nan_fac=1,
    lambda_tau=0, lambda_k=0, bootstrap=F,
    SSpenalty_fac=0, SScontrolPenalty_fac=0, boot_seed=sample(1:10000,1)
  )
  {
    
    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
    if(!require(MEIGOR)) stop("MEIGOR (essR) package not found.
                              SSm not available. Install the MEIGOR package and load it or try the Genetic Algorithm
                              optimiser instead.");
    
    
    checkSignals(CNOlist=cnolist,model=model)
    
    adjMat=incidence2Adjacency(model);
    if(is.null(ode_parameters)){
      ode_parameters=createLBodeContPars(model,random=TRUE);
    }
    if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
    
    #Check if essR is installed
    dummy_f<-function(x){
      return(0);
    }
    problem<-list(f=dummy_f,x_L=rep(0),x_U=c(1));
    opts<-list();
    opts$maxeval=0;
    opts$maxtime=0;
    
    val=essR(problem,opts)
    
    problem=list();
    problem$f<-getLBodeContObjFunctionWeighted(cnolist=cnolist,
                                               model=model,
                                               ode_parameters=ode_parameters,
                                               indices=indices,
                                               time=time,
                                               verbose=verbose,
                                               transfer_function=transfer_function,
                                               reltol=reltol,
                                               atol=atol,
                                               maxStepSize=maxStepSize,
                                               maxNumSteps=maxNumSteps,
                                               maxErrTestsFails=maxErrTestsFails,
                                               nan_fac=nan_fac,
                                               lambda_tau=lambda_tau,
                                               lambda_k=lambda_k,
                                               bootstrap=bootstrap,
                                               SSpenalty_fac=SSpenalty_fac,
                                               SScontrolPenalty_fac=SScontrolPenalty_fac,
                                               boot_seed=boot_seed);
    problem$x_L <- ode_parameters$LB[ode_parameters$index_opt_pars];
    problem$x_U <- ode_parameters$UB[ode_parameters$index_opt_pars];
    problem$x_0<- ode_parameters$parValues[ode_parameters$index_opt_pars];
    problem$int_var =0;
    problem$bin_var =0;
    opts=list();
    opts$maxeval=maxeval;
    opts$maxtime=maxtime;
    if(!is.null(local_solver))opts$local_solver=local_solver;
    if(!is.null(ndiverse))opts$ndiverse=ndiverse;      
    if(!is.null(dim_refset))opts$dim_refset=dim_refset;  
    results=MEIGOR::essR(problem,opts);
    ode_parameters$parValues[ode_parameters$index_opt_pars]=results$xbest;
    ode_parameters$ssm_results=results;
    return(ode_parameters);	
  }
  
  ##
  # parEstimationLBodeWeighted
  parEstimationLBodeWeighted<-function (cnolist, model, method="ga",
                                        ode_parameters = NULL, indices = NULL, paramsGA=NULL, paramsSSm=NULL)
  {
    
    
    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
    if (method == "essm"){
      if (is.null(paramsSSm)){
        paramsSSm = defaultParametersSSm()
      }
      ode_parameters = parEstimationLBodeSSmWeighted(cnolist, model,
                                                     ode_parameters=ode_parameters, indices=indices,
                                                     maxeval=paramsSSm$maxeval,
                                                     maxtime=paramsSSm$maxtime,
                                                     ndiverse=paramsSSm$ndiverse,
                                                     dim_refset=paramsSSm$dim_refset,
                                                     local_solver=paramsSSm$local_solver,
                                                     time=paramsSSm$time,
                                                     verbose=paramsSSm$verbose,
                                                     transfer_function=paramsSSm$transfer_function,
                                                     reltol=paramsSSm$reltol,
                                                     atol=paramsSSm$atol,
                                                     maxStepSize=paramsSSm$maxStepSize,
                                                     maxNumSteps=paramsSSm$maxNumSteps,
                                                     maxErrTestsFails=paramsSSm$maxErrTestsFails,
                                                     nan_fac=paramsSSm$nan_fac,
                                                     
                                                     # added
                                                     lambda_tau=paramsSSm$lambda_tau,
                                                     lambda_k=paramsSSm$lambda_k,
                                                     bootstrap=paramsSSm$bootstrap,
                                                     SSpenalty_fac=paramsSSm$SSpenalty_fac,
                                                     SScontrolPenalty_fac=paramsSSm$SScontrolPenalty_fac,
                                                     boot_seed=paramsSSm$boot_seed
      )
      
      
    }
    else if(method=="ga"){
      if (is.null(paramsGA)){
        paramsGA = defaultParametersGA()
      }
      ode_parameters = parEstimationLBodeGA(cnolist, model,
                                            ode_parameters=ode_parameters,
                                            indices=indices,
                                            mutationChance=paramsGA$mutationChance,
                                            popSize=paramsGA$popSize,
                                            iters=paramsGA$iters,
                                            elitism=paramsGA$elitism,
                                            time=paramsGA$time,
                                            monitor=paramsGA$monitor,
                                            verbose=paramsGA$verbose,
                                            transfer_function=paramsGA$transfer_function,
                                            reltol=paramsGA$reltol,
                                            atol=paramsGA$atol,
                                            maxStepSize=paramsGA$maxStepSize,
                                            maxNumSteps=paramsGA$maxNumSteps,
                                            maxErrTestsFails=paramsGA$maxErrTestsFails,
                                            nan_fac=paramsGA$nan_fac)
    }
    else{
      stop ("method argument must be either 'ga' or 'essm'." )
    }
    return(ode_parameters)
  }
  
  ##
  # Modelling of the integrated network
  
  #
  model = integratedModel$model
  reacDiff = model$reacID[integratedModel$integLinksIdx]
  speciesDiff = model$namesSpecies[integratedModel$integSpeciesIdx]
  
  #
  lambda_tau = rep(paramsSSm$lambda_tau, length(ode_parameters$index_tau)) 
  for(ii in 1:length(ode_parameters$index_tau)){
    if(ode_parameters$parNames[ode_parameters$index_tau[ii]]%in%paste0("tau_", speciesDiff)){
      lambda_tau[ii] = penFactor_tau
    }
  }
  paramsSSm$lambda_tau = lambda_tau
  
  #
  lambda_k = rep(paramsSSm$lambda_k, length(ode_parameters$index_k))
  ll = strsplit(x = reacDiff, split = "=", fixed = TRUE)
  ss = unlist(lapply(ll, function(ll) ll[[1]]))
  ss = gsub(pattern = "!", replacement = "", x = ss, fixed = TRUE)
  tt = unlist(lapply(ll, function(ll) ll[[2]]))
  reactions = paste0(ss, "_k_", tt)
  for(ii in 1:length(ode_parameters$index_k)){
    if(ode_parameters$parNames[ode_parameters$index_k[ii]]%in%reactions){
      ss = strsplit(x = ode_parameters$parNames[ode_parameters$index_k[ii]], split = "_k_", fixed = TRUE)[[1]][1]
      tt = strsplit(x = ode_parameters$parNames[ode_parameters$index_k[ii]], split = "_k_", fixed = TRUE)[[1]][2]
      idx = which(model$reacID%in%c(paste0(ss, "=", tt), paste0("!", ss, "=", tt)))[1]
      if(is.infinite(integratedModel$databasePenalty[idx])){
        lambda_k[ii] = penFactor_k
      } else {
        lambda_k[ii] = penFactorPIN_k*(1+integratedModel$databasePenalty[idx])
      }
    }
  }
  paramsSSm$lambda_k = lambda_k
  
  #
  opt_pars = parEstimationLBodeWeighted(cnolist,model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm)
  
  res = list()
  res[[1]] <- opt_pars
  res[[2]] <- integratedModel
  res[[3]] <- cnolist
  
  names(res) <- c("Parameters", "Integrated-Model", "CNOList")
  
  return(res)
  
}
