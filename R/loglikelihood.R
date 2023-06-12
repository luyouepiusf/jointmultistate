##' Calculate the log-likelihood
##'
##' @title Calculate the log-likelihood
##' @param model_fit a model fit returned by \code{joint_multistate_model},
##' @param nqmc an integer value. \cr 
##' Number of quasi Monte Carlo points for evaluating numeric integration.
##' @return a list containing the log-likelihood, AIC, and BIC.
loglikelihood<-function(model_fit,nqmc=2000){
  
  result_est<-model_fit$result_est
  yyikj<-model_fit$data$yyikj
  ttikj<-model_fit$data$ttikj
  nind<-model_fit$data$nind
  ndimy<-model_fit$data$ndimy
  nobsik<-model_fit$data$nobsik
  statebegin<-model_fit$data$statebegin
  stateend<-model_fit$data$stateend
  timebegin<-model_fit$data$timebegin
  timeend<-model_fit$data$timeend
  ntransition<-model_fit$data$yyikj
  possible_transitions<-model_fit$data$possible_transitions
  nstates<-model_fit$data$nstates
  zzi<-model_fit$data$zzi
  ndimz<-model_fit$data$ndimz
  xxi<-model_fit$data$xxi
  ndimx<-model_fit$data$ndimx
  basismatrix_long<-model_fit$data$basismatrix_long
  nbasy<-model_fit$data$nbasy
  ntimepoints<-model_fit$data$ntimepoints
  
  alleeieval<-t(randtoolbox::sobol(nqmc,dim=sum(nbasy),normal=TRUE))
  
  result_lik<-lik_covariate(
    yyikj=yyikj,ttikj=ttikj,nind=nind,ndimy=ndimy,nobsik=nobsik,
    statebegin=statebegin-1,stateend=stateend-1,timebegin=timebegin,timeend=timeend,
    ntransition=ntransition,possible_transitions=possible_transitions+0,nstates=nstates,
    zzi=zzi,ndimz=ndimz,xxi=xxi,ndimx=ndimx,
    basismatrix_long=basismatrix_long,nbas_long=nbasy,alleeieval=alleeieval,
    ntimepoints=ntimepoints,result_est)
  
  nparameter<-
    ndimz*sum(possible_transitions)+
    ndimy*sum(possible_transitions)+
    sum(possible_transitions)-sum(rowSums(possible_transitions)>=1)+
    ndimy*ndimx+
    sum(nbasy)+
    (sum(nbasy))^2
  
  loglikelihood<-result_lik[[1]]
  AIC<-2*nparameter-2*loglikelihood
  BIC<-log(nind)*nparameter-2*loglikelihood
  
  return(list(
    loglikelihood=loglikelihood,
    AIC=AIC,
    BIC=BIC))
}
