##' Fit a joint model of multivariate longitudinal data and multistate data
##'
##' @title Fit a joint model of multivariate longitudinal data and multistate data
##' @description The function \code{joint_multistate_model} fits 
##' a joint model of multivariate longitudinal data and multistate data
##' @param yyikj 3-dimensional numeric array of longitudinal data. \cr
##' \code{yyikj[i,k,j]} is the jth observation of the kth longitudinal measure of the ith subject.
##' @param ttikj 3-dimensional integer array of longitudinal observation times. \cr
##' \code{ttikj[i,k,j]} is the jth longitudinal observation time of the kth longitudinal measure of the ith subject,
##' (i.e., the time of observing \code{yyikj[i,k,j]}).
##' @param nobsik number of observations arranged as an integer matrix. \cr
##' \code{nobsik[i]} is the number of observations for the kth longitudinal measure of the ith subject.
##' @param zzi a numeric matrix of covariates in the multistate model. \cr
##' \code{zzi[i,k]} is the kth covariate of the ith subject.
##' @param xxi a numeric matrix of covariates in the longitudinal model. \cr
##' \code{xxi[i,k]} is the kth covariate of the ith subject.
##' @param statebegin an integer matrix of starting states for the event-history intervals. \cr
##' The states should be encoded as consecutive integers ranging from 1 to the total number of states and 0 represents censoring. \cr
##' For example, if the event history for the ith subject is \cr
##' \eqn{[t_1,t_2)} in state \eqn{s_1}, \cr
##' \eqn{[t_2,t_3)} in state \eqn{s_2}, \cr
##' \eqn{[t_3,t_4)} in state \eqn{s_3}, \cr
##' and censored after \eqn{t_4}, then \cr
##' \code{statebegin[i,1]} is \eqn{s_1}, \cr
##' \code{statebegin[i,2]} is \eqn{s_2}, \cr
##' \code{statebegin[i,3]} is \eqn{s_3}, \cr
##' @param stateend an integer matrix of destination states for the event-history intervals. \cr
##' The states should be encoded as consecutive integers ranging from 1 to the total number of states and 0 represents censoring. \cr
##' For example, if the event history for the ith subject is \cr
##' \eqn{[t_1,t_2)} in state \eqn{s_1}, \cr
##' \eqn{[t_2,t_3)} in state \eqn{s_2}, \cr
##' \eqn{[t_3,t_4)} in state \eqn{s_3}, \cr
##' and censored after \eqn{t_4}, then \cr
##' \code{stateend[i,1]} is \eqn{s_2}, \cr
##' \code{stateend[i,2]} is \eqn{s_3}, \cr
##' \code{stateend[i,3]} is 0 (0 indicates censoring state). \cr
##' @param timebegin an integer matrix of lower bounds for the event-history intervals. \cr
##' For example, if the event history for the ith subject is \cr
##' \eqn{[t_1,t_2)} in state \eqn{s_1}, \cr
##' \eqn{[t_2,t_3)} in state \eqn{s_2}, \cr
##' \eqn{[t_3,t_4)} in state \eqn{s_3}, \cr
##' and censored after \eqn{t_4}, then \cr
##' \code{timebegin[i,1]} is \eqn{t_1}, \cr
##' \code{timebegin[i,2]} is \eqn{t_2}, \cr
##' \code{timebegin[i,3]} is \eqn{t_3}. \cr
##' @param timeend  an integer matrix of upper bounds for the event-history intervals. \cr
##' For example, if the event history for the ith subject is \cr
##' \eqn{[t_1,t_2)} in state \eqn{s_1}, \cr
##' \eqn{[t_2,t_3)} in state \eqn{s_2}, \cr
##' \eqn{[t_3,t_4)} in state \eqn{s_3}, \cr
##' and censored after \eqn{t_4}, then \cr
##' \code{timeend[i,1]} is \eqn{t_2}, \cr
##' \code{timeend[i,2]} is \eqn{t_3}, \cr
##' \code{timeend[i,3]} is \eqn{t_4}. \cr
##' @param ntransition an integer vector of the number of transitions.
##' For example, if the event history for the ith subject is \cr
##' \eqn{[t_1,t_2)} in state \eqn{s_1}, \cr
##' \eqn{[t_2,t_3)} in state \eqn{s_2}, \cr
##' \eqn{[t_3,t_4)} in state \eqn{s_3}, \cr
##' and censored after \eqn{t_4}, then \cr
##' \code{ntransition[i]} is 3.
##' @param possible_transitions a logical matrix indicating the possible transitions. \cr
##' \code{nrow(possible_transitions)} and \code{ncol(possible_transitions)} should be equal to the total number of states. \cr
##' \code{possible_transitions[s1,s2]} is \code{TRUE} if the transition from state \code{s1} to \code{s2} is possible. \cr
##' \code{possible_transitions[s1,s2]} is \code{FALSE} if the transition from state \code{s1} to \code{s2} is not possible. \cr
##' @param basismatrix_long a list of basis matrices. \cr
##' \code{basismatrix_long[[k]]} is the basis matrix for the kth longitudinal measure. \cr
##' \code{basismatrix_long[[k]][t,l]} is the lth basis function of the kth longitudinal measure evaluated at time \code{t}.
##' @param nqmc an integer value. \cr 
##' Number of quasi Monte Carlo points for evaluating numeric integration.
##' @return a list containing the model fit result.
##' @export
joint_multistate_model<-function(
  yyikj,ttikj,nobsik,zzi,xxi,
  statebegin,stateend,timebegin,timeend,ntransition,
  possible_transitions,
  basismatrix_long,
  nqmc=2000){
  
  if(!all(dim(yyikj)==dim(ttikj)))stop("dim(yyikj) and dim(ttikj) do not match.")
  if(dim(yyikj)[1]!=dim(nobsik)[1])stop("dim(yyikj) and dim(nobsik) do not match.")
  if(dim(yyikj)[2]!=dim(nobsik)[2])stop("dim(yyikj) and dim(nobsik) do not match.")
  if(dim(yyikj)[1]!=dim(zzi)[1])stop("dim(yyikj) and dim(zzi) do not match.")
  if(dim(yyikj)[1]!=dim(xxi)[1])stop("dim(yyikj) and dim(xxi) do not match.")
  if(dim(yyikj)[2]!=length(basismatrix_long))stop(paste0("basismatrix_long should be a list of length ",dim(yyikj)[2],"."))
  if(dim(possible_transitions)[1]!=dim(possible_transitions)[2])stop("possible_transitions is not square.")
  
  nind<-dim(yyikj)[1]
  ndimy<-dim(yyikj)[2]
  ndimz<-dim(zzi)[2]
  ndimx<-dim(xxi)[2]
  nstates<-dim(ntransition)[1]
  ntimepoints<-dim(basismatrix_long[[1]])[1]
  nbasy<-c()
  for(kk in 1:ndimy){
    if(dim(basismatrix_long[[kk]])[1]!=ntimepoints)stop("elements in basismatrix_long have unequal numbers of rows.")
    basismatrix_long[[kk]]<-rbind(0,basismatrix_long[[kk]])
    nbasy[kk]<-ncol(basismatrix_long[[kk]])
  }
  
  ## check yyikj ttikj
  for(ii in 1:nind){
    for(kk in 1:ndimy){
      if(!all(is.finite(yyikj[ii,kk,1:nobsik[ii,kk]])))stop("yyikj contains NA.")
      if(!all(is.finite(ttikj[ii,kk,1:nobsik[ii,kk]])))stop("ttikj contains NA.")
      if(!all(ttikj[ii,kk,1:nobsik[ii,kk]]%%1==0))stop("invalid values in ttikj.")
      if(!all(ttikj[ii,kk,1:nobsik[ii,kk]]>=1))stop("invalid values in ttikj.")
      if(!all(ttikj[ii,kk,1:nobsik[ii,kk]]<=ntimepoints))stop("invalid values in ttikj.")
    }
  }
  if(!all(nobsik%%1==0))stop("nobsik are not non-negative integers.")
  if(!all(nobsik%%1>=0))stop("nobsik are not non-negative integers.")
  if(!all(is.finite(zzi)))stop("zzi contains NA.")
  if(!all(is.finite(xxi)))stop("xxi contains NA.")
  for(ii in 1:nind){
    if(!all(statebegin[ii,1:ntransition[ii]]%in%(0:nstates)))stop("invalid values in statebegin.")
    if(!all(stateend[ii,1:ntransition[ii]]%in%(0:nstates)))stop("invalid values in stateend.")
    if(!all(timebegin[ii,1:ntransition[ii]]%%1==0))stop("invalid values in timebegin.")
    if(!all(timebegin[ii,1:ntransition[ii]]>=0))stop("invalid values in timebegin.")
    if(!all(timebegin[ii,1:ntransition[ii]]<=ntimepoints))stop("invalid values in timebegin.")
    if(!all(timeend[ii,1:ntransition[ii]]%%1==0))stop("invalid values in timeend.")
    if(!all(timeend[ii,1:ntransition[ii]]>=0))stop("invalid values in timeend.")
    if(!all(timeend[ii,1:ntransition[ii]]<=ntimepoints))stop("invalid values in timeend.")
  }
  possible_transitions<-as.logical(possible_transitions)
  if(any(is.na(possible_transitions)))stop("possible_transitions should be logical.")
  
  alleeieval<-t(randtoolbox::sobol(nqmc,dim=sum(nbasy),normal=TRUE))
  
  result_est<-multi_covariate(
    yyikj=yyikj,ttikj=ttikj,nind=nind,ndimy=ndimy,nobsik=nobsik,
    statebegin=statebegin-1,stateend=stateend-1,timebegin=timebegin,timeend=timeend,ntransition=ntransition,
    possible_transitions=possible_transitions+0,nstates=nstates,
    zzi=zzi,ndimz=ndimz,xxi=xxi,ndimx=ndimx,
    basismatrix_long=basismatrix_long,nbas_long=nbasy,alleeieval=alleeieval,
    ntimepoints=ntimepoints,
    factor_min=0.5,factor_max=0.5,tol=0.001,tol2=0.001,niter=500,csv=FALSE,tag="0")
  
  model_fit<-list(
    result_est=result_est,
    data=list(
      yyikj=yyikj,ttikj=ttikj,nind=nind,ndimy=ndimy,nobsik=nobsik,
      statebegin=statebegin,stateend=stateend,timebegin=timebegin,timeend=timeend,ntransition=ntransition,
      possible_transitions=possible_transitions,nstates=nstates,
      zzi=zzi,ndimz=ndimz,xxi=xxi,ndimx=ndimx,
      basismatrix_long=basismatrix_long,nbasy=nbasy,
      ntimepoints=ntimepoints))
  
  return(model_fit)
}

