##' Calculate the information matrix (based on the profile/partial likelihood)
##'
##' @title Calculate the information matrix (based on the profile/partial likelihood)
##' @param model_fit a model fit returned by \code{joint_multistate_model},
##' @param nqmc an integer value. \cr 
##' Number of quasi Monte Carlo points for evaluating numeric integration.
##' @return a list containing the information matrices
info_partial<-function(model_fit,nqmc=2000){
  
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
  
  result_info_partial<-info_partial_covariate(
    yyikj=yyikj,ttikj=ttikj,nind=nind,ndimy=ndimy,nobsik=nobsik,
    statebegin=statebegin-1,stateend=stateend-1,timebegin=timebegin,timeend=timeend,
    ntransition=ntransition,possible_transitions=possible_transitions+0,nstates=nstates,
    zzi=zzi,ndimz=ndimz,xxi=xxi,ndimx=ndimx,
    basismatrix_long=basismatrix_long,nbas_long=nbasy,alleeieval=alleeieval,
    ntimepoints=ntimepoints,result_est)
  
  possible_transitions_begin<-c(result_info_partial[[2]])+1
  possible_transitions_end<-c(result_info_partial[[3]])+1
  II_beta_beta_partial_expand<-field_to_matrix(result_info_partial[[1]])
  
  list_info_partial<-vector("list",nstates)
  
  for(idx_begin in 1:nstates){
    if(any(possible_transitions[idx_begin,])){
      all_idx_transitions<-which(possible_transitions_begin==idx_begin)
      all_idx_end<-possible_transitions_end[all_idx_transitions]
      
      idx_row<-II_beta_beta_partial_expand$idx_row%in%all_idx_transitions
      idx_col<-II_beta_beta_partial_expand$idx_col%in%all_idx_transitions
      temp_mat<-II_beta_beta_partial_expand$mat[idx_row,idx_col]
      info_mat0<--temp_mat[-1,-1]
      info_mat<-rbind(0,cbind(0,info_mat0))
      list_info_partial[[idx_begin]]<-info_mat
    }
  }
  return(list_info_partial)
}


