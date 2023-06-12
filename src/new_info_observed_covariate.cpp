#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]

List info_observed_covariate(
    arma::cube yyikj,   // longitudinal
    arma::ucube ttikj,  // longitudinal
    int nind,           // longitudinal
    int ndimy,          // longitudinal
    arma::umat nobsik,  // longitudinal
    
    arma::imat statebegin,  // survival
    arma::imat stateend,    // survival
    arma::imat timebegin,   // survival
    arma::imat timeend,     // survival
    arma::ivec ntransition, // survival
    arma::imat possible_transitions,
    int nstates,
    
    arma::mat zzi,      // survival covariate
    int ndimz,          // survival covariate
    
    arma::mat xxi,      // longitudinal covariate
    int ndimx,          // longitudinal covariate
    
    arma::field<arma::mat> basismatrix_long, // basis
    arma::uvec nbas_long,                    // basis
    arma::mat alleeieval,                    // basis
    
    int ntimepoints,
    List result_est){
  
  const int ndim_rand=arma::sum(nbas_long);
  const int ndim1yz=1+ndimy+ndimz;
  
  int ii,jj,ll,ii2,jj2,ll2,ii3,jj3,ll3,tt,tt1,tt2,elem,ttdiff,iter,iter2,iitt,idx,idx2,kk;
  int idx_begin,idx_end,idx_begin2,idx_end2,idx_begin3,idx_end3;
  int idx1,idx_begin1,idx_end1;
  double temp_double,temp_double2;
  arma::vec temp_vec,temp_vec1,temp_vec2;
  arma::mat temp_mat;
  
  double omega=1.0/double(ntimepoints);
  double ttdiff_double;
  
  const int nabscissas=alleeieval.n_cols;
  
  int nmaxobs=yyikj.n_cols;
  int nmaxtransitions=stateend.n_cols;
  
  arma::uvec I0(ndimy),I1(ndimy);
  for(kk=0;kk<ndimy;kk++){
    I0(kk)=arma::sum(nbas_long.head(kk));
    I1(kk)=arma::sum(nbas_long.head(kk+1))-1;
  }
  
  arma::uvec possible_begin(nstates,arma::fill::zeros);
  int npossible_transitions=0;
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    for(idx_end=0;idx_end<nstates;idx_end++){
      if(possible_transitions(idx_begin,idx_end)==1){
        possible_begin(idx_begin)=1;
        npossible_transitions=npossible_transitions+1;
      }
    }
  }
  arma::vec possible_transitions_begin(npossible_transitions);
  arma::vec possible_transitions_end(npossible_transitions);
  idx=0;
  arma::umat beginend2idx(nstates,nstates);
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    for(idx_end=0;idx_end<nstates;idx_end++){
      if(possible_transitions(idx_begin,idx_end)==1){
        possible_transitions_begin(idx)=idx_begin;
        possible_transitions_end(idx)=idx_end;
        beginend2idx(idx_begin,idx_end)=idx;
        idx=idx+1;
      }
    }
  }
  arma::umat F_ref(nstates,nstates,arma::fill::zeros);
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    for(idx_end=0;idx_end<nstates;idx_end++){
      if(possible_transitions(idx_begin,idx_end)==1){
        F_ref(idx_begin,idx_end)=1;
        break;
      }
    }
  }
  
  arma::vec mub=result_est(0);
  arma::mat Sigmab=result_est(1);Sigmab=0.5*(Sigmab+Sigmab.t());
  arma::mat Sigmab_inv=arma::inv_sympd(Sigmab);Sigmab_inv=0.5*(Sigmab_inv+Sigmab_inv.t());
  arma::vec Sigmae2=result_est(2);
  arma::vec Sigmae2_inv=1.0/Sigmae2;
  arma::mat F_alpha=result_est(3);
  arma::field<arma::vec> F_betay_flat=result_est(4);
  arma::field<arma::vec> F_betaz_flat=result_est(5);
  // arma::field<arma::vec> F_betay=result_est(4);
  // arma::field<arma::vec> F_betaz=result_est(5);
  arma::field<arma::vec> F_hazard0=result_est(6);
  arma::field<arma::uvec> F_hazard0_count=result_est(7);
  arma::field<arma::vec> F_xi=result_est(8);
  
  arma::field<arma::vec> F_betay(nstates,nstates);
  arma::field<arma::vec> F_betaz(nstates,nstates);
  
  for(idx=0;idx<npossible_transitions;idx++){
    idx_begin=possible_transitions_begin(idx);
    idx_end=possible_transitions_end(idx);
    F_betay(idx_begin,idx_end)=F_betay_flat(idx);
    F_betaz(idx_begin,idx_end)=F_betaz_flat(idx);
  }
  
  arma::field<arma::vec> F_basismatrix_long(ndimy,ntimepoints+1);
  for(kk=0;kk<ndimy;kk++){
    for(tt=1;tt<ntimepoints+1;tt++){
      F_basismatrix_long(kk,tt)=basismatrix_long(kk).row(tt).t();
    }
  }
  
  arma::field<arma::umat> F_need_Eexp(nstates,nstates);
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    for(idx_end=0;idx_end<nstates;idx_end++){
      if(possible_transitions(idx_begin,idx_end)==1){
        F_need_Eexp(idx_begin,idx_end)=arma::zeros<arma::umat>(nind,ntimepoints+1);
      }
    }
  }
  
  // step 2
  
  arma::field<arma::vec> F_Ebi(nind);
  arma::field<arma::mat> F_Vbi(nind);
  arma::field<arma::mat> F_Vbi_half(nind);
  
  arma::mat allbbieval(ndim_rand,nabscissas);
  
  arma::vec F_scale(nind,arma::fill::zeros);
  arma::vec F_range(nind,arma::fill::zeros);
  
  arma::vec F_Inti0(nind,arma::fill::zeros);
  arma::field<arma::vec> F_Inti1(nind);
  arma::field<arma::mat> F_Inti2(nind);
  arma::field<arma::vec> F_Inti3(nind);
  arma::field<double> F_Inti4(nind,npossible_transitions,ntimepoints+1);
  arma::field<arma::vec> F_Inti5(nind,npossible_transitions,ntimepoints+1);
  arma::field<arma::mat> F_Inti6(nind,npossible_transitions,ntimepoints+1);
  
  arma::field<arma::vec> F_Ebbi(nind);
  arma::field<arma::mat> F_Ebbibbit(nind);
  arma::field<arma::vec> F_Eeee2(nind);
  arma::field<double> F_Eexp(nind,npossible_transitions,ntimepoints+1);
  arma::field<arma::vec> F_Evexp(nind,npossible_transitions,ntimepoints+1);
  arma::field<arma::mat> F_Evvexp(nind,npossible_transitions,ntimepoints+1);
  
  double Inti0_sum;
  arma::vec Inti1_sum(ndim_rand);
  arma::mat Inti2_sum(ndim_rand,ndim_rand);
  arma::vec Inti3_sum(ndimy);
  
  arma::field<double> Inti4_sum(npossible_transitions,ntimepoints+1);
  arma::field<arma::vec> Inti5_sum(npossible_transitions,ntimepoints+1);
  arma::field<arma::mat> Inti6_sum(npossible_transitions,ntimepoints+1);
  
  for(ii=0;ii<nind;ii++){
    for(ll=0;ll<ntransition(ii);ll++){
      idx_begin=statebegin(ii,ll);
      idx_end=stateend(ii,ll);
      if(idx_end==-1)continue;;
      tt=timeend(ii,ll);
      for(ii2=0;ii2<nind;ii2++){
        for(ll2=0;ll2<ntransition(ii2);ll2++){
          if(statebegin(ii2,ll2)!=idx_begin)continue;
          // timebegin < tt <= timeend
          if(timebegin(ii2,ll2)>=tt)continue;
          if(timeend(ii2,ll2)<tt)continue;
          for(idx2=0;idx2<npossible_transitions;idx2++){
            if(possible_transitions_begin(idx2)!=idx_begin)continue;
            idx_begin2=possible_transitions_begin(idx2);
            idx_end2=possible_transitions_end(idx2);
            F_need_Eexp(idx_begin2,idx_end2)(ii2,tt)=1;
          }
        }
      }
    }
  }
  
  arma::field<arma::vec> F_yik(nind,ndimy);
  arma::field<arma::uvec> F_tik(nind,ndimy);
  arma::field<arma::mat> F_Phi(nind,ndimy);
  arma::field<arma::vec> F_zzi(nind);
  arma::field<arma::mat> F_xxik(nind,ndimy);
  arma::mat Phi;
  // arma::field<arma::mat> F_Basis(nind);
  // arma::mat Basis;
  arma::uvec ttsubvec;
  arma::vec yysubvec;
  arma::mat xxsubmat;
  arma::mat PhitxSigmae;
  
  for(ii=0;ii<nind;ii++){
    F_zzi(ii)=zzi.row(ii).t();
    for(kk=0;kk<ndimy;kk++){
      ttsubvec.zeros(nobsik(ii,kk));
      yysubvec.zeros(nobsik(ii,kk));
      xxsubmat.zeros(nobsik(ii,kk),ndimx);
      for(jj=0;jj<nobsik(ii,kk);jj++){
        ttsubvec(jj)=ttikj(ii,kk,jj);
        yysubvec(jj)=yyikj(ii,kk,jj);
        xxsubmat.row(jj)=xxi.row(ii);
      }
      F_tik(ii,kk)=ttsubvec;
      F_yik(ii,kk)=yysubvec;
      F_Phi(ii,kk)=basismatrix_long(kk).rows(ttsubvec);
      F_xxik(ii,kk)=xxsubmat;
    }
  }
  
  ///--- Estep
  arma::vec yyvec;
  arma::mat Phi_temp;
  arma::vec bbitemp(ndim_rand),yytemp(ndimy),vvtemp(ndim1yz);
  double sumterm,tt_double,deltaterm,Tcondb;
  // arma::vec avec(ndim);
  arma::vec ytemp,muhat;
  arma::mat resid,resid2;
  arma::vec bbi_sum(ndim_rand);
  arma::mat bbibbit_sum(ndim_rand,ndim_rand);
  arma::vec eee2_sum(ndimy);
  arma::mat xxitxxi_sum(ndimx,ndimx);
  arma::vec xxityyi_sum(ndimx);
  
  double expterm,sum_risk;
  double logF_term,logh_term,loglik,exploglik;
  
  arma::vec all_loglik(nabscissas);
  arma::vec all_exploglik(nabscissas);
  double scale;
  double range;
  
  ////// start here //////
  
  arma::vec theta(ndim1yz);
  
  arma::vec F_expl(nind);
  arma::field<arma::vec> F_dtheta(nind,npossible_transitions);
  arma::field<arma::mat> F_ddtheta(nind,npossible_transitions);
  arma::field<arma::mat> F_dthetaxdtheta(nind,npossible_transitions,npossible_transitions);
  arma::field<double> F_dhazard(nind,nstates,ntimepoints+1);
  arma::field<double> F_ddhazard(nind,nstates,ntimepoints+1);
  arma::field<arma::mat> F_dhazardxdhazard(nind,nstates);
  arma::field<arma::vec> F_dthetaxdhazard(nind,npossible_transitions,ntimepoints);
  arma::field<arma::vec> F_dthetadhazard(nind,npossible_transitions,ntimepoints);
  
  double sum_expl;
  arma::vec sum_dtheta(ndim1yz);
  arma::mat sum_ddtheta(ndim1yz,ndim1yz);
  arma::mat sum_dthetaxdtheta(ndim1yz,ndim1yz);
  double sum_dhazard;
  double sum_ddhazard;
  arma::mat sum_dhazardxdhazard(ntimepoints+1,ntimepoints+1);
  arma::vec sum_dthetaxdhazard(ndim1yz);
  arma::vec sum_dthetadhazard(ndim1yz);
  
  arma::field<arma::vec> F_dlogf_dbeta(nabscissas,npossible_transitions);
  arma::field<arma::mat> F_ddlogf_ddbeta(nabscissas,npossible_transitions);
  arma::field<double> F_dlogf_dhazard(nabscissas,nstates,ntimepoints+1);
  arma::field<double> F_ddlogf_ddhazard(nabscissas,nstates,ntimepoints+1);
  arma::field<arma::vec> F_ddlogf_dbetadhazard(nabscissas,npossible_transitions,ntimepoints+1);
  
  arma::vec a_dlogf_dbeta(ndim1yz);
  arma::mat a_ddlogf_ddbeta(ndim1yz,ndim1yz);
  double a_dlogf_dhazard;
  double a_ddlogf_ddhazard;
  arma::vec a_ddlogf_dbetadhazard(ndim1yz);
  
  arma::field<arma::vec> F_theta(nstates,nstates);
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    for(idx_end=0;idx_end<nstates;idx_end++){
      if(possible_transitions(idx_begin,idx_end)){
        theta(0)=F_alpha(idx_begin,idx_end);
        theta.subvec(1,ndimy)=F_betay(idx_begin,idx_end);
        theta.subvec(ndimy+1,ndimy+ndimz)=F_betaz(idx_begin,idx_end);
        F_theta(idx_begin,idx_end)=theta;
      }
    }
  }
  
  for(ii=0;ii<nind;ii++){
    
    Rcpp::Rcout<<ii<<" ";
    temp_mat=Sigmab_inv;
    for(kk=0;kk<ndimy;kk++){
      temp_mat.submat(I0(kk),I0(kk),I1(kk),I1(kk))+=Sigmae2_inv(kk)*(F_Phi(ii,kk).t()*F_Phi(ii,kk));
    }
    
    temp_mat=0.5*(temp_mat+temp_mat.t());
    F_Vbi(ii)=arma::inv_sympd(temp_mat);
    F_Vbi(ii)=0.5*(F_Vbi(ii)+F_Vbi(ii).t());
    F_Vbi_half(ii)=arma::sqrtmat_sympd(F_Vbi(ii));
    F_Vbi_half(ii)=0.5*(F_Vbi_half(ii)+F_Vbi_half(ii).t());
    
    temp_vec=arma::zeros(ndim_rand);
    for(kk=0;kk<ndimy;kk++){
      temp_vec.subvec(I0(kk),I1(kk))=Sigmae2_inv(kk)*(F_Phi(ii,kk).t()*
        (F_yik(ii,kk)-F_Phi(ii,kk)*mub.subvec(I0(kk),I1(kk))-F_xxik(ii,kk)*F_xi(kk)));
    }
    F_Ebi(ii)=mub+F_Vbi(ii)*temp_vec;
    
    for(elem=0;elem<nabscissas;elem++){
      allbbieval.col(elem)=F_Vbi_half(ii)*alleeieval.col(elem)+F_Ebi(ii);
      bbitemp=allbbieval.col(elem);
      
      loglik=0.0;
      
      for(ll=0;ll<ntransition(ii);ll++){
        idx_begin=statebegin(ii,ll);
        idx_end=stateend(ii,ll);
        if(idx_end==-1)continue;
        if(possible_transitions(idx_begin,idx_end)!=1)continue;
        tt=timeend(ii,ll);
        for(kk=0;kk<ndimy;kk++){
          yytemp(kk)=
            arma::dot(F_basismatrix_long(kk,tt),bbitemp.subvec(I0(kk),I1(kk)))+
            arma::dot(xxi.row(ii),F_xi(kk));
        }
        vvtemp(0)=1.0;
        vvtemp.subvec(1,ndimy)=yytemp;
        vvtemp.subvec(ndimy+1,ndimy+ndimz)=F_zzi(ii);
        
        sumterm=arma::dot(F_theta(idx_begin,idx_end),vvtemp);
        logh_term=std::log(F_hazard0(idx_begin)(tt))+sumterm;
        loglik=loglik+logh_term;
      }
      
      for(ll=0;ll<ntransition(ii);ll++){
        idx_begin=statebegin(ii,ll);
        
        logF_term=0.0;
        for(ii2=0;ii2<nind;ii2++){
          for(ll2=0;ll2<ntransition(ii2);ll2++){
            tt=timeend(ii2,ll2);
            idx_begin2=statebegin(ii2,ll2);
            idx_end2=stateend(ii2,ll2);
            if(idx_begin!=idx_begin2)continue;
            if(idx_end2==-1)continue;
            // at risk if
            // timebegin < tt <= timeend
            if(timebegin(ii,ll)>=tt)continue;
            if(timeend(ii,ll)<tt)continue;
            yytemp.zeros(ndimy);
            for(kk=0;kk<ndimy;kk++){
              yytemp(kk)=
                arma::dot(F_basismatrix_long(kk,tt),bbitemp.subvec(I0(kk),I1(kk)))+
                arma::dot(xxi.row(ii),F_xi(kk));
            }
            vvtemp(0)=1.0;
            vvtemp.subvec(1,ndimy)=yytemp;
            vvtemp.subvec(ndimy+1,ndimy+ndimz)=F_zzi(ii);
            
            for(idx=0;idx<npossible_transitions;idx++){
              if(possible_transitions_begin(idx)!=idx_begin)continue;
              idx_end=possible_transitions_end(idx);
              expterm=std::exp(arma::dot(F_theta(idx_begin,idx_end),vvtemp));
              logF_term=logF_term-F_hazard0(idx_begin)(tt)*expterm;
            }
          }
        }
        loglik=loglik+logF_term;
      }
      
      all_loglik(elem)=loglik;
    }
    
    scale=all_loglik.max();
    range=all_loglik.max()-all_loglik.min();
    F_scale(ii)=scale;
    F_range(ii)=range;
    all_loglik=all_loglik-scale;
    all_exploglik=arma::exp(all_loglik);
    
    for(elem=0;elem<nabscissas;elem++){
      allbbieval.col(elem)=F_Vbi_half(ii)*alleeieval.col(elem)+F_Ebi(ii);
      bbitemp=allbbieval.col(elem);
      
      // F_dlogf_dbeta, F_ddlogf_ddbeta
      for(idx=0;idx<npossible_transitions;idx++){
        idx_begin=possible_transitions_begin(idx);
        idx_end=possible_transitions_end(idx);
        a_dlogf_dbeta.fill(0.0);
        a_ddlogf_ddbeta.fill(0.0);
        for(ll=0;ll<ntransition(ii);ll++){
          if(statebegin(ii,ll)!=idx_begin)continue;
          if(stateend(ii,ll)!=idx_end)continue;
          tt=timeend(ii,ll);
          for(kk=0;kk<ndimy;kk++){
            yytemp(kk)=
              arma::dot(F_basismatrix_long(kk,tt),bbitemp.subvec(I0(kk),I1(kk)))+
              arma::dot(xxi.row(ii),F_xi(kk));
          }
          vvtemp(0)=1.0;
          vvtemp.subvec(1,ndimy)=yytemp;
          vvtemp.subvec(ndimy+1,ndimy+ndimz)=F_zzi(ii);
          a_dlogf_dbeta=a_dlogf_dbeta+vvtemp;
        }
        for(ll=0;ll<ntransition(ii);ll++){
          if(statebegin(ii,ll)!=idx_begin)continue;
          for(tt=timebegin(ii,ll)+1;tt<=timeend(ii,ll);tt++){
            if(F_hazard0_count(idx_begin)(tt)==0)continue;
            for(kk=0;kk<ndimy;kk++){
              yytemp(kk)=
                arma::dot(F_basismatrix_long(kk,tt),bbitemp.subvec(I0(kk),I1(kk)))+
                arma::dot(xxi.row(ii),F_xi(kk));
            }
            vvtemp(0)=1.0;
            vvtemp.subvec(1,ndimy)=yytemp;
            vvtemp.subvec(ndimy+1,ndimy+ndimz)=F_zzi(ii);
            expterm=std::exp(arma::dot(F_theta(idx_begin,idx_end),vvtemp));
            a_dlogf_dbeta=a_dlogf_dbeta-
              F_hazard0_count(idx_begin)(tt)*F_hazard0(idx_begin)(tt)*expterm*vvtemp;
            a_ddlogf_ddbeta=a_ddlogf_ddbeta-
              F_hazard0_count(idx_begin)(tt)*F_hazard0(idx_begin)(tt)*expterm*vvtemp*vvtemp.t();
          }
        }
        F_dlogf_dbeta(elem,idx)=a_dlogf_dbeta;
        F_ddlogf_ddbeta(elem,idx)=a_ddlogf_ddbeta;
      }
      
      // a_dlogf_dhazard, a_ddlogf_ddhazard
      for(idx_begin=0;idx_begin<nstates;idx_begin++){
        if(!possible_begin(idx_begin))continue;
        for(tt=1;tt<ntimepoints+1;tt++){
          if(F_hazard0_count(idx_begin)(tt)==0)continue;
          for(kk=0;kk<ndimy;kk++){
            yytemp(kk)=
              arma::dot(F_basismatrix_long(kk,tt),bbitemp.subvec(I0(kk),I1(kk)))+
              arma::dot(xxi.row(ii),F_xi(kk));
          }
          vvtemp(0)=1.0;
          vvtemp.subvec(1,ndimy)=yytemp;
          vvtemp.subvec(ndimy+1,ndimy+ndimz)=F_zzi(ii);
          a_dlogf_dhazard=0.0;
          a_ddlogf_ddhazard=0.0;
          for(ll=0;ll<ntransition(ii);ll++){
            if(statebegin(ii,ll)!=idx_begin)continue;
            if(timeend(ii,ll)==-1)continue;
            if(timeend(ii,ll)==tt){
              a_dlogf_dhazard=a_dlogf_dhazard+1.0/F_hazard0(idx_begin)(tt);
              a_ddlogf_ddhazard=a_ddlogf_ddhazard-1.0/F_hazard0(idx_begin)(tt)/F_hazard0(idx_begin)(tt);
            }
          }
          for(ll=0;ll<ntransition(ii);ll++){
            if(statebegin(ii,ll)!=idx_begin)continue;
            if(timebegin(ii,ll)>=tt)continue;
            if(timeend(ii,ll)<tt)continue;
            for(idx_end=0;idx_end<nstates;idx_end++){
              if(!possible_transitions(idx_begin,idx_end))continue;
              expterm=std::exp(arma::dot(F_theta(idx_begin,idx_end),vvtemp));
              a_dlogf_dhazard=a_dlogf_dhazard-F_hazard0_count(idx_begin)(tt)*expterm;
            }
          }
          F_dlogf_dhazard(elem,idx_begin,tt)=a_dlogf_dhazard;
          F_ddlogf_ddhazard(elem,idx_begin,tt)=a_ddlogf_ddhazard;
        }
      }
      
      // a_ddlogf_dbetadhazard
      for(idx=0;idx<npossible_transitions;idx++){
        idx_begin=possible_transitions_begin(idx);
        idx_end=possible_transitions_end(idx);
        for(tt=1;tt<ntimepoints+1;tt++){
          if(F_hazard0_count(idx_begin)(tt)==0)continue;
          for(kk=0;kk<ndimy;kk++){
            yytemp(kk)=
              arma::dot(F_basismatrix_long(kk,tt),bbitemp.subvec(I0(kk),I1(kk)))+
              arma::dot(xxi.row(ii),F_xi(kk));
          }
          vvtemp(0)=1.0;
          vvtemp.subvec(1,ndimy)=yytemp;
          vvtemp.subvec(ndimy+1,ndimy+ndimz)=F_zzi(ii);
          a_ddlogf_dbetadhazard.fill(0.0);
          for(ll=0;ll<ntransition(ii);ll++){
            if(statebegin(ii,ll)!=idx_begin)continue;
            if(timebegin(ii,ll)>=tt)continue;
            if(timeend(ii,ll)<tt)continue;
            expterm=std::exp(arma::dot(F_theta(idx_begin,idx_end),vvtemp));
            a_ddlogf_dbetadhazard=a_ddlogf_dbetadhazard-F_hazard0_count(idx_begin)(tt)*expterm*vvtemp;
          }
          F_ddlogf_dbetadhazard(elem,idx,tt)=a_ddlogf_dbetadhazard;
        }
      }
    }
    
    sum_expl=0.0;
    for(elem=0;elem<nabscissas;elem++){
      sum_expl=sum_expl+all_exploglik(elem);
    }
    F_expl(ii)=sum_expl/double(nabscissas);
    
    for(idx=0;idx<npossible_transitions;idx++){
      idx_begin=possible_transitions_begin(idx);
      idx_end=possible_transitions_end(idx);
      sum_dtheta.fill(0.0);
      sum_ddtheta.fill(0.0);
      for(elem=0;elem<nabscissas;elem++){
        sum_dtheta=sum_dtheta+
          all_exploglik(elem)*F_dlogf_dbeta(elem,idx);
        sum_ddtheta=sum_ddtheta+
          all_exploglik(elem)*F_ddlogf_ddbeta(elem,idx);
      }
      F_dtheta(ii,idx)=sum_dtheta/double(nabscissas);
      F_ddtheta(ii,idx)=sum_ddtheta/double(nabscissas);
    }
    
    for(idx1=0;idx1<npossible_transitions;idx1++){
      idx_begin1=possible_transitions_begin(idx1);
      idx_end1=possible_transitions_end(idx1);
      for(idx2=0;idx2<npossible_transitions;idx2++){
        idx_begin2=possible_transitions_begin(idx2);
        idx_end2=possible_transitions_end(idx2);
        if(idx_begin1!=idx_begin2)continue;
        sum_dthetaxdtheta.fill(0.0);
        for(elem=0;elem<nabscissas;elem++){
          sum_dthetaxdtheta=sum_dthetaxdtheta+
            all_exploglik(elem)*F_dlogf_dbeta(elem,idx1)*F_dlogf_dbeta(elem,idx2).t();
        }
        F_dthetaxdtheta(ii,idx1,idx2)=sum_dthetaxdtheta/double(nabscissas);
      }
    }
    
    for(idx_begin=0;idx_begin<nstates;idx_begin++){
      if(!possible_begin(idx_begin))continue;
      for(tt=1;tt<ntimepoints+1;tt++){
        if(F_hazard0_count(idx_begin)(tt)==0)continue;
        sum_dhazard=0.0;
        sum_ddhazard=0.0;
        for(elem=0;elem<nabscissas;elem++){
          sum_dhazard=sum_dhazard+
            all_exploglik(elem)*F_dlogf_dhazard(elem,idx_begin,tt);
          sum_ddhazard=sum_ddhazard+
            all_exploglik(elem)*F_ddlogf_ddhazard(elem,idx_begin,tt);
        }
        F_dhazard(ii,idx_begin,tt)=sum_dhazard/double(nabscissas);
        F_ddhazard(ii,idx_begin,tt)=sum_ddhazard/double(nabscissas);
      }
    }
    
    for(idx_begin=0;idx_begin<nstates;idx_begin++){
      if(!possible_begin(idx_begin))continue;
      sum_dhazardxdhazard.fill(0.0);
      for(tt1=1;tt1<ntimepoints+1;tt1++){
        if(F_hazard0_count(idx_begin)(tt1)==0)continue;
        for(tt2=1;tt2<ntimepoints+1;tt2++){
          if(F_hazard0_count(idx_begin)(tt2)==0)continue;
          for(elem=0;elem<nabscissas;elem++){
            sum_dhazardxdhazard(tt1,tt2)=sum_dhazardxdhazard(tt1,tt2)+
              all_exploglik(elem)*F_dlogf_dhazard(elem,idx_begin,tt1)*F_dlogf_dhazard(elem,idx_begin,tt2);
          }
        }
      }
      F_dhazardxdhazard(ii,idx_begin)=sum_dhazardxdhazard/double(nabscissas);
    }
    
    for(idx=0;idx<npossible_transitions;idx++){
      idx_begin=possible_transitions_begin(idx);
      idx_end=possible_transitions_end(idx);
      for(tt=1;tt<ntimepoints+1;tt++){
        if(F_hazard0_count(idx_begin)(tt)==0)continue;
        sum_dthetadhazard.fill(0.0);
        sum_dthetaxdhazard.fill(0.0);
        for(elem=0;elem<nabscissas;elem++){
          sum_dthetadhazard=sum_dthetadhazard+
            all_exploglik(elem)*F_ddlogf_dbetadhazard(elem,idx,tt);
          sum_dthetaxdhazard=sum_dthetaxdhazard+
            all_exploglik(elem)*F_dlogf_dhazard(elem,idx_begin,tt)*F_dlogf_dbeta(elem,idx);
        }
        F_dthetadhazard(ii,idx,tt)=sum_dthetadhazard/double(nabscissas);
        F_dthetaxdhazard(ii,idx,tt)=sum_dthetaxdhazard/double(nabscissas);
      }
    }
  }
  
  arma::field<arma::mat> F_Hesian_theta_theta(npossible_transitions,npossible_transitions);
  arma::field<arma::mat> F_Hesian_theta_hazard(npossible_transitions);
  arma::field<arma::mat> F_Hesian_hazard_hazard(nstates);
  
  arma::mat a_Hessian_theta_theta(ndim1yz,ndim1yz);
  arma::mat a_Hessian_theta_hazard(ndim1yz,ntimepoints+1);
  arma::mat a_Hessian_hazard_hazard(ntimepoints+1,ntimepoints+1);
  
  for(idx1=0;idx1<npossible_transitions;idx1++){
    idx_begin1=possible_transitions_begin(idx1);
    idx_end1=possible_transitions_end(idx1);
    for(idx2=0;idx2<npossible_transitions;idx2++){
      idx_begin2=possible_transitions_begin(idx2);
      idx_end2=possible_transitions_end(idx2);
      if(idx_begin1!=idx_begin2)continue;
      
      a_Hessian_theta_theta.fill(0.0);
      if(idx_end1==idx_end2){
        for(ii=0;ii<nind;ii++){
          a_Hessian_theta_theta=a_Hessian_theta_theta+
            F_dthetaxdtheta(ii,idx1,idx2)/F_expl(ii)+
            F_ddtheta(ii,idx1)/F_expl(ii)-
            F_dtheta(ii,idx1)*F_dtheta(ii,idx2).t()/F_expl(ii)/F_expl(ii);
        }
      }else{
        for(ii=0;ii<nind;ii++){
          a_Hessian_theta_theta=a_Hessian_theta_theta+
            F_dthetaxdtheta(ii,idx1,idx2)/F_expl(ii)-
            F_dtheta(ii,idx1)*F_dtheta(ii,idx2).t()/F_expl(ii)/F_expl(ii);
        }
      }
      F_Hesian_theta_theta(idx1,idx2)=a_Hessian_theta_theta;
    }
  }
  
  for(idx=0;idx<npossible_transitions;idx++){
    idx_begin=possible_transitions_begin(idx);
    idx_end=possible_transitions_end(idx);
    a_Hessian_theta_hazard.fill(0.0);
    for(tt=1;tt<ntimepoints+1;tt++){
      if(F_hazard0_count(idx_begin)(tt)==0)continue;
      for(ii=0;ii<nind;ii++){
        a_Hessian_theta_hazard.col(tt)=a_Hessian_theta_hazard.col(tt)+
          F_dthetaxdhazard(ii,idx,tt)/F_expl(ii)+
          F_dthetadhazard(ii,idx,tt)/F_expl(ii)-
          F_dhazard(ii,idx_begin,tt)/F_expl(ii)*F_dtheta(ii,idx)/F_expl(ii);
      }
    }
    F_Hesian_theta_hazard(idx)=a_Hessian_theta_hazard;
  }
  
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    if(!possible_begin(idx_begin))continue;
    a_Hessian_hazard_hazard.fill(0.0);
    for(tt1=1;tt1<ntimepoints+1;tt1++){
      if(F_hazard0_count(idx_begin)(tt1)==0)continue;
      for(tt2=1;tt2<ntimepoints+1;tt2++){
        if(F_hazard0_count(idx_begin)(tt2)==0)continue;
        for(ii=0;ii<nind;ii++){
          a_Hessian_hazard_hazard(tt1,tt2)=a_Hessian_hazard_hazard(tt1,tt2)+
            F_dhazardxdhazard(ii,idx_begin)(tt1,tt2)/F_expl(ii)-
            F_dhazard(ii,idx_begin,tt1)*F_dhazard(ii,idx_begin,tt2)/F_expl(ii)/F_expl(ii);
        }
        if(tt1==tt2){
          for(ii=0;ii<nind;ii++){
            a_Hessian_hazard_hazard(tt1,tt2)=a_Hessian_hazard_hazard(tt1,tt2)+
              F_ddhazard(ii,idx_begin,tt1)/F_expl(ii);
          }
        }
      }
    }
    F_Hesian_hazard_hazard(idx_begin)=a_Hessian_hazard_hazard;
  }
  
  /// --- result --- ///
  
  List result(10);
  
  result(0)=F_Hesian_theta_theta;
  result(1)=F_Hesian_theta_hazard;
  result(2)=F_Hesian_hazard_hazard;
  result(3)=F_hazard0_count;
  result(4)=possible_transitions_begin;
  result(5)=possible_transitions_end;
  
  return(result);
}

