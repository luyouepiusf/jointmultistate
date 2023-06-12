#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]

List lik_covariate(
    arma::cube yyikj,   // longitudinal
    arma::ucube ttikj,  // longitudinal
    int nind,           // longitudinal
    int ndimy,      // longitudinal
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
  
  int ii,jj,ll,ii2,jj2,ll2,tt,elem,iter,idx,idx2,kk;
  int idx_begin,idx_end,idx_begin2,idx_end2;
  arma::vec temp_vec,temp_vec1,temp_vec2;
  arma::mat temp_mat;
  
  const int nabscissas=alleeieval.n_cols;
  
  int nmaxobs=yyikj.n_cols;
  int nmaxtransitions=stateend.n_cols;
  
  arma::uvec I0(ndimy),I1(ndimy);
  for(kk=0;kk<ndimy;kk++){
    I0(kk)=arma::sum(nbas_long.head(kk));
    I1(kk)=arma::sum(nbas_long.head(kk+1))-1;
  }
  
  int npossible_transitions=0;
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    for(idx_end=0;idx_end<nstates;idx_end++){
      if(possible_transitions(idx_begin,idx_end)==1){
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
  
  Rcpp::Rcout<<"Step 1 (a2)"<<std::endl;
  
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
  
  Rcpp::Rcout<<"Step 1 (d)"<<std::endl;
  
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
  
  Rcpp::Rcout<<"Step 1 (g)"<<std::endl;
  
  ///--- Estep
  arma::vec yyvec;
  arma::mat Phi_temp;
  arma::vec bbitemp(ndim_rand),yytemp(ndimy),vvtemp(ndim1yz);
  double sumterm,deltaterm;
  // arma::vec avec(ndim);
  arma::vec ytemp,muhat;
  arma::mat resid,resid2;
  arma::vec bbi_sum(ndim_rand);
  arma::mat bbibbit_sum(ndim_rand,ndim_rand);
  arma::vec eee2_sum(ndimy);
  arma::mat xxitxxi_sum(ndimx,ndimx);
  arma::vec xxityyi_sum(ndimx);
  
  double expterm;
  double logF_term,logh_term,loglik,exploglik;
  
  arma::vec all_loglik(nabscissas);
  double scale;
  double range;
  
  ////// start here //////
  
  arma::vec theta(ndim1yz);
  arma::field<arma::vec> F_theta(nstates,nstates);
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    for(idx_end=0;idx_end<nstates;idx_end++){
      if(possible_transitions(idx_begin,idx_end)){
        theta(0)=1.0;
        theta.subvec(1,ndimy)=F_betay(idx_begin,idx_end);
        theta.subvec(ndimy+1,ndimy+ndimz)=F_betaz(idx_begin,idx_end);
        F_theta(idx_begin,idx_end)=theta;
      }
    }
  }
  
  for(ii=0;ii<nind;ii++){

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
  }
  
  /// --- loglikelihood --- ///
  
  arma::vec alllogLLi(nabscissas),allLLi(nabscissas);
  arma::vec allloglikelihood(nind);
  double sumloglikelihood,maxlogLLi,LLi,logy_term,logb_term,logmc_term;
  double logdet_Sigmab,logdet_signb,logdet_Vbi,logdet_signV;
  arma::log_det(logdet_Sigmab,logdet_signb,2.0*3.14159265359*Sigmab);
  
  sumloglikelihood=0.0;
  for(ii=0;ii<nind;ii++){
    for(elem=0;elem<nabscissas;elem++){
      allbbieval.col(elem)=F_Vbi_half(ii)*alleeieval.col(elem)+F_Ebi(ii);
      bbitemp=allbbieval.col(elem);
      
      loglik=0.0;
      logy_term=0.0;
      for(kk=0;kk<ndimy;kk++){
        muhat=F_Phi(ii,kk)*bbitemp.subvec(I0(kk),I1(kk))+F_xxik(ii,kk)*F_xi(kk);
        resid=F_yik(ii,kk)-muhat;
        for(jj=0;jj<nobsik(ii,kk);jj++){
          logy_term=logy_term-0.5*resid(jj)*resid(jj)/Sigmae2(kk)-0.5*std::log(2.0*arma::datum::pi*Sigmae2(kk));
        }
      }
      loglik=loglik+logy_term;
      
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
      logb_term=-0.5*arma::dot(bbitemp-mub,Sigmab_inv*(bbitemp-mub))-0.5*logdet_Sigmab;
      loglik=loglik+logb_term;
      arma::log_det(logdet_Vbi,logdet_signV,2.0*3.14159265359*F_Vbi(ii));
      logmc_term=-0.5*arma::dot(bbitemp-F_Ebi(ii),arma::solve(F_Vbi(ii),bbitemp-F_Ebi(ii)))-0.5*logdet_Vbi;
      alllogLLi(elem)=loglik-logmc_term;
    }
    maxlogLLi=arma::max(alllogLLi);
    allLLi=arma::exp(alllogLLi-maxlogLLi);
    LLi=arma::sum(allLLi)/double(nabscissas);
    allloglikelihood(ii)=std::log(LLi)+maxlogLLi;
    sumloglikelihood=sumloglikelihood+allloglikelihood(ii);
  }
  
  /// --- result --- ///
  
  List result(20);
  result(0)=sumloglikelihood;
  result(1)=allloglikelihood;
  result(2)=F_Vbi;
  result(3)=F_Ebi;
  // result(2)=F_II_partial;
  
  return(result);
}

