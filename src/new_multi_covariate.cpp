#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List multi_covariate(
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
    double factor_min=0.2,
    double factor_max=0.8,
    double tol=1e-5,
    double tol2=0.01,
    int niter=500,
    bool csv=false,
    std::string tag="0"){
  
  arma::wall_clock timer;
  double timing;
  timer.tic();
  
  const int ndim_rand=arma::sum(nbas_long);
  
  int ii,jj,ll,ii2,jj2,ll2,tt,elem,ttdiff,iter,iter2,iitt,idx,idx2,kk;
  int idx_begin,idx_end,idx_begin2,idx_end2;
  double temp_double,temp_double2;
  arma::vec temp_vec,temp_vec1,temp_vec2;
  arma::mat temp_mat;
  
  arma::uvec I0(ndimy),I1(ndimy);
  for(kk=0;kk<ndimy;kk++){
    I0(kk)=arma::sum(nbas_long.head(kk));
    I1(kk)=arma::sum(nbas_long.head(kk+1))-1;
  }
  
  // int n_tt_eval=tt_eval.n_elem;
  double omega=1.0/double(ntimepoints);
  double ttdiff_double;
  
  const int nabscissas=alleeieval.n_cols;
  
  double epsilon_abs=100.0,epsilon_rel=100.0;
  double epsilon_abs_coef=100.0,epsilon_rel_coef=100.0;
  double factor=factor_min;
  
  ///--- Convert yyij to a field
  timer.tic();
  int nmaxobs=yyikj.n_cols;
  int nmaxtransitions=stateend.n_cols;
  
  arma::vec mub(ndim_rand,arma::fill::zeros);
  arma::mat Sigmab(ndim_rand,ndim_rand,arma::fill::eye);Sigmab=100.0*Sigmab;
  arma::mat Sigmab_inv=arma::inv_sympd(Sigmab);
  arma::vec Sigmae2(ndimy,arma::fill::ones);
  arma::vec Sigmae2_inv=1.0/Sigmae2;
  arma::field<arma::vec> F_xi(ndimy);
  for(kk=0;kk<ndimy;kk++)F_xi(kk)=arma::ones(ndimx);
  
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
  arma::field<arma::vec> F_betay(nstates,nstates);
  arma::field<arma::vec> F_betaz(nstates,nstates);
  arma::mat F_alpha(nstates,nstates);
  arma::vec betay_new(ndimy,arma::fill::zeros);
  arma::vec betaz_new(ndimz,arma::fill::zeros);
  arma::vec betay_old(ndimy,arma::fill::zeros);
  arma::vec betaz_old(ndimz,arma::fill::zeros);
  double alpha_new;
  double alpha_old;
  
  // relative change
  arma::vec mub_old(ndim_rand,arma::fill::zeros);
  arma::mat Sigmab_old(ndim_rand,ndim_rand,arma::fill::eye);
  arma::vec Sigmae2_old(ndimy,arma::fill::ones);
  
  arma::vec rel_mub(ndim_rand);
  arma::mat rel_Sigmab(ndim_rand,ndim_rand);
  arma::vec rel_Sigmae2(ndimy);
  arma::mat F_rel_alpha(nstates,nstates);
  arma::field<arma::vec> F_rel_betay(nstates,nstates);
  arma::field<arma::vec> F_rel_betaz(nstates,nstates);
  arma::field<arma::vec> F_rel_hazard0(nstates);
  arma::field<arma::vec> F_rel_xi(ndimy);
  
  arma::vec abs_mub(ndim_rand);
  arma::mat abs_Sigmab(ndim_rand,ndim_rand);
  arma::vec abs_Sigmae2(ndimy);
  arma::mat F_abs_alpha(nstates,nstates);
  arma::field<arma::vec> F_abs_betay(nstates,nstates);
  arma::field<arma::vec> F_abs_betaz(nstates,nstates);
  arma::field<arma::vec> F_abs_hazard0(nstates);
  arma::field<arma::vec> F_abs_xi(ndimy);
  
  arma::vec xi;
  arma::vec xi_old;
  
  idx=0;
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    for(idx_end=0;idx_end<nstates;idx_end++){
      if(possible_transitions(idx_begin,idx_end)==1){
        F_betay(idx_begin,idx_end)=betay_new;
        F_betaz(idx_begin,idx_end)=betaz_new;
        possible_transitions_begin(idx)=idx_begin;
        possible_transitions_end(idx)=idx_end;
        idx=idx+1;
      }
    }
  }
  
  arma::field<arma::vec> F_hazard0(nstates);
  arma::field<arma::uvec> F_hazard0_count(nstates);
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    for(idx_end=0;idx_end<nstates;idx_end++){
      if(possible_transitions(idx_begin,idx_end)==1){
        F_hazard0_count(idx_begin)=arma::zeros<arma::uvec>(ntimepoints+1);
        temp_vec=arma::zeros<arma::vec>(ntimepoints+1);
        temp_vec.fill(arma::datum::nan);
        F_hazard0(idx_begin)=temp_vec;
      }
    }
  }
  
  for(ii=0;ii<nind;ii++){
    for(ll=0;ll<ntransition(ii);ll++){
      idx_begin=statebegin(ii,ll);
      idx_end=stateend(ii,ll);
      if(idx_end==-1)continue;
      tt=timeend(ii,ll);
      F_hazard0(idx_begin)(tt)=0.01;
      F_hazard0_count(idx_begin)(tt)=F_hazard0_count(idx_begin)(tt)+1;
    }
  }
  
  F_alpha.fill(arma::datum::nan);
  for(idx_begin=0;idx_begin<nstates;idx_begin++){
    for(idx_end=0;idx_end<nstates;idx_end++){
      if(possible_transitions(idx_begin,idx_end)==1){
        F_alpha(idx_begin,idx_end)=0.0;
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
  
  arma::field<arma::vec> F_Ebi(nind);
  arma::field<arma::mat> F_Vbi(nind);
  arma::field<arma::mat> F_Vbi_half(nind);
  arma::field<arma::vec> F_muhat(nind,kk);
  arma::mat test_mat;
  arma::vec test_vec;
  
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
  arma::field<arma::vec> F_Eyexp(nind,npossible_transitions,ntimepoints+1);
  arma::field<arma::mat> F_Eyyexp(nind,npossible_transitions,ntimepoints+1);
  
  double Inti0_sum;
  arma::vec Inti1_sum(ndim_rand);
  arma::mat Inti2_sum(ndim_rand,ndim_rand);
  arma::vec Inti3_sum(ndimy);
  
  arma::field<double> Inti4_sum(npossible_transitions,ntimepoints+1);
  arma::field<arma::vec> Inti5_sum(npossible_transitions,ntimepoints+1);
  arma::field<arma::mat> Inti6_sum(npossible_transitions,ntimepoints+1);
  
  arma::field<arma::vec> F_basismatrix_long(ndimy,ntimepoints+1);
  // arma::field<arma::vec> F_basismatrix_coef(ntimepoints);
  
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
  
  List result_temp(10);
  
  ///--- Estep
  arma::vec yyvec;
  arma::mat Phi_temp;
  arma::vec bbitemp(ndim_rand),yytemp(ndimy);
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
  double scale;
  double range;
  
  double denominatory;
  arma::vec numerator1y(ndimy);
  arma::mat numerator2y(ndimy,ndimy);
  
  double denominatorz;
  arma::vec numerator1z(ndimz);
  arma::mat numerator2z(ndimz,ndimz);
  
  double denominatora;
  double numerator1a;
  double numerator2a;
  
  
  arma::vec SSy_new(ndimy);
  arma::mat IIy_new(ndimy,ndimy);
  arma::vec SSz_new(ndimz);
  arma::mat IIz_new(ndimz,ndimz);
  double SSa_new;
  double IIa_new;
  
  arma::field<arma::vec> F_SSy(nstates,nstates);
  arma::field<arma::mat> F_IIy(nstates,nstates);
  
  arma::field<arma::vec> F_SSz(nstates,nstates);
  arma::field<arma::mat> F_IIz(nstates,nstates);
  
  arma::mat F_SSa(nstates,nstates,arma::fill::zeros);
  arma::mat F_IIa(nstates,nstates,arma::fill::zeros);
  
  // output csv
  
  
  std::ofstream mub_csv;
  std::ofstream Sigmab_csv;
  std::ofstream Sigmae2_csv;
  std::ofstream F_xi_csv;
  std::ofstream F_betay_csv;
  std::ofstream F_betaz_csv;
  std::ofstream F_alpha_csv;
  
  
  if(csv){
    mub_csv.open("mub"+tag+".csv");
    Sigmab_csv.open("Sigmab"+tag+".csv");
    Sigmae2_csv.open("Sigmae2"+tag+".csv");
    F_xi_csv.open("F_xi"+tag+".csv");
    F_betay_csv.open("F_betay"+tag+".csv");
    F_betaz_csv.open("F_betaz"+tag+".csv");
    F_alpha_csv.open("F_alpha"+tag+".csv");
  }
  
  for(iter=0;iter<niter;iter++){
    factor=double(iter)/double(niter)*factor_max+double(niter-iter)/double(niter)*factor_min;
    Rcpp::Rcout<<"iter: "<<iter<<", factor: "<<factor<<std::endl;
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
      
      Inti0_sum=0.0;
      Inti1_sum.fill(0.0);
      Inti2_sum.fill(0.0);
      Inti3_sum.fill(0.0);
      for(idx=0;idx<npossible_transitions;idx++){
        idx_begin=possible_transitions_begin(idx);
        idx_end=possible_transitions_end(idx);
        for(tt=1;tt<ntimepoints+1;tt++){
          if(F_need_Eexp(idx_begin,idx_end)(ii,tt)!=1)continue;
          Inti4_sum(idx,tt)=0.0;
          Inti5_sum(idx,tt)=arma::zeros(ndimy);
          Inti6_sum(idx,tt)=arma::zeros(ndimy,ndimy);
        }
      }
      
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
          sumterm=
            arma::dot(F_betay(idx_begin,idx_end),yytemp)+
            arma::dot(F_betaz(idx_begin,idx_end),F_zzi(ii));
          logh_term=std::log(F_hazard0(idx_begin)(tt))+F_alpha(idx_begin,idx_end)+sumterm;
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
              for(idx=0;idx<npossible_transitions;idx++){
                if(possible_transitions_begin(idx)!=idx_begin)continue;
                idx_end=possible_transitions_end(idx);
                expterm=std::exp(
                  arma::dot(F_betay(idx_begin,idx_end),yytemp)+
                    arma::dot(F_betaz(idx_begin,idx_end),F_zzi(ii)));
                logF_term=logF_term-F_hazard0(idx_begin)(tt)*std::exp(F_alpha(idx_begin,idx_end))*expterm;
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
      
      // if(all_loglik.max()-scale>10.0){
      //   Rcpp::Rcout<<"{"<<ii<<"}"<<"["<<scale<<"]"<<all_loglik.t();
      // }
      
      temp_vec=arma::zeros<arma::vec>(nabscissas);
      temp_vec.fill(scale-10.0);
      all_loglik=arma::max(all_loglik,temp_vec);
      temp_vec=arma::zeros<arma::vec>(nabscissas);
      temp_vec.fill(scale+10.0);
      all_loglik=arma::min(all_loglik,temp_vec);
      
      for(elem=0;elem<nabscissas;elem++){
        
        allbbieval.col(elem)=F_Vbi_half(ii)*alleeieval.col(elem)+F_Ebi(ii);//////
        bbitemp=allbbieval.col(elem);
        loglik=all_loglik(elem);
        
        exploglik=std::exp(loglik-scale);
        
        Inti0_sum=Inti0_sum+exploglik;
        Inti1_sum=Inti1_sum+exploglik*bbitemp;
        Inti2_sum=Inti2_sum+exploglik*(bbitemp-mub)*(bbitemp-mub).t();
        
        // yyvec=yyij.submat(ii,0,ii,nobs(ii)-1).t();
        // ytemp=yyvecij.submat(ii,0,ii,nobs(ii)*ndim-1).t();
        eee2_sum.fill(0.0);
        for(kk=0;kk<ndimy;kk++){
          muhat=F_Phi(ii,kk)*bbitemp.subvec(I0(kk),I1(kk))+F_xxik(ii,kk)*F_xi(kk);
          resid=F_yik(ii,kk)-muhat;
          for(jj=0;jj<nobsik(ii,kk);jj++){
            eee2_sum(kk)=eee2_sum(kk)+resid(jj)*resid(jj);
          }
        }
        Inti3_sum=Inti3_sum+exploglik*eee2_sum;
        
        for(idx=0;idx<npossible_transitions;idx++){
          idx_begin=possible_transitions_begin(idx);
          idx_end=possible_transitions_end(idx);
          for(tt=1;tt<ntimepoints+1;tt++){
            if(F_need_Eexp(idx_begin,idx_end)(ii,tt)!=1)continue;
            yytemp.zeros(ndimy);
            for(kk=0;kk<ndimy;kk++){
              yytemp(kk)=
                arma::dot(F_basismatrix_long(kk,tt),bbitemp.subvec(I0(kk),I1(kk)))+
                arma::dot(xxi.row(ii),F_xi(kk));
            }
            expterm=std::exp(
              arma::dot(F_betay(idx_begin,idx_end),yytemp)+
                arma::dot(F_betaz(idx_begin,idx_end),F_zzi(ii)));
            Inti4_sum(idx,tt)=Inti4_sum(idx,tt)+exploglik*expterm;
            Inti5_sum(idx,tt)=Inti5_sum(idx,tt)+exploglik*expterm*yytemp;
            Inti6_sum(idx,tt)=Inti6_sum(idx,tt)+exploglik*expterm*yytemp*yytemp.t();
          }
        }
      }
      
      F_Inti0(ii)=Inti0_sum/double(nabscissas);
      F_Inti1(ii)=Inti1_sum/double(nabscissas);
      F_Inti2(ii)=Inti2_sum/double(nabscissas);
      F_Inti3(ii)=Inti3_sum/double(nabscissas);
      
      // Rcpp::Rcout<<"{"<<ii<<"}***"<<F_Inti0(ii)<<std::endl;
      
      F_Ebbi(ii)=F_Inti1(ii)/F_Inti0(ii);
      F_Ebbibbit(ii)=F_Inti2(ii)/F_Inti0(ii);
      F_Eeee2(ii)=F_Inti3(ii)/F_Inti0(ii);
      
      for(idx=0;idx<npossible_transitions;idx++){
        idx_begin=possible_transitions_begin(idx);
        idx_end=possible_transitions_end(idx);
        for(tt=1;tt<ntimepoints+1;tt++){
          if(F_need_Eexp(idx_begin,idx_end)(ii,tt)!=1)continue;
          F_Inti4(ii,idx,tt)=Inti4_sum(idx,tt)/double(nabscissas);
          F_Inti5(ii,idx,tt)=Inti5_sum(idx,tt)/double(nabscissas);
          F_Inti6(ii,idx,tt)=Inti6_sum(idx,tt)/double(nabscissas);
          
          F_Eexp(ii,idx,tt)=F_Inti4(ii,idx,tt)/F_Inti0(ii);
          F_Eyexp(ii,idx,tt)=F_Inti5(ii,idx,tt)/F_Inti0(ii);
          F_Eyyexp(ii,idx,tt)=F_Inti6(ii,idx,tt)/F_Inti0(ii);
        }
      }
    }
    
    bbi_sum.fill(0.0);
    bbibbit_sum.fill(0.0);
    eee2_sum.fill(0.0);
    for(ii=0;ii<nind;ii++){
      bbi_sum=bbi_sum+F_Ebbi(ii);
      bbibbit_sum=bbibbit_sum+F_Ebbibbit(ii);
      eee2_sum=eee2_sum+F_Eeee2(ii);
    }
    
    mub_old=mub;
    mub=bbi_sum/double(nind);
    Sigmab_old=Sigmab;
    Sigmab=bbibbit_sum/double(nind);
    Sigmab=0.5*(Sigmab+Sigmab.t());
    Sigmab_inv=arma::inv_sympd(Sigmab);
    Sigmab_inv=0.5*(Sigmab_inv+Sigmab_inv.t());
    
    Sigmae2_old=Sigmae2;
    for(kk=0;kk<ndimy;kk++){
      Sigmae2(kk)=eee2_sum(kk)/double(arma::sum(nobsik.col(kk)));
      // xxitxxi_sum.fill(0.0);
      // xxityyi_sum.fill(0.0);
      // for(ii=0;ii<nind;ii++){
      //   xxitxxi_sum=xxitxxi_sum+F_xxik(ii,kk).t()*F_xxik(ii,kk);
      //   xxityyi_sum=xxityyi_sum+F_xxik(ii,kk).t()*(F_yik(ii,kk)-F_Phi(ii,kk)*F_Ebbi(ii).subvec(I0(kk),I1(kk)));
      //   // xxityyi_sum=xxityyi_sum+F_xxik(ii,kk).t()*(F_yik(ii,kk)-F_Phi(ii,kk)*mub.subvec(I0(kk),I1(kk)));
      // }
      // xi=arma::solve(xxitxxi_sum,xxityyi_sum);
      // F_xi(kk)=xi;
    }
    Sigmae2_inv=1.0/Sigmae2;
    
    // adjust
    for(kk=0;kk<ndimy;kk++){
      test_mat.zeros(nbas_long(kk)+ndimx,nbas_long(kk)+ndimx);
      test_vec.zeros(nbas_long(kk)+ndimx);
      for(ii=0;ii<nind;ii++){
        F_muhat(ii,kk)=F_Phi(ii,kk)*F_Ebbi(ii).subvec(I0(kk),I1(kk))+F_xxik(ii,kk)*F_xi(kk);
        test_mat.submat(0,0,nbas_long(kk)-1,nbas_long(kk)-1)+=F_Phi(ii,kk).t()*F_Phi(ii,kk);
        test_mat.submat(0,nbas_long(kk),nbas_long(kk)-1,nbas_long(kk)+ndimx-1)+=F_Phi(ii,kk).t()*F_xxik(ii,kk);
        test_mat.submat(nbas_long(kk),0,nbas_long(kk)+ndimx-1,nbas_long(kk)-1)+=F_xxik(ii,kk).t()*F_Phi(ii,kk);
        test_mat.submat(nbas_long(kk),nbas_long(kk),nbas_long(kk)+ndimx-1,nbas_long(kk)+ndimx-1)+=F_xxik(ii,kk).t()*F_xxik(ii,kk);
        test_vec.subvec(0,nbas_long(kk)-1)+=F_Phi(ii,kk).t()*F_muhat(ii,kk);
        test_vec.subvec(nbas_long(kk),nbas_long(kk)+ndimx-1)+=F_xxik(ii,kk).t()*F_muhat(ii,kk);
      }
      test_vec=arma::solve(test_mat,test_vec);
      mub.subvec(I0(kk),I1(kk))=test_vec.subvec(0,nbas_long(kk)-1);
      xi_old=F_xi(kk);
      xi=test_vec.subvec(nbas_long(kk),nbas_long(kk)+ndimx-1);
      F_xi(kk)=xi;
      F_abs_xi(kk)=arma::abs(xi-xi_old);
      F_rel_xi(kk)=F_abs_xi(kk)/(arma::abs(xi)+tol2);
    }
    
    abs_mub=arma::abs(mub-mub_old);
    rel_mub=abs_mub/(arma::abs(mub)+tol2);
    abs_Sigmab=arma::abs(Sigmab-Sigmab_old);
    rel_Sigmab=abs_Sigmab/(arma::abs(Sigmab)+tol2);
    abs_Sigmae2=arma::abs(Sigmae2-Sigmae2_old);
    rel_Sigmae2=abs_Sigmae2/(arma::abs(Sigmae2)+tol2);
    
    
    
    epsilon_abs_coef=0.0;
    epsilon_rel_coef=0.0;
    for(idx=0;idx<npossible_transitions;idx++){
      idx_begin=possible_transitions_begin(idx);
      idx_end=possible_transitions_end(idx);
      
      // if(possible_transitions(idx_begin,idx_end)!=1)continue;
      // Rcpp::Rcout<<"("<<idx_begin<<idx_end<<")";
      
      // yy part
      
      betay_new=F_betay(idx_begin,idx_end);
      betay_old=betay_new;
      // for(iter2=0;iter2<1;iter2++){
      SSy_new.fill(0.0);
      IIy_new.fill(0.0);
      for(ii=0;ii<nind;ii++){
        for(ll=0;ll<ntransition(ii);ll++){
          if(statebegin(ii,ll)!=idx_begin)continue;
          if(stateend(ii,ll)==-1)continue;
          tt=timeend(ii,ll);
          if(stateend(ii,ll)==idx_end){
            yytemp.zeros(ndimy);
            for(kk=0;kk<ndimy;kk++){
              yytemp(kk)=
                arma::dot(F_basismatrix_long(kk,tt),F_Ebbi(ii).subvec(I0(kk),I1(kk)))+
                arma::dot(xxi.row(ii),F_xi(kk));
            }
            SSy_new=SSy_new+yytemp;
          }
          denominatory=0.0;
          numerator1y.zeros(ndimy);
          numerator2y.zeros(ndimy,ndimy);
          for(ii2=0;ii2<nind;ii2++){
            for(ll2=0;ll2<ntransition(ii2);ll2++){
              if(statebegin(ii2,ll2)!=idx_begin)continue;
              if(timebegin(ii2,ll2)>=tt)continue;
              if(timeend(ii2,ll2)<tt)continue;
              for(idx2=0;idx2<npossible_transitions;idx2++){
                if(possible_transitions_begin(idx2)!=idx_begin)continue;
                idx_begin2=possible_transitions_begin(idx2);
                idx_end2=possible_transitions_end(idx2);
                if(idx_end==idx_end2){
                  denominatory=denominatory+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
                  numerator1y=numerator1y+F_Eyexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
                  numerator2y=numerator2y+F_Eyyexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
                }else{
                  denominatory=denominatory+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
                }
              }
            }
          }
          SSy_new=SSy_new-numerator1y/denominatory;
          IIy_new=IIy_new-numerator2y/denominatory+
            (numerator1y/denominatory)*(numerator1y/denominatory).t();
        }
      }
      betay_new=betay_new-factor*arma::solve(IIy_new,SSy_new);
      // }
      F_betay(idx_begin,idx_end)=betay_new;
      F_SSy(idx_begin,idx_end)=SSy_new;
      F_IIy(idx_begin,idx_end)=IIy_new;
      F_abs_betay(idx_begin,idx_end)=arma::abs(betay_new-betay_old);
      F_rel_betay(idx_begin,idx_end)=F_abs_betay(idx_begin,idx_end)/(arma::abs(betay_new)+tol2);
      
      // zz part
      
      betaz_new=F_betaz(idx_begin,idx_end);
      betaz_old=betaz_new;
      // for(iter2=0;iter2<1;iter2++){
      SSz_new.fill(0.0);
      IIz_new.fill(0.0);
      for(ii=0;ii<nind;ii++){
        for(ll=0;ll<ntransition(ii);ll++){
          if(statebegin(ii,ll)!=idx_begin)continue;
          if(stateend(ii,ll)==-1)continue;
          tt=timeend(ii,ll);
          if(stateend(ii,ll)==idx_end){
            temp_vec=F_zzi(ii);
            SSz_new=SSz_new+temp_vec;
          }
          denominatorz=0.0;
          numerator1z.zeros(ndimz);
          numerator2z.zeros(ndimz,ndimz);
          for(ii2=0;ii2<nind;ii2++){
            for(ll2=0;ll2<ntransition(ii2);ll2++){
              if(statebegin(ii2,ll2)!=idx_begin)continue;
              if(timebegin(ii2,ll2)>=tt)continue;
              if(timeend(ii2,ll2)<tt)continue;
              for(idx2=0;idx2<npossible_transitions;idx2++){
                if(possible_transitions_begin(idx2)!=idx_begin)continue;
                idx_begin2=possible_transitions_begin(idx2);
                idx_end2=possible_transitions_end(idx2);
                if(idx_end==idx_end2){
                  temp_vec=F_zzi(ii2);
                  denominatorz=denominatorz+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
                  numerator1z=numerator1z+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2))*temp_vec;
                  numerator2z=numerator2z+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2))*temp_vec*temp_vec.t();
                }else{
                  denominatorz=denominatorz+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
                }
              }
            }
          }
          SSz_new=SSz_new-numerator1z/denominatorz;
          IIz_new=IIz_new-numerator2z/denominatorz+
            (numerator1z/denominatorz)*(numerator1z/denominatorz).t();
        }
      }
      betaz_new=betaz_new-factor*arma::solve(IIz_new,SSz_new);
      // }
      F_betaz(idx_begin,idx_end)=betaz_new;
      F_SSz(idx_begin,idx_end)=SSz_new;
      F_IIz(idx_begin,idx_end)=IIz_new;
      F_abs_betaz(idx_begin,idx_end)=arma::abs(betaz_new-betaz_old);
      F_rel_betaz(idx_begin,idx_end)=F_abs_betaz(idx_begin,idx_end)/(arma::abs(betaz_new)+tol2);
      
      // alpha part
      
      if(F_ref(idx_begin,idx_end)==1){
        alpha_new=0.0;
        alpha_old=0.0;
      }else{
        alpha_new=F_alpha(idx_begin,idx_end);
        alpha_old=alpha_new;
        // for(iter2=0;iter2<1;iter2++){
        SSa_new=0.0;
        IIa_new=0.0;
        for(ii=0;ii<nind;ii++){
          for(ll=0;ll<ntransition(ii);ll++){
            if(statebegin(ii,ll)!=idx_begin)continue;
            if(stateend(ii,ll)==-1)continue;
            tt=timeend(ii,ll);
            if(stateend(ii,ll)==idx_end){
              SSa_new=SSa_new+1.0;
            }
            denominatora=0.0;
            numerator1a=0.0;
            numerator2a=0.0;
            for(ii2=0;ii2<nind;ii2++){
              for(ll2=0;ll2<ntransition(ii2);ll2++){
                if(statebegin(ii2,ll2)!=idx_begin)continue;
                if(timebegin(ii2,ll2)>=tt)continue;
                if(timeend(ii2,ll2)<tt)continue;
                for(idx2=0;idx2<npossible_transitions;idx2++){
                  if(possible_transitions_begin(idx2)!=idx_begin)continue;
                  idx_begin2=possible_transitions_begin(idx2);
                  idx_end2=possible_transitions_end(idx2);
                  if(idx_end==idx_end2){
                    denominatora=denominatora+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
                    numerator1a=numerator1a+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
                    numerator2a=numerator2a+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
                  }else{
                    denominatora=denominatora+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
                  }
                }
              }
            }
            SSa_new=SSa_new-numerator1a/denominatora;
            IIa_new=IIa_new-numerator2a/denominatora+
              (numerator1a/denominatora)*(numerator1a/denominatora);
          }
        }
        alpha_new=alpha_new-factor*SSa_new/IIa_new;
        // }
      }
      F_alpha(idx_begin,idx_end)=alpha_new;
      F_SSa(idx_begin,idx_end)=SSa_new;
      F_IIa(idx_begin,idx_end)=IIa_new;
      F_abs_alpha(idx_begin,idx_end)=std::abs(alpha_new-alpha_old);
      F_rel_alpha(idx_begin,idx_end)=F_abs_alpha(idx_begin,idx_end)/(std::abs(alpha_new)+tol2);
      
      epsilon_abs_coef=std::max(epsilon_abs_coef,F_abs_betay(idx_begin,idx_end).max());
      epsilon_abs_coef=std::max(epsilon_abs_coef,F_abs_betaz(idx_begin,idx_end).max());
      epsilon_abs_coef=std::max(epsilon_abs_coef,F_abs_alpha(idx_begin,idx_end));
      
      epsilon_rel_coef=std::max(epsilon_rel_coef,F_rel_betay(idx_begin,idx_end).max());
      epsilon_rel_coef=std::max(epsilon_rel_coef,F_rel_betaz(idx_begin,idx_end).max());
      epsilon_rel_coef=std::max(epsilon_rel_coef,F_rel_alpha(idx_begin,idx_end));
    }
    
    epsilon_abs=epsilon_abs_coef;
    epsilon_abs=std::max(epsilon_abs,abs_mub.max());
    epsilon_abs=std::max(epsilon_abs,abs_Sigmab.diag().max());
    epsilon_abs=std::max(epsilon_abs,abs_Sigmae2.max());
    epsilon_rel=epsilon_rel_coef;
    epsilon_rel=std::max(epsilon_rel,rel_mub.max());
    epsilon_rel=std::max(epsilon_rel,rel_Sigmab.diag().max());
    epsilon_rel=std::max(epsilon_rel,rel_Sigmae2.max());
    // Rcpp::Rcout<<"rel_mub"<<rel_mub.max()<<"rel_Sigmab"<<rel_Sigmab.max()<<"rel_Sigmae2"<<rel_Sigmae2.max()<<std::endl;
    for(kk=0;kk<ndimy;kk++){
      epsilon_abs=std::max(epsilon_abs,F_abs_xi(kk).max());
      epsilon_rel=std::max(epsilon_rel,F_rel_xi(kk).max());
    }
    
    if(epsilon_rel_coef<tol)break;
    timing=timer.toc();
    // Rcpp::Rcout<<"  epsrel: "<<epsilon_rel<<", epsabs: "<<epsilon_abs<<")";
    Rcpp::Rcout<<"  epsrelcoef: "<<epsilon_rel_coef<<", epsabscoef: "<<epsilon_abs_coef<<")";
    Rcpp::Rcout<<"  timing: "<<timing<<")"<<std::endl;
    
    for(ii=0;ii<nind;ii++){
      for(ll=0;ll<ntransition(ii);ll++){
        idx_begin=statebegin(ii,ll);
        idx_end=stateend(ii,ll);
        if(idx_end==-1)continue;
        tt=timeend(ii,ll);
        
        sum_risk=0.0;
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
              sum_risk=sum_risk+F_Eexp(ii2,idx2,tt)*std::exp(F_alpha(idx_begin2,idx_end2));
            }
          }
        }
        F_hazard0(idx_begin)(tt)=1.0/sum_risk;
      }
    }
    
    // output csv
    if(csv){
      mub_csv<<iter<<" ";
      mub.t().raw_print(mub_csv<<std::scientific);
      temp_mat=Sigmab;temp_mat.insert_cols(0,double(iter)*arma::ones(ndim_rand));
      temp_mat.raw_print(Sigmab_csv<<std::scientific);
      Sigmae2_csv<<iter<<" ";
      Sigmae2.t().raw_print(Sigmae2_csv<<std::scientific);
      for(kk=0;kk<ndimy;kk++){
        F_xi_csv<<iter<<" "<<kk<<" ";
        F_xi(kk).t().raw_print(F_xi_csv<<std::scientific);
      }
      for(idx=0;idx<npossible_transitions;idx++){
        idx_begin=possible_transitions_begin(idx);
        idx_end=possible_transitions_end(idx);
        F_betay_csv<<iter<<" "<<idx_begin<<" "<<idx_end<<" ";
        F_betay(idx_begin,idx_end).t().raw_print(F_betay_csv<<std::scientific);
        F_betaz_csv<<iter<<" "<<idx_begin<<" "<<idx_end<<" ";
        F_betaz(idx_begin,idx_end).t().raw_print(F_betaz_csv<<std::scientific);
        F_alpha_csv<<iter<<" "<<idx_begin<<" "<<idx_end<<" ";
        F_alpha_csv<<std::scientific<<F_alpha(idx_begin,idx_end)<<std::endl;
      }
    }
  }
  
  if(csv){
    mub_csv.close();
    Sigmab_csv.close();
    Sigmae2_csv.close();
    F_xi_csv.close();
    F_betay_csv.close();
    F_betaz_csv.close();
    F_alpha_csv.close();
  }
  
  arma::field<arma::vec> F_betay_flat(npossible_transitions);
  arma::field<arma::vec> F_betaz_flat(npossible_transitions);
  for(idx=0;idx<npossible_transitions;idx++){
    idx_begin=possible_transitions_begin(idx);
    idx_end=possible_transitions_end(idx);
    F_betay_flat(idx)=F_betay(idx_begin,idx_end);
    F_betaz_flat(idx)=F_betaz(idx_begin,idx_end);
  }
  
  List result(30);
  
  result(0)=mub;
  result(1)=Sigmab;
  result(2)=Sigmae2;
  result(3)=F_alpha;
  result(4)=F_betay_flat;
  result(5)=F_betaz_flat;
  result(6)=F_hazard0;
  result(7)=F_hazard0_count;
  result(8)=F_xi;
  
  result(10)=F_betay;
  result(11)=F_betaz;
  result(12)=F_Ebbi;
  
  result(20)=iter;
  result(21)=epsilon_abs;
  result(22)=epsilon_rel;
  
  return(result);
}
