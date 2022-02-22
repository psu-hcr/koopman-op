#ifndef ERRORCOST_HPP
#define ERRORCOST_HPP
#include<armadillo>
#include <math.h>

template <class system>
class errorcostfl {
    system* sys;
    double eps = pow(10,-5);
    double infw = 100.;
    
  public:
    arma::mat Q;
    arma::mat R;
    std::function<arma::vec(double)> xd;
	arma::mat Sigma;
    
    errorcostfl(arma::mat _Q, arma::mat _R,std::function<arma::vec(double)> _xd, system *_sys, const arma::vec SigDiag){
      Q=_Q; R=_R; sys=_sys; xd = _xd;// initialize with Q, R, sys, xd
      Sigma = arma::diagmat(SigDiag);// initialize with Sigma
      };
	  
    inline double l (const arma::vec& x,const arma::vec& u,double ti){
      arma::vec xproj = sys->proj_func(x);
      arma::vec dfdk = sys->zfuncs->zx(x);
	  double fishTopt = arma::as_scalar(dfdk.t()*Sigma*dfdk);
      return arma::as_scalar(((xproj.t()-xd(ti).t())*Q*(xproj-xd(ti))+u.t()*R*u)/2) + (infw*1.)/(fishTopt+eps);
      }
	  
    inline arma::vec dldx (const arma::vec& x,const arma::vec& u,double ti){
      arma::vec xproj = sys->proj_func(x);
      arma::vec dfdk = sys->zfuncs->zx(x);
	  double fishTopt = arma::as_scalar(dfdk.t()*Sigma*dfdk);
      return Q*(xproj-xd(ti)) + ((-2.*infw)/pow(fishTopt+eps,2))*Sigma*dfdk;
      }
	  
    double calc_cost(const arma::mat& x,const arma::mat& u){
      arma::vec xproj;
      double J1 = 0.0;
      for (int i = 0; i<x.n_cols; i++){
        xproj = sys->proj_func(x.col(i));
        J1+=l(xproj,u.col(i),sys->tcurr+(double)i*sys->dt);
        }
      return J1;
      }
    
};

#endif