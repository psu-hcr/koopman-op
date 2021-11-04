#ifndef FISHERCOST_HPP
#define FISHERCOST_HPP
#include<armadillo>
#include <math.h>
#include"src/koopsys.hpp"

template <class system>//system type must be KoopSys
class fishcost {
    system* sys;
	arma::mat Sigma;
	double eps = 1./10000.;
  public:
    fishcost(system *_sys, const arma::vec SigDiag){
		sys=_sys;
		Sigma = arma::diagmat(SigDiag);

      // initialize with sys, weight matrices, and reference
      };
    inline double l (const arma::vec& x,const arma::vec& u,double ti){
      arma::vec xproj = sys->proj_func(x);
	  arma::vec dfdk = sys->zfuncs->zx(x,u);
	  double fishTopt = arma::as_scalar(dfdk.t()*Sigma*dfdk);
      return 1./(fishTopt+eps);//not required for SAC but useful in most cases
      }
    inline arma::vec dldx (const arma::vec& x,const arma::vec& u,double ti){//REQUIRED
      arma::vec xproj = sys->proj_func(x);
	  arma::vec dfdk = sys->zfuncs->zx(x,u);
	  double fishTopt = arma::as_scalar(dfdk.t()*Sigma*dfdk);
	  
      return (-2./pow(fishTopt+eps,2))*dfdk.t()*Sigma;//return the derivative of the incremental cost at a particular time and state
      }
    double calc_cost(const arma::mat& x,const arma::mat& u){//REQUIRED:calculate the total cost
      arma::vec xproj;
      double J1 = 0.0;//may need calculate cost of barrier functions before main l calculation
      for (int i = 0; i<x.n_cols; i++){
        xproj = sys->proj_func(x.col(i));
        J1+=l(xproj,u.col(i),sys->tcurr+(double)i*sys->dt);
        }
      return J1;
      }
    
};

#endif