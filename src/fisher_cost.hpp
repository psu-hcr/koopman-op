#ifndef FISHERCOST_HPP
#define FISHERCOST_HPP
#include<armadillo>
#include <math.h>
//#include"src/koopsys.hpp"

template <class system, class objective>//system type must be KoopSys<>
class fishcost {
    system* sys;//uses dfdx,zfuncs,hx
	objective* task;//uses task->l, task->dldx,task->dldu,task->dmudz
	arma::mat Sigma;
	double eps = pow(10,-5);
	double infw = 100.;
  public:
    fishcost(system *_sys, objective *_task,const arma::vec SigDiag){
		sys=_sys;task = _task;
		Sigma = arma::diagmat(SigDiag);

      // initialize with sys, weight matrices, and reference
      };
    inline double l (const arma::vec& x,const arma::vec& u,double ti){
      arma::vec dfdk = sys->zfuncs->zx(x,u);
	  double fishTopt = arma::as_scalar(dfdk.t()*Sigma*dfdk);
      return (infw*1.)/(fishTopt+eps)+task->l(x,u,ti);
      }
    inline arma::vec dldx (const arma::vec& x,const arma::vec& u,double ti){//REQUIRED
       arma::vec dfdk = sys->zfuncs->zx(x);
	  double fishTopt = arma::as_scalar(dfdk.t()*Sigma*dfdk);
	  arma::vec dldz = ((-2.*infw)/pow(fishTopt+eps,2))*Sigma*dfdk+task->dldx(x,u,ti);
	  return dldz+task->dmudz(x).t()*task->dldu(x,u,ti);
      }
	inline arma::vec dldu (const arma::vec& x,const arma::vec& u,double ti){//REQUIREd 
    return task->dldu(x,u,ti);
      }
	inline arma::mat dfdx (const arma::vec& x,const arma::vec& u){//required
		
		return sys->dfdx(x,u)+sys->hx(x)*task->dmudz(x);}
    double calc_cost(const arma::mat& x,const arma::mat& u){//REQUIRED:calculatetotal cost
      
      double J1 = 0.0;
      for (int i = 0; i<x.n_cols; i++){
        
        J1+=l(x,u.col(i),sys->tcurr+(double)i*sys->dt);
        }
      return J1;
      }
    
};

#endif