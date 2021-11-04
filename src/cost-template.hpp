#ifndef templateCOST_HPP
#define templateCOST_HPP
#include<armadillo>

template <class system>
class templatecost {
    system* sys;
  public:
    emplatetcost(system *_sys){
      // initialize with sys, weight matrices, and reference
      };
    inline double l (const arma::vec& x,const arma::vec& u,double ti){
      arma::vec xproj = sys->proj_func(x);
      return arma::as_scalar();//not required for SAC but useful in most cases
      }
    inline arma::vec dldx (const arma::vec& x,const arma::vec& u,double ti){//REQUIRED
      arma::vec xproj = sys->proj_func(x);
      return ;//return the derivative of the incremental cost at a particular time and state
      }
	  inline arma::vec dldu(const arma::vec& x,const arma::vec& u,double ti){//REQUIRED for active learning only
      arma::vec xproj = sys->proj_func(x);
      return ;//return the derivative of the incremental cost at a particular time and state
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