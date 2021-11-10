#ifndef ALK_HPP
#define ALK_HPP
#include<armadillo>
#include <iostream>
#include"xustruct.hpp"


template <class system, class objective,class policy>
class alk {
  system* sys; //from sys use sys->f, sys->proj_func, sys->dfdx, sys->hx,sys->dt
  policy* pol;
  std::function<arma::vec(double)> unom;
  double T;  
  arma::vec umax;
  arma::mat Rinv;
  //saturation function for the control        
  arma::vec saturation(const arma::vec& u){
    arma::vec usat; usat.zeros(u.n_rows);
    for (int i = 0; i<u.n_rows; i++){
      if(u(i) > umax(i)) usat(i) = umax(i);
        else if(u(i) < -umax(i)){ usat(i) = -umax(i);}
          else usat(i) = u(i);};
    return usat;}
  
    
  public:
  objective* cost; //from cost use cost->dldx, cost->calc_cost,cost->dfdx
  //bool iterative=false;
  int T_index;
  //arma::mat ulist;
    
  alk(system *_sys, objective *_cost,policy *_pol, double _T,const arma::vec& _umax,const arma::mat& _R){
    sys = _sys; cost=_cost; pol=_pol;T=_T;umax = _umax;
    T_index = T/sys->dt; Rinv = _R.i();
    //ulist = arma::zeros(umax.n_rows,T_index); 
    //for(int i = 0;i<ulist.n_cols-1;i++) ulist.col(i) = unom(sys->tcurr + (double)i*sys->dt);
  };
    
  void ustar_calc();//main function for calculating the current action
  arma::mat xforward();//forward simulation of x
  arma::mat rhoback(const arma::mat& xsol); //backward simulation of the adjoint
  inline arma::vec f(const arma::vec& rho, xupair pair){
  	return -cost->dldx(pair.x,pair.u,pair.t);}// - cost->dfdx(pair.x,pair.u).t()*rho;}//f for rho backwards sim
	};

//main function for calculating a single control vector
template <class system, class objective,class policy>
void alk<system,objective,policy>::ustar_calc(){ 
  arma::vec ustar;
  ustar = arma::zeros<arma::vec>(size(sys->Ucurr));
  arma::mat xsol,rhosol;    
  xsol = xforward();
  rhosol = rhoback(xsol);  
  ustar = -Rinv*(sys->Ku*sys->zfuncs->dvdu(sys->Xcurr)).t()*rhosol.col(0)+pol->mu(sys->Xcurr,sys->tcurr);
  ustar=saturation(ustar);
return;}
    
//forward simulation of x
template <class system, class objective,class policy>
arma::mat alk<system,objective,policy>::xforward(){
  arma::mat xsol = arma::zeros<arma::mat>(sys->Xcurr.n_rows,T_index);
  arma::vec x0 = sys->Xcurr; double ti;
  for(int i = 0; i<T_index;i++){
  	ti = sys->tcurr+(double)i*sys->dt;
    xsol.col(i)=x0;
    x0 = RK4_step<system,const arma::vec&>(sys,x0,pol->mu(x0,ti),sys->dt);
  }    
return xsol;}

//backward simulaiton of the adjoint
template <class system, class objective,class policy>
arma::mat alk<system,objective,policy>::rhoback(const arma::mat& xsol){
  arma::mat rhosol = arma::zeros<arma::mat>(sys->Xcurr.n_rows,T_index);
  arma::vec rho0 = sys->Xcurr;
  xupair current;
  rho0.zeros();
  for(int i = T_index-1; i>=0;i--){
    rhosol.col(i)=rho0;
    current.x =xsol.col(i);
    current.t = sys->tcurr+(double)i*sys->dt;
	current.u = pol->mu(current.x,current.t);
	rho0 = RK4_step<alk,xupair>(this,rho0,current,-1.0*sys->dt);
  } 
return rhosol;}


#endif