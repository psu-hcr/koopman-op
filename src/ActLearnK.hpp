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
  int T_index;
      
  alk(system *_sys, objective *_cost,policy *_pol, double _T,const arma::vec& _umax,const arma::mat& _R){
    sys = _sys; cost=_cost; pol=_pol;T=_T;umax = _umax;
    T_index = T/sys->dt; Rinv = _R.i(); 
    };
    
  arma::vec ustar_calc();//main function for calculating the current action
  arma::mat xforward();//forward simulation of x
  arma::mat rhoback(const arma::mat& xsol); //backward simulation of the adjoint
  inline arma::vec f(const arma::vec& rho, xupair pair){
    arma::vec x = pair.x; arma::vec u = pair.u; double ti = pair.t;
  	return -(sys->dfdx(x,u)+sys->hx(x)*pol->dmudz(x,ti)).t()*rho
		-(cost->dldz(x,u,ti)+pol->dmudz(x,ti).t()*pol->dldu(x,u,ti)) ;
	}//f for rho backwards sim
	};

//main function for calculating a single control vector
template <class system, class objective,class policy>
arma::vec alk<system,objective,policy>::ustar_calc(){ 
  arma::vec ustar;
  ustar = arma::zeros<arma::vec>(size(sys->Ucurr));
  arma::mat xsol,rhosol;    
  xsol = xforward();
  rhosol = rhoback(xsol);  
  ustar = -Rinv*(sys->Ku*sys->zfuncs->dvdu(sys->Xcurr)).t()*rhosol.col(0)
  		+pol->mu(sys->Xcurr,sys->tcurr);
  
return saturation(ustar);}
    
//forward simulation of x
template <class system, class objective,class policy>
arma::mat alk<system,objective,policy>::xforward(){
  arma::mat xsol = arma::zeros<arma::mat>(sys->Xcurr.n_rows,T_index);
  arma::vec x0 = sys->Xcurr; double ti=sys->tcurr;
  xsol.col(0)=x0;
  for(int i = 1; i<T_index;i++){
 	x0 = x0+sys->f(x0,pol->mu(x0,ti))*sys->dt;
	ti = sys->tcurr+(double)i*sys->dt;
	xsol.col(i)=x0;
    //x0 = RK4_step<system,const arma::vec&>(sys,x0,pol->mu(x0,ti),sys->dt);
  }    
return xsol;}

//backward simulaiton of the adjoint
template <class system, class objective,class policy>
arma::mat alk<system,objective,policy>::rhoback(const arma::mat& xsol){
  arma::mat rhosol = arma::zeros<arma::mat>(sys->Xcurr.n_rows,T_index);
  arma::vec rho0 = sys->Xcurr;
  xupair current;
  rho0.zeros();
  rhosol.col(T_index-1)=rho0;
  for(int i = T_index-2; i>=0;i--){
    current.x =xsol.col(i);
    current.t = sys->tcurr+(double)i*sys->dt;
	current.u = pol->mu(current.x,current.t);
	//rho0 = RK4_step<alk,xupair>(this,rho0,current,-1.0*sys->dt);
	rho0 = rho0-f(rho0,current)*sys->dt;
	rhosol.col(i)=rho0;
  } 
return rhosol;}


#endif