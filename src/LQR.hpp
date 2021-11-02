#ifndef LQR_HPP
#define LQR_HPP
#include<armadillo>
#include <iostream>

//template <class system, class objective>
class lqr {
  arma::mat A, B, Q, R, Rinv,Qf;//qf=arma::zeros()
  int horizon;// horizon = 20
  arma::vec xd;
  arma::vec umax;
  void calc_gains();
  double dt;
  //saturation function for the control        
  arma::vec saturation(const arma::vec& u){
    arma::vec usat; usat.zeros(u.n_rows);
    for (int i = 0; i<u.n_rows; i++){
      if(u(i) > umax(i)) usat(i) = umax(i);
        else if(u(i) < -umax(i)){ usat(i) = -umax(i);}
          else usat(i) = u(i);};
    return usat;}
      
  public: 
  lqr(const arma::mat& _A,const arma::mat& _B,const arma::mat& _Q,const arma::mat& _R,const arma::mat& _Qf,int _horizon,const arma::vec& _umax, double _dt){
    A=_A; B=_B; Q=_Q; R=_R; Rinv = _R.i(); Qf=_Qf; horizon = _horizon; umax = _umax; dt = _dt;
    
  };
  void set_xd(const arma::vec& xd);
     
  inline arma::mat f(const arma::mat& P, int NA){
	  return A.t()*P + P*A - P*B*Rinv*B.t()*P + Q;
  	}//f for backwards sim in Ricatti eqtn
};


void lqr::calc_gains(){ 
	arma::mat P = Qf; int na = 1.;
   for(int i = horizon-1;i>=0;i--){
	   P = RK4_step<lqr,int>(this,P,na,-1.0*dt);
   	}
	arma::mat K=Rinv*B.t()*P;
return;}

#endif