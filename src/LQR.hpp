#ifndef LQR_HPP
#define LQR_HPP
#include<armadillo>
#include <iostream>

class lqr {
  arma::mat A, B, Q, R, Rinv,Qf,K;//qf=arma::zeros()
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
  lqr(const arma::mat& _Q,const arma::mat& _R,const arma::mat& _Qf,int _horizon,const arma::vec& _umax, double _dt){
    Q=_Q; R=_R; Rinv = _R.i(); Qf=_Qf; 
	horizon = _horizon; umax = _umax; dt = _dt;
    };
  void set_xd(const arma::vec& _xd){xd=_xd;};
  void calc_gains(const arma::mat& _A,const arma::mat& _B);
  arma::vec mu(const arma::vec& _x);
  inline arma::mat f(const arma::mat& Pvec, int NA){
	  arma::mat P = arma::reshape(Pvec,arma::size(Qf));
	  arma::mat Pdot = A.t()*P + P*A - P*B*Rinv*B.t()*P + Q;
	  return Pdot.as_col();
  	}//f for backwards sim in Ricatti eqtn
};


void lqr::calc_gains(const arma::mat& _A,const arma::mat& _B){ 
	A = _A; B=_B;//update Kx, Ku before recalculating gain
	arma::vec Pflat = Qf.as_col(); int na = 1.;
   for(int i = horizon-1;i>=0;i--){
	   Pflat = RK4_step<lqr,int>(this,Pflat,na,-1.0*dt);
   	}
	arma::mat P = arma::reshape(Pflat,arma::size(Qf));
	K=Rinv*B.t()*P;
return;}

arma::vec lqr::mu(const arma::vec& _x){
	return K*(_x-xd);
}

#endif