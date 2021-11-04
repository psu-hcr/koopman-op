#ifndef LQR_HPP
#define LQR_HPP
#include<armadillo>
#include <iostream>

class lqr {
  arma::mat A, B, Q, R, Rinv,Qf,K;//qf=arma::zeros()
  int horizon;// horizon = 20
  std::function<arma::vec(double)> xd;
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
  lqr(const arma::mat& _Q,const arma::mat& _R,const arma::mat& _Qf,int _horizon,const arma::vec& _umax, std::function<arma::vec(double)> _xd, double _dt){
    Q=_Q; R=_R; Rinv = _R.i(); Qf=_Qf; 
	horizon = _horizon; umax = _umax; xd=_xd, dt = _dt;
    };
  
  void calc_gains(const arma::mat& _A,const arma::mat& _B);
  arma::vec mu(const arma::vec& _x,double ti);
  arma::vec dmudz(const arma::vec& _x);
  inline arma::mat f(const arma::mat& Pvec, int NA){
	  arma::mat P = arma::reshape(Pvec,arma::size(Qf));
	  arma::mat Pdot = A.t()*P + P*A - P*B*Rinv*B.t()*P + Q;
	  return Pdot.as_col();
  	}//f for backwards sim in Ricatti eqtn
	//function below required for active learning integration with this policy
  double l (const arma::vec& x,const arma::vec& u,double ti){
  return arma::as_scalar((x-xd(ti)).t()*Q*(x-xd(ti))+u.t()*R*u);}
  arma::vec dldx (const arma::vec& x,const arma::vec& u,double ti){
    return Q*(x-xd(ti));}
  arma::vec dldu (const arma::vec& x,const arma::vec& u,double ti){
    return R*u;}
  
  
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

arma::vec lqr::mu(const arma::vec& _x, double ti){
	return K*(_x-xd(ti));
}
arma::vec lqr::dmudz(const arma::vec& _x){
  return K;
}

#endif