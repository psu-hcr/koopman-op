#ifndef LQR_HPP
#define LQR_HPP
#include<armadillo>
#include <iostream>

class lqr {
  arma::mat A, B, Q, R, Rinv,Qf;//qf=arma::zeros()
  int horizon;// horizon = 20
  double tcurr;
  
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
	std::function<arma::vec(double)> xd;
  	lqr(const arma::mat& _Q,const arma::mat& _R,const arma::mat& _Qf,int _horizon,const arma::vec& _umax, std::function<arma::vec(double)> _xd, double _dt){
    	Q=_Q; R=_R; Rinv = _R.i(); Qf=_Qf; 
		horizon = _horizon; umax = _umax; xd=_xd, dt = _dt;
		xdim = xd(0).n_rows;udim = umax.n_rows;
		Klist.zeros(udim*(horizon),xdim);
    };
  
  arma::mat K(double ti);
  arma::mat Klist;int xdim,udim;
  void calc_gains(const arma::mat& _A,const arma::mat& _B, double tc);
  arma::vec mu(const arma::vec& _x,double ti);
  arma::mat dmudz(const arma::vec& _x,double ti);
  inline arma::mat f(const arma::mat& Pvec, int NA){
	  arma::mat P = arma::reshape(Pvec,arma::size(Qf));
	  arma::mat Pdot = -(A.t()*P + P*A - ((P*B)*Rinv)*(B.t()*P) + Q); //std::cout<<Pdot;
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


void lqr::calc_gains(const arma::mat& _A,const arma::mat& _B, double _tcurr){ 
	A = _A; B=_B;//update Kx, Ku before recalculating gain
	tcurr = _tcurr;//update current time to get correct K from list
	arma::vec Pflat = Qf.as_col(); int na = 1.;
	arma::mat P;
   for(int i = horizon;i>0;i--){
	   P = arma::reshape(Pflat,arma::size(Qf)); 
	   Klist.submat(udim*i-udim,0,udim*i-1,xdim-1)=Rinv*B.t()*P;  
	   //Pflat = RK4_step<lqr,int>(this,Pflat,na,-1.0*dt);
	   Pflat = Pflat - f(Pflat,i)*dt;
	   }
return;}

arma::mat lqr::K(double ti){
	int i = round((ti-tcurr)/dt);
	arma::mat Ki = Klist.submat(udim*i,0,udim*i+(udim-1),xdim-1);
	return Ki;}

arma::vec lqr::mu(const arma::vec& _x, double ti){
	return saturation(-K(tcurr)*(_x-xd(ti)));}
arma::mat lqr::dmudz(const arma::vec& _x,double ti){
  return -K(ti);
  }

#endif