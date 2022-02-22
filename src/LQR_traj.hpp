#ifndef LQR_TRAJ_HPP
#define LQR_TRAJ_HPP
#include<armadillo>
#include <iostream>

template <class trajectory>
class LQR_traj {
  	arma::mat A, B, Q, R, Rinv,Qf;//qf=arma::zeros()
  	int horizon;// horizon = 20
  	double tcurr;
  	arma::vec umax;
  	void calc_gains();
  	double dt;
	
  	// saturation function for the control        
  	arma::vec saturation(const arma::vec& u){
    	arma::vec usat; usat.zeros(u.n_rows);
    	for (int i = 0; i<u.n_rows; i++){
      		if(u(i) > umax(i)) usat(i) = umax(i);
        	else if(u(i) < -umax(i)){ usat(i) = -umax(i);}
          	else usat(i) = u(i);
		};
    	return usat;
	};
      
	public: 
	trajectory* traj;
	arma::mat K(double ti);
  	arma::mat Klist;int xdim,udim;
	void calc_gains(const arma::mat& _A,const arma::mat& _B, double tc);
  	arma::vec mu(const arma::vec& _x,double ti);
  	arma::mat dmudz(const arma::vec& _x,double ti);
	
  	LQR_traj(const arma::mat& _Q,const arma::mat& _R,const arma::mat& _Qf,int _horizon,const arma::vec& _umax, trajectory* _traj, double _dt){
		// initialize twoLink class
    	Q=_Q; R=_R; Rinv = _R.i(); Qf=_Qf; 
		horizon = _horizon; umax = _umax; traj = _traj; dt = _dt;
		xdim = traj->desire_jointstate(0).n_rows;udim = umax.n_rows;
		Klist.zeros(udim*(horizon),xdim);
    };
  
  	inline arma::mat f(const arma::mat& Pvec, int NA){
	  	arma::mat P = arma::reshape(Pvec,arma::size(Qf));
		//std::cout<<"resize P"<<std::endl;
		//std::cout<<P<<std::endl;
	  	arma::mat Pdot = -(A.t()*P + P*A - ((P*B)*Rinv)*(B.t()*P) + Q);
		//std::cout<<"A.t()*P"<<std::endl;
		//std::cout<<A.t()*P<<std::endl;
		//std::cout<<"P*A"<<std::endl;
		//std::cout<<P*A<<std::endl;
		//std::cout<<"((P*B)*Rinv)*(B.t()*P)"<<std::endl;
		//std::cout<<((P*B)*Rinv)*(B.t()*P)<<std::endl;
		//std::cout<<"Q"<<std::endl;
		//std::cout<<Q<<std::endl;
		//std::cout<<"Pdot"<<std::endl;
		//std::cout<<Pdot<<std::endl;
	  	return Pdot.as_col();
  	}//f for backwards sim in Ricatti eqtn
	//function below required for active learning integration with this policy
	
  	double l (const arma::vec& x,const arma::vec& u,double ti){
		arma::vec xd = traj->desire_jointstate(ti);
  		return arma::as_scalar((x-xd).t()*Q*(x-xd)+u.t()*R*u);
	}
	
  	arma::vec dldx (const arma::vec& x,const arma::vec& u,double ti){
		arma::vec xd = traj->desire_jointstate(ti);
    	return Q*(x-xd);
	}
	
  	arma::vec dldu (const arma::vec& x,const arma::vec& u,double ti){
    	return R*u;
	}
};

template <class trajectory>
void LQR_traj<trajectory>::calc_gains(const arma::mat& _A,const arma::mat& _B, double _tcurr){ 
	A = _A; B = _B;	// update Kx, Ku before recalculating gain
	tcurr = _tcurr;	// update current time to get correct K from list
	arma::vec Pflat = Qf.as_col(); int na = 1.;
	arma::mat P;
   	for(int i = horizon;i>0;i--){
		//std::cout<<"i "<<i<<std::endl;
	   	P = arma::reshape(Pflat,arma::size(Qf));
	   	//std::cout<<"P"<<std::endl;
	   	//std::cout<<P<<std::endl;
	   	Klist.submat(udim*i-udim,0,udim*i-1,xdim-1)=Rinv*B.t()*P;
	   	//std::cout<<"f(Pflat,i)"<<std::endl;
	   	//std::cout<<f(Pflat,i)<<std::endl; 
		//Pflat = RK4_step<lqr,int>(this,Pflat,na,-1.0*dt);
	   	Pflat = Pflat - f(Pflat,i)*dt; 
	}
	return;
}

template <class trajectory>
arma::mat LQR_traj<trajectory>::K(double ti){
	int i = round((ti-tcurr)/dt);
	//std::cout<<"LQR_traj::K i "<<i<<std::endl;
	//std::cout<<"udim"<<udim<<std::endl;
	//std::cout<<"size(Klist)"<<size(Klist)<<endl;
	arma::mat Ki = Klist.submat(udim*i,0,udim*i+(udim-1),xdim-1);	
	return Ki;
}

template <class trajectory>
arma::vec LQR_traj<trajectory>::mu(const arma::vec& _x, double ti){
	arma::vec xd = traj->desire_jointstate(ti);
	return saturation(-K(tcurr)*(_x-xd));
}

template <class trajectory>
arma::mat LQR_traj<trajectory>::dmudz(const arma::vec& _x,double ti){
  return -K(ti);
  }

#endif