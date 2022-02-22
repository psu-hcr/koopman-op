#ifndef ERRORCOSTIK_HPP
#define ERRORCOSTIK_HPP
#include<armadillo>
#include <math.h>

template <class system, class trajectory>
class errorcostIK {
    system* sys;
    double eps = pow(10,-5);
    double infw = 100.;
    
  public:
    arma::mat Q;
    arma::mat R;
    arma::mat Sigma;
    arma::mat B;
	arma::vec BoundaryDiag;
	arma::vec Joint_limit;
	double QR;
	double fisher;
	double boundary;
	int NumOfJoint;
	trajectory* traj;
    
    // initialization
    errorcostIK(arma::mat _Q, arma::mat _R,trajectory* _traj, system *_sys, const arma::vec SigDiag, const arma::vec _JoLimit){
      // Q: gain for x; R: gain for u; 
	  // traj:desire trajectory clas; sys: system for simulation;
      // SigDiag: Diag for fisher cost; 
	  // JoLimit: vector of Joint limit;
      Q=_Q; R=_R; sys=_sys; traj=_traj; Sigma = arma::diagmat(SigDiag);
	  Joint_limit = _JoLimit;
	  BoundaryDiag = arma::zeros<arma::vec>(sys->zfuncs->xdim);
	  NumOfJoint = Joint_limit.n_cols;
	  QR = 0.0;
	  fisher = 0.0;
	  boundary = 0.0;
      //std::cout<<B<<std::endl;
      };
      
    inline double l (const arma::vec& x,const arma::vec& u,double ti){
      return QR_cost(x,u,ti) + fisher_cost(x,u,ti) + boundary_cost(x,u,ti);
      }
	  
	 inline double QR_cost (const arma::vec& x,const arma::vec& u,double ti){
	 	//arma::vec xproj = sys->proj_func(x);
	 	arma::vec xd = traj->desire_jointstate(ti);
		//return arma::as_scalar(((xproj.t()-xd.t())*Q*(xproj-xd)+u.t()*R*u)/2);
		return arma::as_scalar(((x.t()-xd.t())*Q*(x-xd)+u.t()*R*u)/2);
	 }
	 
	 inline double fisher_cost (const arma::vec& x,const arma::vec& u,double ti){
	 	arma::vec dfdk = sys->zfuncs->zx(x);
      	double fishTopt = arma::as_scalar(dfdk.t()*Sigma*dfdk);
	 	return (infw*1.)/(fishTopt+eps);
	 }
	 
	 inline double boundary_cost (const arma::vec& x,const arma::vec& u,double ti){
	 	//arma::vec xproj = sys->proj_func(x);
		for (int i=0; i<NumOfJoint; i++){
			//BoundaryDiag(i) = pow((xproj(i)/(2*Joint_limit(i))), 8);
			BoundaryDiag(i) = pow((x(i)/(2*Joint_limit(i))), 8);
		}
		//std::cout<<BoundaryDiag<<std::endl;
		B = arma::diagmat(BoundaryDiag);
		return arma::as_scalar(x.t()*B*x);
	}
      
    inline arma::vec dldx (const arma::vec& x,const arma::vec& u,double ti){
      //arma::vec xproj = sys->proj_func(x);
      arma::vec xd = traj->desire_jointstate(ti);
      arma::vec dfdk = sys->zfuncs->zx(x);
      double fishTopt = arma::as_scalar(dfdk.t()*Sigma*dfdk);
	  for (int i=0; i<NumOfJoint; i++){
			//BoundaryDiag(i) = pow((xproj(i)/(2*Joint_limit(i))), 8);
			BoundaryDiag(i) = pow((x(i)/(2*Joint_limit(i))), 8);
	  }
	  B = arma::diagmat(BoundaryDiag);
      //return Q*(xproj-xd) + ((-2.*infw)/pow(fishTopt+eps,2))*Sigma*dfdk + 10*B*xproj;
	  return Q*(x-xd) + ((-2.*infw)/pow(fishTopt+eps,2))*Sigma*dfdk + 10*B*x;
      }
      
    double calc_cost(const arma::mat& x,const arma::mat& u){
      //arma::vec xproj;
      double J1 = 0.0;
	  QR = 0.0;
	  fisher = 0.0;
	  boundary = 0.0;
      for (int i = 0; i<x.n_cols; i++){
        //xproj = sys->proj_func(x.col(i));
		//std::cout<<i<<std::endl;
		QR+=QR_cost(x.col(i),u.col(i),sys->tcurr+(double)i*sys->dt);
		fisher+=fisher_cost(x.col(i),u.col(i),sys->tcurr+(double)i*sys->dt);
		boundary+=boundary_cost(x.col(i),u.col(i),sys->tcurr+(double)i*sys->dt);
        J1+=l(x.col(i),u.col(i),sys->tcurr+(double)i*sys->dt);
      }
      return J1;
      }
	 
    
};

#endif