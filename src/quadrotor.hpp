#ifndef QUADROTOR_HPP
#define QUADROTOR_HPP
#include<armadillo>
#include"rk4_int.hpp"



class QuadRotor {
    //private class variables like mass, gravity, damping coefficients
	double kt=0.6,km=0.15,arml=0.2,m=0.6,g=9.81;//0.6,0.15,0.2,0.6,9.81
	arma::mat J,Jinv; //arma::vec Jvec = {0.04,0.0375,0.0675};J=arma::diagmat(Jvec)
	arma::vec e3 = {0.,0.,1.};
    public:
        double dt;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr;
        arma::mat xdlist;
        QuadRotor (double _dt);//argument should be any required parameters
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        //inline arma::mat dfdx(const arma::vec& x, const arma::vec& u);
        //inline arma::vec hx(const arma::vec& x);
        void step(void);
        inline arma::vec get_measurement(const arma::vec& x);
        
};

QuadRotor::QuadRotor (double _dt){//add any additional system parameters
    dt = _dt;
	arma::vec Jvec = {0.04,0.0375,0.0675};
	J=arma::diagmat(Jvec);
	Jinv = J.i();
    
}

arma::vec QuadRotor::proj_func (const arma::vec& x){//anglewrapping function if needed
    arma::vec xwrap=x;
    //do nothing if no angle wrapping is required
    return xwrap;
}
inline arma::vec QuadRotor::f(const arma::vec& x, const arma::vec& u){//control affine required
	arma::mat h = arma::reshape(x.subvec(0,15),4,4);
	arma::vec omega = x.subvec(16,18);
	arma::vec v = x.subvec(19,21);
	arma::mat twhat =  {{0.,-omega(2),omega(1),v(0)},
						{-omega(2),0.,-omega(0),v(1)},
						{-omega(1),omega(0),0.,v(2)},
						{0.,0.,0.,0.}};
	arma::mat R = h.submat(0,0,2,2);
	double F = kt*(u(0)+u(1)+u(2)+u(3));
	arma::vec M = {kt*arml*(u(1)-u(3)),
				   kt*arml*(u(2)-u(0)),
				   km*arml*(u(0)-u(1)+u(2)-u(3))
					};
	arma::mat hdot = h*twhat;
	arma::vec omegadot = Jinv*(M + arma::cross(J*omega,omega));
	arma::vec vdot = F/m*e3 - arma::cross(omega,v) - g*R.t()*e3; 
    arma::vec xdot = arma::join_cols(h.as_col(),omega,v);
    return xdot;
}; 
inline arma::vec QuadRotor::get_measurement(const arma::vec& x){
	arma::mat h = arma::reshape(x.subvec(0,15),4,4);
	arma::vec twist = x.subvec(16,21);
	arma::mat R = h.submat(0,0,2,2);
	arma::vec ag = R*(g*e3);
	arma::vec z = arma::join_cols(ag,twist);
	return z;
};
/*
inline arma::mat QuadRotor::dfdx(const arma::vec& x, const arma::vec& u){
    arma::mat A = {//add more rows as needed
        {},
        {},
        {},
        {}
    };
    return A;
}; 

inline arma::vec QuadRotor::hx(const arma::vec& x){
    arma::vec H = {//add more rows as needed
        {},
        {},
        {},
        {}
    };
    return H;
}; 
*/
void QuadRotor::step(){
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};


#endif