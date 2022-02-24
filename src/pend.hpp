#ifndef pend_HPP
#define pend_HPP
#include<armadillo>
#include"rk4_int.hpp"

const double PI = 3.1415926535987;

class pend {
	double L = 3.0; 
	double m = 1.5;
	double g = 9.81;
    public:
        double dt;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr;
        arma::mat xdlist;
        pend (double _dt);//argument should be any required parameters
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        inline arma::mat dfdx(const arma::vec& x, const arma::vec& u);
        inline arma::vec hx(const arma::vec& x);
        void step(void);
        inline arma::vec get_measurement(const arma::vec& x);
		arma::vec B = {0, -1/(m*L) };
        
};

pend::pend (double _dt){//add any additional system parameters
    dt = _dt; 
}

arma::vec pend::proj_func (const arma::vec& x){//anglewrapping function if needed
    arma::vec xwrap=x;
    xwrap(0) = fmod(x(0)+PI, 2*PI);
    if (xwrap(0) < 0.0) xwrap(0) = xwrap(0) + 2*PI;
    xwrap(0) = xwrap(0) - PI;
    return xwrap;
}
inline arma::vec pend::f(const arma::vec& x, const arma::vec& u){//control affine required
    arma::vec xdot ={x(1), 
								   -g/L*sin(x(0)) +1/(m*L)*u(0) };
    return xdot;
}; 
inline arma::vec pend::get_measurement(const arma::vec& x){
	return x;
};

inline arma::mat pend::dfdx(const arma::vec& x, const arma::vec& u){
    arma::mat A = {//add more rows as needed
        {0, 1},
        {-g/L*cos(x(0)), 0}
    };
    return A;
}; 

inline arma::vec pend::hx(const arma::vec& x){
    arma::vec H = {//add more rows as needed
        {0},
        {-1/(m*L)}
    };
    return H;
}; 

void pend::step(){
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};


#endif