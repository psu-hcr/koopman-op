#ifndef twolink_HPP
#define twolink_HPP
#include<armadillo>
#include"rk4_int.hpp"

const double PI = 3.1415926535987;

class twolink {
    // private class variables
	
	// mass, length and Radius
	double m1 = 1, m2 = 1, L1 = 0.5, L2 = 0.3, R1 = 0.1, R2 = 0.1;
	
	// simplified coefficient
	double a = 0.166667*m2*pow(L2,2) + 0.5*m2*pow(R2,2);
	double b = 0.25*m1*pow(R1,2);
	double c = 0.41667*m2*pow(L2,2);
	
    public:
		// public variable
        double dt;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr;
        twolink (double _dt);//argument should be any required parameters
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        inline arma::mat dfdx(const arma::vec& x);
        inline arma::mat hx(const arma::vec& x);
        void step(void);
		inline arma::vec get_measurement(const arma::vec& x);
};

twolink::twolink (double _dt){
	//initialize twolink class 
    dt = _dt;	// time step
}

arma::vec twolink::proj_func (const arma::vec& x){
	// anglewrapping function if needed
    arma::vec xwrap=x;
    xwrap(0) = fmod(x(0)+PI, 2*PI);
    if (xwrap(0) < 0.0) xwrap(0) = xwrap(0) + 2*PI;
    xwrap(0) = xwrap(0) - PI;
	xwrap(2) = fmod(x(2)+PI, 2*PI);
    if (xwrap(2) < 0.0) xwrap(2) = xwrap(2) + 2*PI;
    xwrap(2) = xwrap(2) - PI;
    return xwrap;
}

inline arma::vec twolink::f(const arma::vec& x, const arma::vec& u){
    arma::vec xdot = {
		x[1],
		//-(a*sin(2*x[2])/(a*pow(sin(x[2]),2)+b))*x[1] + 1/(a*pow(sin(x[2]), 2)+b)*u[0],
		-(a*sin(2*x[2])/(a*pow(sin(x[2]),2)+b))*x[1] + u[0],
		x[3],
		-a*sin(2*x[2])/c*x[1] + 1/c*u[1]
	};
    return xdot;
}; 

inline arma::vec twolink::get_measurement(const arma::vec& x){
	return proj_func(x);
};

inline arma::mat twolink::dfdx(const arma::vec& x){
	double ax = -a*sin(2*x[2])/(a*pow(sin(x[2]),2)+b);
	double bx = -x[1]*(2*a*cos(2*x[2])*(a*pow(sin(x[2]),2)+1) - a*sin(2*x[2])*(2*a*sin(x[2])*cos(x[2])))/
		pow((a*pow(sin(x[2]),2) + b),2);
	double cx = -a*sin(2*x[2])/c;
	double dx = -x[1]*2*a*cos(2*x[2])/c;
    arma::mat A = {
        {0, 1, 0, 0},
        {0, ax, bx, 0},
        {0, 0, 0, 1},
        {0, cx, dx, 0}
    };
    return A;
}; 

inline arma::mat twolink::hx(const arma::vec& x){
    arma::mat B = {//add more rows as needed
        {0, 0},
        //{1/(a*pow(sin(x[2]),2) + b), 0},
		{1, 0},
        {0, 0},
        {0, 1/c}
    };
    return B;
}; 


void twolink::step(){
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};


#endif