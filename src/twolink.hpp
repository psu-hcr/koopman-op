#ifndef TWOLINK_HPP
#define TWOLINK_HPP
#include<armadillo>
#include"rk4_int.hpp"



class twoLink {
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
        twoLink (double _dt, arma::vec _Xinit);//argument should be any required parameters
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        inline arma::mat dfdx(const arma::vec& x);
        inline arma::mat hx(const arma::vec& x);
        void step(void);
};

twoLink::twoLink (double _dt, arma::vec _Xinit){
	//initialize twoLink class 
    dt = _dt;	// time step
	Xcurr = _Xinit;	// initial condition, it should be a 4 elements vector
	//std::cout<<Xcurr<<std::endl;
	Ucurr = {0, 0};	
	//std::cout<<a<<"\n"<<b<<"\n"<<c<<"\n"<<std::endl;
}

arma::vec twoLink::proj_func (const arma::vec& x){
	// anglewrapping function if needed
    arma::vec xwrap=x; // do nothing if no angle wrapping is required
    return xwrap;
}

inline arma::vec twoLink::f(const arma::vec& x, const arma::vec& u){
	// differentail equation for system. control affine required
	// x is a 4 elements vector, u is a 2 elements vector
    arma::vec xdot = {x[1],
					  -(a*sin(2*x[3])/(a*pow(sin(x[3]),2)+b))*x[2] + 1/(a*pow(sin(x[3]), 2)+b)*u[0],
					  x[3],
					  -a*sin(2*x[3])/c*x[2] + 1/c*u[1]
	};
	//std::cout<<xdot<<std::endl;
    return xdot;
}; 


inline arma::mat twoLink::dfdx(const arma::vec& x){
	double a = 8*sin(2*x[2])/(8*pow(sin(x[2]),2)+1);
	double b = x[1]*(16*cos(2*x[2])*(8*pow(sin(x[2]),2)+1) - 8*sin(2*x[2])*(16*sin(x[2])*cos(x[2])))/
		pow((8*pow(sin(x[2]),2) + 1),2);
	double c = 0.533*sin(2*x[2]);
	double d = x[1]*1.066*cos(2*x[2]);
    arma::mat A = {//add more rows as needed
        {0, 1, 0, 0},
        {0, a, b, 0},
        {0, 0, 0, 1},
        {0, c, d, 0}
    };
    return A;
}; 

inline arma::mat twoLink::hx(const arma::vec& x){
    arma::mat B = {//add more rows as needed
        {0, 0},
        {20000/(800*pow(sin(x[2]),2) + 50), 0},
        {0, 0},
        {0, 26.67}
    };
    return B;
}; 


void twoLink::step(){
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};


#endif