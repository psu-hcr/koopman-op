#ifndef LINEARTL_HPP
#define LINEARTL_HPP
#include<armadillo>
#include"rk4_int.hpp"



class LinearTL {
	
    public:
		// public variable
        double dt;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr;
        LinearTL (double _dt, arma::vec _Xinit);//argument should be any required parameters
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        inline arma::mat dfdx(const arma::vec& x);
        inline arma::mat hx(const arma::vec& x, const arma::vec& u);
        void step(void);
		arma::mat A = {{0, 1, 0, 0},{0, 0, 0, 0},{0, 0, 0, 1}, {0, 0, 0, 0}}; 
 		arma::mat B = {{0, 0}, {400, 0}, {0, 0}, {0, 26.67}}; 
};

LinearTL::LinearTL (double _dt, arma::vec _Xinit){
	//initialize LinearTL class 
    dt = _dt;	// time step
	Xcurr = _Xinit;	// initial condition, it should be a 4 elements vector
	//std::cout<<Xcurr<<std::endl;
	Ucurr = {0, 0};	
	//std::cout<<a<<"\n"<<b<<"\n"<<c<<"\n"<<std::endl;
}

arma::vec LinearTL::proj_func (const arma::vec& x){
	// anglewrapping function if needed
    arma::vec xwrap=x; // do nothing if no angle wrapping is required
    return xwrap;
}

inline arma::vec LinearTL::f(const arma::vec& x, const arma::vec& u){
	// differentail equation for system. control affine required
	// x is a 4 elements vector, u is a 2 elements vector
    arma::vec xdot = A*x+B*u;
	//std::cout<<xdot<<std::endl;
    return xdot;
}; 


inline arma::mat LinearTL::dfdx(const arma::vec& x){
    return A;
}; 

inline arma::mat LinearTL::hx(const arma::vec& x, const arma::vec& u){
    return B;
}; 


void LinearTL::step(){
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};


#endif