#ifndef twomass_HPP
#define twomass_HPP
#include<armadillo>
#include"rk4_int.hpp"



class twomass {
	double k = 3.0; 
	double m = 1.5;
    public:
        double dt;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr;
        arma::mat xdlist;
        twomass (double _dt);//argument should be any required parameters
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        //inline arma::mat dfdx(const arma::vec& x, const arma::vec& u);
        //inline arma::vec hx(const arma::vec& x);
        void step(void);
        inline arma::vec get_measurement(const arma::vec& x);
		
		arma::mat A = {{0, 1, 0, 0},{-(k+k)/m, 0, k/m, 0},{0, 0, 0, 1}, {k/m, 0, -k/m, 0}}; 
 		arma::mat B = {{0, 0}, {1/m, 0}, {0, 0}, {0, 1/m}}; 
		
        
};

twomass::twomass (double _dt){//add any additional system parameters
    dt = _dt; 
}

arma::vec twomass::proj_func (const arma::vec& x){//anglewrapping function if needed
    arma::vec xwrap=x;
    //do nothing if no angle wrapping is required
    return xwrap;
}
inline arma::vec twomass::f(const arma::vec& x, const arma::vec& u){//control affine required
    arma::vec xdot = A*x + B*u;
    return xdot;
}; 
inline arma::vec twomass::get_measurement(const arma::vec& x){
	return x;
};
/*
inline arma::mat twomass::dfdx(const arma::vec& x, const arma::vec& u){
    arma::mat A = {//add more rows as needed
        {},
        {},
        {},
        {}
    };
    return A;
}; 

inline arma::vec twomass::hx(const arma::vec& x){
    arma::vec H = {//add more rows as needed
        {},
        {},
        {},
        {}
    };
    return H;
}; 
*/
void twomass::step(){
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};


#endif