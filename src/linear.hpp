#ifndef linear_HPP
#define linear_HPP
#include<armadillo>
#include"rk4_int.hpp"



class linear {
	double k = 1.0; 
	double m = 1.0;
	double b = 0.5;
    public:
        double dt;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr;
        arma::mat xdlist;
        linear (double _dt);//argument should be any required parameters
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        //inline arma::mat dfdx(const arma::vec& x, const arma::vec& u);
        //inline arma::vec hx(const arma::vec& x);
        void step(void);
        inline arma::vec get_measurement(const arma::vec& x);
		
		
		arma::mat A = {{0, 1},{-k/m, -b/m}}; 
 		arma::mat B = {{0, 1},{1/m, 0}}; 
		//arma::mat A = {{0, 1, 0, 0},{0, 0, 0, 0},{0, 0, 0, 1}, {0, 0, 0, 0}}; 
 		//arma::mat B = {{0, 0}, {400, 0}, {0, 0}, {0, 26.67}}; 
		
        
};

linear::linear (double _dt){//add any additional system parameters
    dt = _dt; 
}

arma::vec linear::proj_func (const arma::vec& x){//anglewrapping function if needed
    arma::vec xwrap=x;
    //do nothing if no angle wrapping is required
    return xwrap;
}
inline arma::vec linear::f(const arma::vec& x, const arma::vec& u){//control affine required
    arma::vec xdot = A*x + B*u;
    return xdot;
}; 
inline arma::vec linear::get_measurement(const arma::vec& x){
	return x;
};
/*
inline arma::mat linear::dfdx(const arma::vec& x, const arma::vec& u){
    arma::mat A = {//add more rows as needed
        {},
        {},
        {},
        {}
    };
    return A;
}; 

inline arma::vec linear::hx(const arma::vec& x){
    arma::vec H = {//add more rows as needed
        {},
        {},
        {},
        {}
    };
    return H;
}; 
*/
void linear::step(){
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};


#endif