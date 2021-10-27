#ifndef systname_HPP
#define systname_HPP
#include<armadillo>
#include"rk4_int.hpp"



class SystName {
    //private class variables like mass, gravity, damping coefficients
    public:
        double dt;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr;
        arma::mat xdlist;
        SystName (double _dt);//argument should be any required parameters
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        inline arma::mat dfdx(const arma::vec& x, const arma::vec& u);//optional if model not used for control
        inline arma::vec hx(const arma::vec& x);//optional if model not used for control
        void step(void);
        
        
};

SystName::SystName (double _dt){//add any additional system parameters
    dt = _dt;//save system parameters to private variables
    
}

arma::vec SystName::proj_func (const arma::vec& x){//anglewrapping function if needed
    arma::vec xwrap=x;
    //do nothing if no angle wrapping is required
    return xwrap;
}
inline arma::vec SystName::f(const arma::vec& x, const arma::vec& u){//control affine required
    arma::vec xdot = {//makes a column vector
                      };
    return xdot;
}; 

inline arma::mat SystName::dfdx(const arma::vec& x, const arma::vec& u){
    arma::mat A = {//add more rows as needed
        {},
        {},
        {},
        {}
    };
    return A;
}; 

inline arma::vec SystName::hx(const arma::vec& x){
    arma::vec H = {//add more rows as needed
        {},
        {},
        {},
        {}
    };
    return H;
}; 

void SystName::step(){
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};


#endif