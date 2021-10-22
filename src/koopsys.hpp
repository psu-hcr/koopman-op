#ifndef KOOPSYS_HPP
#define KOOPSYS_HPP
#include<armadillo>
#include"rk4_int.hpp"

//const double PI = 3.1415926535987;

template<class basis>
class KoopSys {
    double m, B, g, h;
	basis* zfuncs;
    arma::mat A;
    arma::mat G;
    arma::mat K;
    void calc_K(const arma::vec& x, const arma::vec& u);
    public:
        double dt;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr;
        arma::mat xdlist;
        KoopSys (double _dt, basis *_zfuncs);
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        inline arma::mat dfdx(const arma::vec& x, const arma::vec& u);
        inline arma::vec hx(const arma::vec& x);
        void step(void);
        
        
};

template<class basis>
KoopSys<basis>::KoopSys (double _dt, basis *_zfuncs){
    zfuncs=_zfuncs;//basis functions/functions of the observables
    dt = _dt;//step size
    A = arma::zeros(zfuncs->zdim,zfuncs->zdim);
    G = arma::zeros(zfuncs->zdim,zfuncs->zdim);
}

template<class basis>
arma::vec KoopSys<basis>::proj_func (const arma::vec& x){
    arma::vec xwrap=x;
    xwrap(0) = fmod(x(0)+PI, 2*PI);
    if (xwrap(0) < 0.0) xwrap(0) = xwrap(0) + 2*PI;
    xwrap(0) = xwrap(0) - PI;
    return xwrap;
}
template<class basis>
arma::vec KoopSys<basis>::f(const arma::vec& x, const arma::vec& u){
    arma::vec xdot = {x(1),
                      g/h*sin(x(0)) + B*x(1)/(m*h*h)-u(0)*cos(x(0))/h,
                      x(3),
                      u(0)};;
    return xdot;
}; 
template<class basis>
inline arma::mat KoopSys<basis>::dfdx(const arma::vec& x, const arma::vec& u){
    arma::mat A = {
        {0,1,0,0},
        {g/h*cos(x(0))+ u(0)*sin(x(0))/h, B/(m*h*h),0,0},
        {0,0,0,1},
        {0,0,0,0}
    };
    return A;
}; 
template<class basis>
inline arma::vec KoopSys<basis>::hx(const arma::vec& x){
    arma::vec H = {
        {0},
        {-cos(x(0))/h},
        {0},
        {1}
    };
    return H;
}; 
template<class basis>
void KoopSys<basis>::step(){
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};
template<class basis>
void KoopSys<basis>::calc_K(const arma::vec& x, const arma::vec& u){
    arma::vec ztplus1 = zfuncs->zxu(x,u);
    arma::vec zt = zfuncs->zxu(Xcurr,Ucurr);
    A = A + ztplus1*zt.t();
    G = G + zt*zt.t();
    K=A*arma::pinv(G);

};

#endif