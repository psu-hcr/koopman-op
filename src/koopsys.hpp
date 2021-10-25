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
    double Mindex=1;
    arma::mat Kdisc;
    arma::mat Kx,Ku;
    public:
        arma::mat K;
        double dt;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr,Xprev,Uprev,Zcurr;
        arma::mat xdlist;
        KoopSys (double _dt, basis *_zfuncs);
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        inline arma::mat dfdx(const arma::vec& x, const arma::vec& u);
        inline arma::vec hx(const arma::vec& x);
        void step(void);
        void update_XU(const arma::vec& x,const arma::vec& u);
        void calc_K(void);
        
};

template<class basis>
KoopSys<basis>::KoopSys (double _dt, basis *_zfuncs){
    zfuncs=_zfuncs;//basis functions/functions of the observables
    dt = _dt;//step size
    A = arma::zeros(zfuncs->zdim,zfuncs->zdim);
    G = arma::zeros(zfuncs->zdim,zfuncs->zdim);
    K = arma::ones<arma::mat>(zfuncs->zdim,zfuncs->zdim);
    //Kx = arma::eye(zfuncs->xdim,zfuncs->xdim);
    //Ku = arma::eye(zfuncs->udim,zfuncs->udim);
    
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
arma::vec KoopSys<basis>::f(const arma::vec& z, const arma::vec& u){
    arma::vec zx = z.subvec(0,zfuncs->xdim-1);
    arma::vec zu = z.subvec(zfuncs->xdim,zfuncs->zdim-1);
    arma::vec zdot = Kx*zx+Ku*zu;
    return zdot;
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
    Zcurr = zfuncs->zxu(Xcurr,Ucurr);
    Zcurr = RK4_step(this,Zcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};

template<class basis>
void KoopSys<basis>::update_XU(const arma::vec& x,const arma::vec& u){
    Uprev = Ucurr;
    Ucurr = u;
    Xprev = Xcurr;
    Xcurr = x;
    Mindex++;
};
template<class basis>
void KoopSys<basis>::calc_K(){
    arma::vec ztplus1 = zfuncs->zxu(Xcurr,Ucurr);
    arma::vec zt = zfuncs->zxu(Xprev,Uprev);
    A = ((Mindex-1)*A + ztplus1*zt.t())/Mindex;
    G = ((Mindex-1)*G + zt*zt.t())/Mindex;
    Kdisc=A*arma::pinv(G);
    arma::cx_mat Ktemp;
    try{
    Ktemp=arma::logmat(Kdisc);
    K=arma::real(Ktemp);
    K=K/dt;//cout<<"NO Error here!"<<endl;
    }
    catch (...){
    //cout<<"Error here!"<<endl;
    K = arma::ones<arma::mat>(zfuncs->zdim,zfuncs->zdim);
    //cout<<"Error"<<arma::as_scalar(K(0,0))<<endl;
    }
    Kx = K.submat(0,0,K.n_cols-1, zfuncs->xdim-1);
    Ku = K.submat(0,zfuncs->xdim,K.n_cols-1,K.n_cols-1);
};

#endif