#ifndef KOOPSYS_HPP
#define KOOPSYS_HPP
#include<armadillo>
#include"rk4_int.hpp"

//const double PI = 3.1415926535987;

template<class basis>
class KoopSys {
    double m, B, g, h;
	
    arma::mat A;
    arma::mat G;
    double Mindex=1;
    arma::mat Kdisc;
    //arma::mat Kx,Ku;
    public:
		basis* zfuncs;
        arma::mat K,Kx,Ku;
        double dt;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr,Xprev,Uprev,Zcurr;
        arma::mat xdlist;
        KoopSys (double _dt, basis *_zfuncs);
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        inline arma::mat dfdx(const arma::vec& x, const arma::vec& u);
        inline arma::mat hx(const arma::vec& x);
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
    K = arma::randn<arma::mat>(zfuncs->zdim,zfuncs->zdim);
    //Kx = arma::eye(zfuncs->xdim,zfuncs->xdim);
    //Ku = arma::eye(zfuncs->udim,zfuncs->udim);
    
}

template<class basis>
arma::vec KoopSys<basis>::proj_func (const arma::vec& x){
   return zfuncs->proj_func(x);
}
template<class basis>
arma::vec KoopSys<basis>::f(const arma::vec& zx, const arma::vec& u){
    //arma::vec zx = z.subvec(0,zfuncs->xdim-1);
    arma::vec zu = zfuncs->zu(zx,u);
    arma::vec zdot = Kx*zx+Ku*zu;
    return zdot;
}; 
template<class basis>
inline arma::mat KoopSys<basis>::dfdx(const arma::vec& x, const arma::vec& u){
    arma::mat A = Kx;
    return A;
}; 
template<class basis>
inline arma::mat KoopSys<basis>::hx(const arma::vec& z){
    arma::mat B = Ku*zfuncs->dvdu(z);//or just Ku?
    
    return B;
}; 
template<class basis>
void KoopSys<basis>::step(){ 
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};

template<class basis>
void KoopSys<basis>::update_XU(const arma::vec& x,const arma::vec& u){
    Uprev = Ucurr;
    Ucurr = u;
    Xprev = Xcurr;
    Xcurr = zfuncs->zx(x);
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
    }
    catch (...){
    K = 0.1*arma::ones<arma::mat>(zfuncs->zdim,zfuncs->zdim);
    }
    Kx = K.submat(0,0,zfuncs->xdim-1, zfuncs->xdim-1);
    Ku = K.submat(0,zfuncs->xdim,zfuncs->xdim-1,K.n_cols-1);
};

#endif