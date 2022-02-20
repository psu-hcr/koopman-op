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
    double Mindex=0;
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
        void calc_K(const arma::vec& x,const arma::vec& u);
        
};

template<class basis>
KoopSys<basis>::KoopSys (double _dt, basis *_zfuncs){
    zfuncs=_zfuncs;//basis functions/functions of the observables
    dt = _dt;//step size
    A = arma::zeros(zfuncs->zdim,zfuncs->zdim);
    G = arma::zeros(zfuncs->zdim,zfuncs->zdim);
	K = 0.1*arma::ones<arma::mat>(zfuncs->zdim,zfuncs->zdim);
    //K = arma::randn<arma::mat>(zfuncs->zdim,zfuncs->zdim);
    //Kx = arma::eye(zfuncs->xdim,zfuncs->xdim);
    //Ku = arma::eye(zfuncs->udim,zfuncs->udim);
    
}

template<class basis>
arma::vec KoopSys<basis>::proj_func (const arma::vec& x){
   return zfuncs->proj_func(x);
}
template<class basis>
arma::vec KoopSys<basis>::f(const arma::vec& zx, const arma::vec& u){
    arma::vec zu = zfuncs->zu(zx,u);
	arma::vec zdot = Kx*zx+Ku*zu;
    return zdot;
}; 
template<class basis>
inline arma::mat KoopSys<basis>::dfdx(const arma::vec& x, const arma::vec& u){
    arma::mat Ax = Kx;
    return Ax;
}; 
template<class basis>
inline arma::mat KoopSys<basis>::hx(const arma::vec& z){
    arma::mat Bx = Ku*zfuncs->dvdu(z);//or just Ku?
    
    return Bx;
}; 
template<class basis>
void KoopSys<basis>::step(){ 
    Xcurr = Xcurr+f(Xcurr,Ucurr)*dt;//RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};

template<class basis>
void KoopSys<basis>::update_XU(const arma::vec& x,const arma::vec& u){
    Uprev = Ucurr;
    Ucurr = u;
    Xprev = Xcurr;
    Xcurr = zfuncs->zx(x);
	tcurr = tcurr+dt;
    
};
template<class basis>
void KoopSys<basis>::calc_K(const arma::vec& x,const arma::vec& u){
	Mindex++;
	update_XU(x,u);
    arma::vec ztplus1 = zfuncs->zxu(Xcurr,Uprev);
    arma::vec zt = zfuncs->zxu(Xprev,Uprev);//cout<<ztplus1.t()<<endl;
    A = A + (ztplus1*zt.t()-A)/Mindex;//cout<<A<<endl;
    G = G + (zt*zt.t()-G)/Mindex;//cout<<G<<endl;
	try{
    Kdisc=A*arma::pinv(G);
	K = 0.1*arma::ones<arma::mat>(zfuncs->zdim,zfuncs->zdim);
	//K = arma::randn<arma::mat>(zfuncs->zdim,zfuncs->zdim);
	arma::cx_mat Ktemp;
    Ktemp=arma::logmat(Kdisc);//dt;
    K=arma::real(Ktemp);//cout<<K<<endl;
    }
    catch (...){//cout<<"This is a problem."<<endl;
    K = 0.1*arma::ones<arma::mat>(zfuncs->zdim,zfuncs->zdim);
	//K = arma::randn<arma::mat>(zfuncs->zdim,zfuncs->zdim);
    }
    Kx = K.submat(0,0,zfuncs->xdim-1, zfuncs->xdim-1);//cout<<Kx<<endl;
    Ku = K.submat(0,zfuncs->xdim,zfuncs->xdim-1,K.n_cols-1);
};

#endif