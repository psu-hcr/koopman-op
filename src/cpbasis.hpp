#ifndef CPBASIS_HPP
#define CPBASIS_HPP
#include<armadillo>
#include"rk4_int.hpp"

//const double PI = 3.1415926535987;


class CPBASIS {
	arma::vec zx(const arma::vec& x);
	arma::vec zu(const arma::vec& x, const arma::vec& u);
	public:
        arma::vec zxu(const arma::vec& x, const arma::vec& u);
		int zdim = 10;
		//int xdim = 5;
		//int udim = 5;
        
        
};

//KoopSys::KoopSys (double _dt, basis *_zfuncs){
//    zfuncs=_zfuncs;//basis functions/functions of the observables
//    dt = _dt;//step size
//}

arma::vec CPBASIS::zx (const arma::vec& x){//th,thdot,xc,xcdot,xcdot^2
    arma::vec psix = {x(0),
                      x(1),
					  x(2),
					  x(3),
					  x(3)*x(3)};
	return psix;
}
arma::vec CPBASIS::zu(const arma::vec& x, const arma::vec& u){
    arma::vec psiu = {u(0),
                      u(0)*cos(x(0)),
					  u(0)*cos(x(1)),
                      20.*cos(x(0)*PI/20.)*cos(x(0)*PI/20.),
                      1.};
    return psiu;
}; 

inline arma::vec CPBASIS::zxu(const arma::vec& x,const arma::vec& u){
    arma::vec psixu = arma::join_cols(zx(x),zu(x,u));
    return psixu;
}; 

#endif