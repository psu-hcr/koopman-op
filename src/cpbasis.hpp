#ifndef CPBASIS_HPP
#define CPBASIS_HPP
#include<armadillo>
#include"rk4_int.hpp"

//const double PI = 3.1415926535987;


class CPBASIS {
	
	public:
        arma::vec zxu(const arma::vec& x, const arma::vec& u);
		arma::vec zx(const arma::vec& x);
		arma::vec zu(const arma::vec& x, const arma::vec& u);
		arma::vec dvdu(const arma::vec& z);
		int zdim = 7;
		int xdim = 6;
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
					  x(3)*x(3),
					 1.0};
	return psix;
}
arma::vec CPBASIS::zu(const arma::vec& x, const arma::vec& u){
    arma::vec psiu = {u(0)};/*,	make it just the control so Ku = B								
                      u(0)*cos(x(0)),						
					  u(0)*cos(x(1)),
                      20.*cos(x(0)*PI/20.)*cos(x(0)*PI/20.)};*/
    return psiu;
}; 
/*
arma::vec CPBASIS::dvdu(const arma::vec& z){
    arma::vec psiu = {{1.},
					  {cos(z(0))},
					  {cos(z(1))},
					  {(-40.*PI/20.)*sin(z(0)*PI/20.)}};
    return psiu;
};
*/

inline arma::vec CPBASIS::zxu(const arma::vec& x,const arma::vec& u){
    arma::vec psixu = arma::join_cols(zx(x),zu(x,u));
    return psixu;
}; 

#endif