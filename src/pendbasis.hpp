#ifndef pendBasis_HPP
#define pendBasis_HPP
#include<armadillo>



class pendBasis {
	
	public:
		arma::vec proj_func (const arma::vec& x);
        arma::vec zxu(const arma::vec& x, const arma::vec& u);
		arma::vec zx(const arma::vec& x);
		arma::vec zu(const arma::vec& x, const arma::vec& u);
		arma::mat dvdu(const arma::vec& z);
		int zdim = 5;//set to length of z including v(x,u)
		int xdim = 4;//set to length of z(x(t)) only
        
};

arma::vec pendBasis::proj_func (const arma::vec& x){//angle wrapping function
    arma::vec xwrap=x;
    xwrap(0) = fmod(x(0)+PI, 2*PI);
    if (xwrap(0) < 0.0) xwrap(0) = xwrap(0) + 2*PI;
    xwrap(0) = xwrap(0) - PI;
    return xwrap;
}

arma::vec pendBasis::zx (const arma::vec& x){//x is the measurement ag,omega,v
	arma::vec psix = { x(0),
									x(1),
									sin(x(0)),
									cos(x(0))
									};
	return psix;
}
arma::vec pendBasis::zu(const arma::vec& x, const arma::vec& u){//recommend psiu=u;
    arma::vec psiu = u;
    return psiu;
}; 

arma::mat pendBasis::dvdu(const arma::vec& z){//if psiu=u this is the identity
    arma::mat psiu = {1};
    return psiu;
};


inline arma::vec pendBasis::zxu(const arma::vec& x,const arma::vec& u){
    arma::vec psixu = arma::join_cols(zx(x),zu(x,u));
    return psixu;
}; 

#endif