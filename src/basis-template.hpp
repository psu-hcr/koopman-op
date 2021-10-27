#ifndef BASISName_HPP
#define BASISName_HPP
#include<armadillo>


class BasisName {
	
	public:
		arma::vec proj_func (const arma::vec& x);
        arma::vec zxu(const arma::vec& x, const arma::vec& u);
		arma::vec zx(const arma::vec& x);
		arma::vec zu(const arma::vec& x, const arma::vec& u);
		arma::vec dvdu(const arma::vec& z);
		int zdim = ;//set to length of z including v(x,u)
		int xdim = ;//set to length of z(x(t)) only        
        
};

arma::vec BasisName::proj_func (const arma::vec& x){//angle wrapping function
    arma::vec xwrap=x;
    //do nothing if not angle wrapping is required
    return xwrap;
}

arma::vec BasisName::zx (const arma::vec& x){//note the order of x with a comment here
    arma::vec psix = {};
	return psix;
}
arma::vec BasisName::zu(const arma::vec& x, const arma::vec& u){//recommend psiu=u;
    arma::vec psiu = {};
    return psiu;
}; 

arma::vec BasisName::dvdu(const arma::vec& z){//if psiu=u these are all 1s
    arma::vec psiu = {{1.},
					  {1.},
					  {1.},
					  {1.}};
    return psiu;
};


inline arma::vec BasisName::zxu(const arma::vec& x,const arma::vec& u){
    arma::vec psixu = arma::join_cols(zx(x),zu(x,u));
    return psixu;
}; 

#endif