#ifndef twomassBasis_HPP
#define twomassBasis_HPP
#include<armadillo>



class twomassBasis {
	
	public:
		arma::vec proj_func (const arma::vec& x);
        arma::vec zxu(const arma::vec& x, const arma::vec& u);
		arma::vec zx(const arma::vec& x);
		arma::vec zu(const arma::vec& x, const arma::vec& u);
		arma::mat dvdu(const arma::vec& z);
		int zdim = 6;//set to length of z including v(x,u)
		int xdim = 4;//set to length of z(x(t)) only
        
};

arma::vec twomassBasis::proj_func (const arma::vec& x){//angle wrapping function
    arma::vec xwrap=x;
    //do nothing if not angle wrapping is required
    return xwrap;
}

arma::vec twomassBasis::zx (const arma::vec& x){//x is the measurement ag,omega,v
	arma::vec psix = x;
	return psix;
}
arma::vec twomassBasis::zu(const arma::vec& x, const arma::vec& u){//recommend psiu=u;
    arma::vec psiu = u;
    return psiu;
}; 

arma::mat twomassBasis::dvdu(const arma::vec& z){//if psiu=u this is the identity
    arma::mat psiu = {{1, 0},{0, 1}};
    return psiu;
};


inline arma::vec twomassBasis::zxu(const arma::vec& x,const arma::vec& u){
    arma::vec psixu = arma::join_cols(zx(x),zu(x,u));
    return psixu;
}; 

#endif