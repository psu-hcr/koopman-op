#ifndef linearBasis_HPP
#define linearBasis_HPP
#include<armadillo>



class linearBasis {
	
	public:
		arma::vec proj_func (const arma::vec& x);
        arma::vec zxu(const arma::vec& x, const arma::vec& u);
		arma::vec zx(const arma::vec& x);
		arma::vec zu(const arma::vec& x, const arma::vec& u);
		arma::mat dvdu(const arma::vec& z);
		//int zdim = 6;//set to length of z including v(x,u)
		//int xdim = 4;//set to length of z(x(t)) only
		int zdim = 4;
		int xdim = 2;
        
};

arma::vec linearBasis::proj_func (const arma::vec& x){//angle wrapping function
    arma::vec xwrap=x;
    //do nothing if not angle wrapping is required
    return xwrap;
}

arma::vec linearBasis::zx (const arma::vec& x){//x is the measurement ag,omega,v
	arma::vec psix = x;
	return psix;
}
arma::vec linearBasis::zu(const arma::vec& x, const arma::vec& u){//recommend psiu=u;
    arma::vec psiu = u;
    return psiu;cout<<arma::size(psiu)<<endl;
}; 

arma::mat linearBasis::dvdu(const arma::vec& z){//if psiu=u this is the identity
    arma::mat psiu = {{1, 0},{0, 1}};
    return psiu;
};


inline arma::vec linearBasis::zxu(const arma::vec& x,const arma::vec& u){
    arma::vec psixu = arma::join_cols(zx(x),zu(x,u));
    return psixu;
}; 

#endif