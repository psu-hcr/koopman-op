#ifndef QUADBASIS_HPP
#define QUADBASIS_HPP
#include<armadillo>



class QuadBasis {
	
	public:
		arma::vec proj_func (const arma::vec& x);
        arma::vec zxu(const arma::vec& x, const arma::vec& u);
		arma::vec zx(const arma::vec& x);
		arma::vec zu(const arma::vec& x, const arma::vec& u);
		arma::mat dvdu(const arma::vec& z);
		int zdim = 22;//set to length of z including v(x,u)
		int xdim = 18;//set to length of z(x(t)) only        
        
};

arma::vec QuadBasis::proj_func (const arma::vec& x){//angle wrapping function
    arma::vec xwrap=x;
    //do nothing if not angle wrapping is required
    return xwrap;
}

arma::vec QuadBasis::zx (const arma::vec& x){//x is the measurement ag,omega,v
    arma::vec psix = {x(8)*x(4),
					  x(7)*x(5),
					  x(8)*x(3),
					  x(6)*x(5),
					  x(7)*x(3),
					  x(6)*x(4),
					  x(4)*x(5),
					  x(3)*x(5),
					  x(3)*x(4)};
	 psix = arma::join_cols(x.subvec(0,8),psix);
	return psix;
}
arma::vec QuadBasis::zu(const arma::vec& x, const arma::vec& u){//recommend psiu=u;
    arma::vec psiu = u;
    return psiu;cout<<arma::size(psiu)<<endl;
}; 

arma::mat QuadBasis::dvdu(const arma::vec& z){//if psiu=u this is the identity
    arma::mat psiu = {{1.,0.,0.,0.},
					  {0.,1.,0.,0.},
					  {0.,0.,1.,0.},
					  {0.,0.,0.,1.}};
    return psiu;
};


inline arma::vec QuadBasis::zxu(const arma::vec& x,const arma::vec& u){
    arma::vec psixu = arma::join_cols(zx(x),zu(x,u));
    return psixu;
}; 

#endif