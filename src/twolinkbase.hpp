#ifndef twolinkbase_HPP
#define twolinkbase_HPP
#include<armadillo>
#include<cmath>


class twolinkbase {
	
	public:
		arma::vec proj_func (const arma::vec& x);
        arma::vec zxu(const arma::vec& x, const arma::vec& u);
		arma::vec zx(const arma::vec& x);
		arma::vec zu(const arma::vec& x, const arma::vec& u);
		arma::mat dvdu(const arma::vec& z);
		double angle_wrap(double input);
		int zdim = 6;// set to length of z including v(x,u)
		int xdim = 4;// set to length of z(x(t)) only 
		int udim = zdim-xdim;
        
};

double twolinkbase::angle_wrap (double input){
	// angle wrapping function
    double output = fmod(input+M_PI, 2*M_PI);
    if (output < 0.0) output = output + 2*M_PI;
    output = output - M_PI;
    return output;
}

arma::vec twolinkbase::proj_func (const arma::vec& x){
	// project function
    arma::vec xwrap=x;
    return xwrap;
}

arma::vec twolinkbase::zx (const arma::vec& x){// note the order of x with a comment here
	arma::vec pos = {x(0),
					 x(1)};
	arma::vec vel = {x(2),
					 x(3)};
    /*arma::vec psix = {pos(0),
					  pos(1),
					  vel(0),
					  vel(1),
					  pos(0)*pos(1),
					  pow(pos(0), 2)*pow(pos(1), 2),
					  pow(pos(0), 3)*pow(pos(1), 3),
					  vel(0)*vel(1),
					  pow(vel(0), 2)*pow(vel(1), 2),
					  pow(vel(0), 3)*pow(vel(1), 3),
					  1,
					  8*sin(2*pos(1))/(8*pow(sin(pos(1)), 2)+1)*vel(0),
					  sin(2*pos(1))*vel(0)
					  };
	return psix;*/
	return x;
}
arma::vec twolinkbase::zu(const arma::vec& x, const arma::vec& u){//recommend psiu=u;
    arma::vec psiu = {u};
    return psiu;
}; 

arma::mat twolinkbase::dvdu(const arma::vec& z){//if psiu=u this is the identity
     arma::mat dvu= 1.*arma::eye(udim,udim);
    return dvu;
};


inline arma::vec twolinkbase::zxu(const arma::vec& x,const arma::vec& u){
    arma::vec psixu = arma::join_cols(zx(x),zu(x,u));
    return psixu;
}; 

#endif