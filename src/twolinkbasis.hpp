#ifndef twolinkbasis_HPP
#define twolinkbasis_HPP
#include<armadillo>
#include<cmath>


class twolinkbasis {
	
	// mass, length and Radius
	double m1 = 1, m2 = 1, L1 = 0.5, L2 = 0.3, R1 = 0.1, R2 = 0.1;
	
	// simplified coefficient
	double a = 0.166667*m2*pow(L2,2) + 0.5*m2*pow(R2,2);
	double b = 0.25*m1*pow(R1,2);
	double c = 0.41667*m2*pow(L2,2);
	
	public:
		arma::vec proj_func (const arma::vec& x);
        arma::vec zxu(const arma::vec& x, const arma::vec& u);
		arma::vec zx(const arma::vec& x);
		arma::vec zu(const arma::vec& x, const arma::vec& u);
		arma::mat dvdu(const arma::vec& z);
		double angle_wrap(double input);
		int zdim = 10;// set to length of z including v(x,u)
		int xdim = 8;// set to length of z(x(t)) only 
		int udim = zdim-xdim;
        
};

double twolinkbasis::angle_wrap (double input){
	// angle wrapping function
    double output = fmod(input+M_PI, 2*M_PI);
    if (output < 0.0) output = output + 2*M_PI;
    output = output - M_PI;
    return output;
}

arma::vec twolinkbasis::proj_func (const arma::vec& x){
	arma::vec xwrap=x;
    return xwrap;
}

arma::vec twolinkbasis::zx (const arma::vec& x){// note the order of x with a comment here
	arma::vec pos = {
		x(0),
		x(2)
	};
	arma::vec vel = {
		x(1),
		x(3)
	};
    arma::vec psix = {
		pos(0),
		vel(0),
		pos(1),
		vel(1),
		/*pos(0)*pos(1),
		pow(pos(0)*pos(1), 2),
		pow(pos(0)*pos(1), 3),
		vel(0)*vel(1),
		pow(vel(0)*vel(1), 2),
		pow(vel(0)*vel(1), 3),
		1*/
		sin(2*x[2])/(pow(sin(x[2]),2)+b/a),
		x[1]*(2*a*cos(2*x[2])*(a*pow(sin(x[2]),2)+1) - a*sin(2*x[2])*(2*a*sin(x[2])*cos(x[2])))/
			pow((a*pow(sin(x[2]),2) + b),2),
		sin(2*x[2]),
		x[1]*cos(2*x[2])
	};
	return psix;
}
arma::vec twolinkbasis::zu(const arma::vec& x, const arma::vec& u){//recommend psiu=u;
    arma::vec psiu = {u};
    return psiu;
}; 

arma::mat twolinkbasis::dvdu(const arma::vec& zx){//if psiu=u this is the identity
     arma::mat dvu= arma::eye(udim, udim);
    return dvu;
};


inline arma::vec twolinkbasis::zxu(const arma::vec& x,const arma::vec& u){
    arma::vec psixu = arma::join_cols(zx(x),zu(x,u));
    return psixu;
}; 

#endif