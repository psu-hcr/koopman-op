#ifndef iiwaBasis3_HPP
#define iiwaBasis3_HPP
#include<armadillo>
#include<cmath>


class iiwaBasis3 {
	
	public:
		arma::vec proj_func (const arma::vec& x);
        arma::vec zxu(const arma::vec& x, const arma::vec& u);
		arma::vec zx(const arma::vec& x);
		arma::vec zu(const arma::vec& x, const arma::vec& u);
		arma::mat dvdu(const arma::vec& z);
		double angle_wrap(double input);
		int zdim = 58;//set to length of z including v(x,u)
		int xdim = 51;//set to length of z(x(t)) only        
        
};

double iiwaBasis3::angle_wrap (double input){
    double output = fmod(input+M_PI, 2*M_PI);
    if (output < 0.0) output = output + 2*M_PI;
    output = output - M_PI;
    return output;
}

arma::vec iiwaBasis3::proj_func (const arma::vec& x){//angle wrapping function
    arma::vec xwrap=x;
	xwrap(0) = angle_wrap(x(0));
	xwrap(1) = angle_wrap(x(1));
	xwrap(2) = angle_wrap(x(2));
	xwrap(3) = angle_wrap(x(3));
	xwrap(4) = angle_wrap(x(4));
	xwrap(5) = angle_wrap(x(5));
	xwrap(6) = angle_wrap(x(6));
	xwrap(14) = xwrap(0)*xwrap(1);
	xwrap(15) = xwrap(1)*xwrap(2);
	xwrap(16) = xwrap(2)*xwrap(3);
	xwrap(17) = xwrap(3)*xwrap(4);
	xwrap(18) = xwrap(4)*xwrap(5);
	xwrap(19) = xwrap(5)*xwrap(6);
	xwrap(20) = pow(xwrap(0), 2)*pow(xwrap(1), 2);
	xwrap(21) = pow(xwrap(1), 2)*pow(xwrap(2), 2);
	xwrap(22) = pow(xwrap(2), 2)*pow(xwrap(3), 2);
	xwrap(23) = pow(xwrap(3), 2)*pow(xwrap(4), 2);
	xwrap(24) = pow(xwrap(4), 2)*pow(xwrap(5), 2);
	xwrap(25) = pow(xwrap(5), 2)*pow(xwrap(6), 2);
	xwrap(26) = pow(xwrap(0), 3)*pow(xwrap(1), 3);
	xwrap(27) = pow(xwrap(1), 3)*pow(xwrap(2), 3);
	xwrap(28) = pow(xwrap(2), 3)*pow(xwrap(3), 3);
	xwrap(29) = pow(xwrap(3), 3)*pow(xwrap(4), 3);
	xwrap(30) = pow(xwrap(4), 3)*pow(xwrap(5), 3);
	xwrap(31) = pow(xwrap(5), 3)*pow(xwrap(6), 3);
    
    return xwrap;
}

arma::vec iiwaBasis3::zx (const arma::vec& x){//note the order of x with a comment here
	arma::vec pos = {x(0),
					 x(1),
					 x(2),
					 x(3),
					 x(4),
					 x(5),
					 x(6)};
	arma::vec vel = {x(7),
					 x(8),
					 x(9),
					 x(10),
					 x(11),
					 x(12),
					 x(13)};
    arma::vec psix = {pos(0),
					  pos(1),
					  pos(2),
					  pos(3),
					  pos(4),
					  pos(5),
					  pos(6),
					  vel(0),
					  vel(1),
					  vel(2),
					  vel(3),
					  vel(4),
					  vel(5),
					  vel(6),
					  pos(0)*pos(1),
					  pos(1)*pos(2),
					  pos(2)*pos(3),
					  pos(3)*pos(4),
					  pos(4)*pos(5),
					  pos(5)*pos(6),
					  pow(pos(0), 2)*pow(pos(1), 2),
					  pow(pos(1), 2)*pow(pos(2), 2),
					  pow(pos(2), 2)*pow(pos(3), 2),
					  pow(pos(3), 2)*pow(pos(4), 2),
					  pow(pos(4), 2)*pow(pos(5), 2),
					  pow(pos(5), 2)*pow(pos(6), 2),
					  pow(pos(0), 3)*pow(pos(1), 3),
					  pow(pos(1), 3)*pow(pos(2), 3),
					  pow(pos(2), 3)*pow(pos(3), 3),
					  pow(pos(3), 3)*pow(pos(4), 3),
					  pow(pos(4), 3)*pow(pos(5), 3),
					  pow(pos(5), 3)*pow(pos(6), 3),
					  vel(0)*vel(1),
					  vel(1)*vel(2),
					  vel(2)*vel(3),
					  vel(3)*vel(4),
					  vel(4)*vel(5),
					  vel(5)*vel(6),
					  pow(vel(0), 2)*pow(vel(1), 2),
					  pow(vel(1), 2)*pow(vel(2), 2),
					  pow(vel(2), 2)*pow(vel(3), 2),
					  pow(vel(3), 2)*pow(vel(4), 2),
					  pow(vel(4), 2)*pow(vel(5), 2),
					  pow(vel(5), 2)*pow(vel(6), 2),
					  pow(vel(0), 3)*pow(vel(1), 3),
					  pow(vel(1), 3)*pow(vel(2), 3),
					  pow(vel(2), 3)*pow(vel(3), 3),
					  pow(vel(3), 3)*pow(vel(4), 3),
					  pow(vel(4), 3)*pow(vel(5), 3),
					  pow(vel(5), 3)*pow(vel(6), 3),
					  1
					  };
	return psix;
}
arma::vec iiwaBasis3::zu(const arma::vec& x, const arma::vec& u){//recommend psiu=u;
    arma::vec psiu = {u};
    return psiu;
}; 

arma::mat iiwaBasis3::dvdu(const arma::vec& z){//if psiu=u this is the identity
     arma::mat dvu= 1.*arma::eye(7,7);
    return dvu;
};


inline arma::vec iiwaBasis3::zxu(const arma::vec& x,const arma::vec& u){
    arma::vec psixu = arma::join_cols(zx(x),zu(x,u));
    return psixu;
}; 

#endif