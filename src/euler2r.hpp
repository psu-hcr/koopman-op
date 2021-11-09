#ifndef EULER2R_HPPul
#define EULER2R_HPP
#include<armadillo>

arma::mat euler2R(const arma::vec& ang){
arma::mat Rz = {{cos(ang(2)),-sin(ang(2)),0.},
				{sin(ang(2)),cos(ang(2)),0.},
				{0.,0.,1.}};
arma::mat Rx= {{1.,0.,0.},
				{0.,cos(ang(0)),-sin(ang(0))},
				{0.,sin(ang(0)),cos(ang(0))}};
arma::mat Ry = {{cos(ang(1)),0.,sin(ang(1))},
				{0.,1.,0.},
			   {-sin(ang(1)),0.,cos(ang(1))}};
return Rx*Ry*Rz;
}

#endif