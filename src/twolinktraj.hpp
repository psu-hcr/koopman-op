#ifndef twolinktraj_HPP
#define twolinktraj_HPP
#include<armadillo>

class twolinktraj{

	public:
	int xdim;
	arma::vec dtraj;
	int total_step;				// Total number of step for one period
	double period = 2*3.1415;	// Trajectory period
	double dt;					// Time step
	arma::mat djsmat;			// Matrix storing desire joint state
	double t_curr;				// current time
	
	twolinktraj(double _dt, int _xdim){
		// class initialization
		xdim = _xdim; 
		dt = _dt;
		total_step = period/dt;
		djsmat = arma::zeros(4,total_step);
		
	}
	
	arma::vec desire_traj(int step){
		// desire trajectory
		return dtraj;
	}
	
	
	arma::vec desire_jointstate(double t){
		// desire joint state
		arma::vec djs = arma::zeros(xdim);
		//arma::vec djs = arma::ones(xdim);
		djs(0) = sin(t); djs(1) = cos(t); //std::cout<<djs<<std::endl;
		djs(2) = 0.5*cos(t); djs(3) = -0.5*sin(t);
		return djs;
	}
		
};
#endif