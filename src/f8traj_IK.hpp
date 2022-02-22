#ifndef F8TRAJIK_HPP
#define F8TRAJIK_HPP
#include<armadillo>

template <class system>
class f8traj_IK{

	system* sys;
	public:
	int xdim;
	arma::vec dtraj = arma::zeros(7);
	int total_step;				// Total number of step for one period
	double period = 2*3.1415;	// Trajectory period
	double dt;					// Time step
	arma::mat djsmat;			// Matrix storing desire joint state
	double t_curr;				// current time
	
	// class initialization
	f8traj_IK(system *_sys){
		sys = _sys;
		xdim = sys->zfuncs->xdim; 
		dt = sys->dt;
		total_step = period/dt;
		djsmat = arma::zeros(7,total_step);
		
	}
	
	// desire trajectory
	arma::vec desire_traj(int step){
		double t = step*dt;
		arma::vec dtraj = arma::zeros(7);
		dtraj[0] = 0.6;						//x
		dtraj[1] = 0.25*cos(1*t);			//y
		dtraj[2] = 0.1*sin(2*t)+0.7;		//z
		dtraj[3] = 1;						//wx
		dtraj[4] = 0;						//wy
		dtraj[5] = 0;						//wz
		dtraj[6] = 0;						//w
		return dtraj;
	}
	
	// desire joint state
	arma::vec desire_jointstate(double t){
		arma::vec djs = arma::zeros(xdim);
		int current_step = int(t/dt) % total_step;
		//djs.subvec(0, 6) = djsmat.col(current_step);
		djs.subvec(0, 6) = djsmat.col(current_step);	// debugging: move with only joint1
		djs(1) = 0; djs(2) = 0; djs(3) = 0; djs(4) = 0; djs(5) = 0; djs(6) = 0;
		//std::cout<<djs<<std::endl;
		return djs;
	}
		
};
#endif