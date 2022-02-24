#include <iostream>
#include <fstream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"src/quadrotor.hpp"
#include"src/error_cost.hpp"
//#include"src/SAC.hpp"
#include"src/rk4_int.hpp"
#include"src/LQR.hpp"
#include"src/koopsys.hpp"
#include"src/fisher_cost.hpp"
#include"src/ActLearnK.hpp"
#include"src/euler2r.hpp"
#include"src/twomass.hpp"
#include"src/twomassbasis.hpp"

arma::vec xdk(double t){//should match xdim defined in basis
	//arma::vec ref = arma::zeros(2);
	arma::vec ref = arma::zeros(4);
	//ref(0) = sin(t); ref(1) = cos(t);
	//ref(2) = sin(t); ref(3) = cos(t);
	return ref;
};

arma::vec unom(double t){
        //return {0.1*sin(t), 0.1*sin(t)};
		return {0.1, 0.1};
};

int main(){   
	arma::arma_rng::set_seed(40);//set seed for reproducibility
	//arma::arma_rng::set_seed_random();
 	
	ofstream myfile;
    myfile.open ("/home/lu/koopman-op/twomass.csv");
 
	double DT = 1e-3;
	double T = 0.1;
	//initialize Koopman system object and simulated quadrotor object
 	twomassBasis basisobj;
 	KoopSys<twomassBasis> systK (DT,&basisobj);
 	twomass syst1 (DT);
	//initialize states and control for both systems
    syst1.Ucurr = {0.77, 0.77};//arma::randn(4); 
    systK.Ucurr = syst1.Ucurr;
    //syst1.Xcurr = {0.1, 0.1};
	syst1.Xcurr = {1.0, 0, 1.0, 0};
 	systK.Xcurr = basisobj.zx(syst1.get_measurement(syst1.Xcurr));
	
	//set values for Q,R,Qf,umax,noisecov,Regularization
 	arma::mat R = arma::eye(syst1.Ucurr.n_rows,syst1.Ucurr.n_rows);
 	arma::mat Qk = arma::zeros(basisobj.xdim,basisobj.xdim);
	//arma::vec Qvec = {1,1};
 	//Qk.submat(0,0,1,1)=10*arma::diagmat(Qvec);
	arma::vec Qvec = {1,1,1,1}; 
 	Qk.submat(0,0,3,3)=100*arma::diagmat(Qvec);
	arma::mat Qf = arma::zeros<arma::mat>(size(Qk));
    arma::vec umax(size(syst1.Ucurr)); umax.fill(6);
 	arma::vec noisecov = 1.0*arma::ones(basisobj.xdim);
	arma::mat Rtil = 10.*arma::eye(systK.Ucurr.n_rows,systK.Ucurr.n_rows);
	//initialize lqr policy, fisher informaiton cost, and active learning controller
	lqr lqrK(Qk, R,Qf,round(T/DT),umax,xdk, DT);
	fishcost<KoopSys<twomassBasis>,lqr> costFI (&systK,&lqrK,noisecov);
 	alk<KoopSys<twomassBasis>,fishcost<KoopSys<twomassBasis>,lqr>,lqr> 		
							ALpol(&systK,&costFI,&lqrK,T,umax,Rtil);
 	
	//set up file to store data
 	myfile<<"time,x1,x2,x3,x4,u1,u2,mu1,mu2,lqr\n";
 	arma::vec measure,agK,mu;
 
 	//add initial conditions to state sample
 	measure = syst1.get_measurement(syst1.Xcurr); 
 	systK.update_XU(measure,syst1.Ucurr); 
 	
 	systK.Kx = arma::randn<arma::mat>(basisobj.xdim,basisobj.xdim); 
 	systK.Ku = arma::randn<arma::mat>(basisobj.xdim,basisobj.zdim-basisobj.xdim); 
	lqrK.calc_gains(systK.Kx,systK.Ku,systK.tcurr); 
 	mu = lqrK.mu(systK.Xcurr,syst1.tcurr); 
 
	while (syst1.tcurr<60.){
		myfile<<syst1.tcurr<<",";
		myfile<<measure(0)<<","<<measure(1)<<","<<measure(2)<<","<<measure(3)<<",";
		myfile<<syst1.Ucurr(0)<<","<<syst1.Ucurr(1)<<",";
		myfile<<mu(0)<<","<<mu(1)<<",";
		myfile<<arma::norm(basisobj.zx(measure)-xdk(systK.tcurr))-20.<<"\n";
		syst1.step();
		measure = syst1.get_measurement(syst1.Xcurr);//sample state
		systK.calc_K(measure,syst1.Ucurr);//add to data set and update Kx, Ku
		lqrK.calc_gains(systK.Kx,systK.Ku,systK.tcurr);//update lqr gain
		mu = lqrK.mu(systK.Xcurr,systK.tcurr);//this is just to record mu
		syst1.Ucurr = ALpol.ustar_calc(); //compute ustar
		//syst1.Ucurr =  lqrK.mu(systK.Xcurr,systK.tcurr);
		if(syst1.Ucurr(0)!=syst1.Ucurr(0)){cout<<"returned a nan"<<endl;
			syst1.Ucurr = unom(syst1.tcurr);
		}
		//systK.step();//this is just to record the model accuracy
		if(fmod(syst1.tcurr,2)<syst1.dt){
			cout<<"Time: "<<syst1.tcurr<<endl<<
			(systK.Xcurr).t()<<"\n"<<lqrK.dmudz(systK.Xcurr,systK.tcurr)<<"\n";
		}
		
		// switch informaiton cost gain
		if(systK.tcurr<20){
			costFI.infw = 100.;
		}
		else{
			costFI.infw = 0.;
		}
		//costFI.infw = 100 * pow(0.8,systK.tcurr); //0.0f;
	}
	myfile.close();

	 ofstream coeff;
	 coeff.open("/home/lu/koopman-op/twomass-koopman.csv");
	 systK.K.save(coeff,arma::csv_ascii);
	 coeff.close();
}