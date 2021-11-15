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
#include"src/quadbasis.hpp"
#include"src/fisher_cost.hpp"
#include"src/ActLearnK.hpp"
#include"src/euler2r.hpp"

arma::vec xdk(double t){//should match xdim defined in basis
		arma::vec ref = arma::zeros(18);
		ref(2) = -9.81;
        return ref;};
arma::vec unom(double t){
        return arma::randn(4);};

int main()
{   arma::arma_rng::set_seed(50);//set seed for reproducibility
	//arma::arma_rng::set_seed_random();
 	
	ofstream myfile;
    myfile.open ("test.csv");
 
	double DT = 1./200.;
	double T = 0.1;
	//initialize Koopman system object and simulated quadrotor object
 	QuadBasis basisobj;
 	KoopSys<QuadBasis> systK (DT,&basisobj);
 	QuadRotor syst1 (DT);
	//initialize states and control for both systems
    syst1.Ucurr = {0.77,-0.77,-0.04,1.08}; systK.Ucurr = syst1.Ucurr;
 	arma::vec anginit = {0.95,0.78,0.82};//(2*arma::randu<arma::vec>(3))-1;
	arma::mat Rinit = euler2R(anginit);
 	arma::vec pinit = {0.,0.,0.,1.};
 	arma::vec Twistinit = {0.91,0.85,0.73,-0.52,0.52,0.74};//(2*arma::randu<arma::vec>(6))-1;
 	arma::mat hinit; hinit.zeros(4,4);
 	hinit.submat(0,0,2,2) = Rinit; hinit.submat(0,3,3,3)=pinit;
	
    syst1.Xcurr=arma::join_cols(hinit.as_col(),Twistinit);
 	systK.Xcurr = basisobj.zx(syst1.get_measurement(syst1.Xcurr));
 	
	//set values for Q,R,Qf,umax,noisecov,Regularization
 	arma::mat R = arma::eye(syst1.Ucurr.n_rows,syst1.Ucurr.n_rows);
 	arma::mat Qk = arma::zeros(basisobj.xdim,basisobj.xdim);
	arma::vec Qvec = {1,1,1,1,1,1,5,5,5};
 	Qk.submat(0,0,8,8)=10*arma::diagmat(Qvec);
	arma::mat Qf = arma::zeros<arma::mat>(size(Qk));
    arma::vec umax(size(syst1.Ucurr)); umax.fill(6);
 	arma::vec noisecov = 1.0*arma::ones(basisobj.xdim);
	arma::mat Rtil = 1.*arma::eye(systK.Ucurr.n_rows,systK.Ucurr.n_rows);
	//initialize lqr policy, fisher informaiton cost, and active learning controller
	lqr lqrK(Qk, R,Qf,round(T/DT),umax,xdk, DT);
	fishcost<KoopSys<QuadBasis>,lqr> costFI (&systK,&lqrK,noisecov);
 	alk<KoopSys<QuadBasis>,fishcost<KoopSys<QuadBasis>,lqr>,lqr> 		
							ALpol(&systK,&costFI,&lqrK,T,umax,Rtil);
 	
//set up file to store data
 myfile<<"time,q1,q2,q3,ag1,ag2,ag3,u1,u2,mu1,lqr\n";
 arma::vec measure,agK,mu;
 
 	//add initial conditions to state sample
 	measure = syst1.get_measurement(syst1.Xcurr);
 	systK.calc_K(measure,syst1.Ucurr);
 /*arma::mat kx = 5.*arma::eye(size(systK.Kx));
 arma::mat ku = 5.*arma::eye(size(systK.Ku));
 systK.Kx=kx; systK.Ku=ku;
 lqrK.calc_gains(kx,ku,systK.tcurr);
 mu = ALpol.ustar_calc();
 //cout<<mu<<endl;*/
 
	lqrK.calc_gains(systK.Kx,systK.Ku,systK.tcurr);
 	mu = lqrK.mu(systK.Xcurr,syst1.tcurr);
 
  
while (syst1.tcurr<60.){
    myfile<<syst1.tcurr<<",";
    agK=systK.Xcurr.subvec(0,2);
    myfile<<measure(0)<<","<<measure(1)<<","<<measure(2)<<",";
	myfile<<agK(0)<<","<<agK(1)<<","<<agK(2)<<",";
    myfile<<syst1.Ucurr(0)<<","<<syst1.Ucurr(1)<<","<<mu(0)<<",";
	myfile<<arma::norm(basisobj.zx(measure)-xdk(systK.tcurr))-10.<<"\n";
	syst1.step();
	measure = syst1.get_measurement(syst1.Xcurr);//sample state
	systK.calc_K(measure,syst1.Ucurr);//add to data set and update Kx, Ku
	lqrK.calc_gains(systK.Kx,systK.Ku,systK.tcurr);//update lqr gain
	mu = lqrK.mu(systK.Xcurr,systK.tcurr);//this is just to record mu
	syst1.Ucurr = ALpol.ustar_calc(); //compute ustar
	//systK.step();//this is just to record the model accuracy
    if(fmod(syst1.tcurr,2)<syst1.dt)cout<<"Time: "<<syst1.tcurr<<endl<<
		(systK.Xcurr).t()<<"\n";
    }
 	costFI.infw = 100.0 * pow(0.99,systK.tcurr);
       
    myfile.close();
 
 ofstream coeff;
 coeff.open("CP-koopman.csv");
 systK.K.save(coeff,arma::csv_ascii);
 coeff.close();
}

