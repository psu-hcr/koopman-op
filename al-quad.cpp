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
    syst1.Ucurr = {0.,0.,0.,0.}; systK.Ucurr = syst1.Ucurr;
 	arma::vec anginit = (2*arma::randu<arma::vec>(3))-1;cout<<anginit<<endl;
	arma::mat Rinit = euler2R(anginit);
 	arma::vec pinit = {0.,0.,0.,1.};
 	arma::vec Twistinit = (2*arma::randu<arma::vec>(6))-1;
 	arma::mat hinit; hinit.zeros(4,4);
 	hinit.submat(0,0,2,2) = Rinit; hinit.submat(0,3,3,3)=pinit;
	
    syst1.Xcurr=arma::join_cols(hinit.as_col(),Twistinit);
 	systK.Xcurr = basisobj.zx(syst1.get_measurement(syst1.Xcurr));
 
	//set values for Q,R,Qf,umax
 	arma::mat R = arma::eye(syst1.Ucurr.n_rows,syst1.Ucurr.n_rows);
 	arma::mat Qk = arma::zeros(basisobj.xdim,basisobj.xdim);
	arma::vec Qvec = {1,1,1,1,1,1,5,5,5};
 	Qk.submat(0,0,8,8)=arma::diagmat(Qvec);
	arma::mat Qf = arma::zeros<arma::mat>(size(Qk));
    arma::vec umax(size(syst1.Ucurr)); umax.fill(6);
 	//errorcost<KoopSys<QuadBasis>> costK (Qk,R,xdk,&systK);
	arma::vec noisecov = 1.0*arma::ones(basisobj.xdim);
	arma::mat Rtil = 0.1*arma::eye(systK.Ucurr.n_rows,systK.Ucurr.n_rows);
	
    //sac<KoopSys<QuadBasis>,errorcost<KoopSys<QuadBasis>>> sacsysK (&systK,&costK,0.,1.0,umax,unom);
	lqr lqrK(Qk, R,Qf,20,umax,xdk, DT);
	fishcost<KoopSys<QuadBasis>,lqr> costFI (&systK,&lqrK,noisecov);
 	alk<KoopSys<QuadBasis>,fishcost<KoopSys<QuadBasis>,lqr>,lqr> ALpol(&systK,&costFI,&lqrK,T,umax,Rtil);
 myfile<<"time,q1,q2,q3,ag1,ag2,ag3,u1,u2,mu1,lqr\n";
 arma::vec measure,agK,mu;
 mu = arma::zeros(arma::size(umax));
 
	while (syst1.tcurr<5.0){
    myfile<<syst1.tcurr<<",";
    measure = syst1.get_measurement(syst1.Xcurr);
	agK=systK.Xcurr.subvec(0,2);
    myfile<<measure(0)<<","<<measure(1)<<","<<measure(2)<<",";
	myfile<<agK(0)<<","<<agK(1)<<","<<agK(2)<<",";
    myfile<<syst1.Ucurr(0)<<","<<syst1.Ucurr(1)<<","<<mu(0)<<",";
	myfile<<0.01*lqrK.l(systK.Xcurr,systK.Ucurr,systK.tcurr)<<"\n";
	syst1.step();
	systK.update_XU(measure,syst1.Ucurr);
	systK.calc_K();
	lqrK.calc_gains(systK.Kx,systK.Ku);mu = lqrK.mu(systK.Xcurr,syst1.tcurr);
	systK.step();
	syst1.Ucurr = ALpol.ustar_calc(); //unom(systK.tcurr);
    //sacsysK.unom_shift();
    if(fmod(syst1.tcurr,5)<syst1.dt)cout<<"Time: "<<syst1.tcurr<<"\n";
    } 
       
    myfile.close();
 
 ofstream coeff;
 coeff.open("CP-koopman.csv");
 systK.K.save(coeff,arma::csv_ascii);
 coeff.close();
}

