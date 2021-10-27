#include <iostream>

#include <fstream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"src/quadrotor.hpp"
#include"src/error_cost.hpp"
#include"src/SAC.hpp"
#include"src/rk4_int.hpp"

#include"src/koopsys.hpp"
#include"src/quadbasis.hpp"

arma::vec xdk(double t){//should match xdim defined in basis
		arma::vec ref = arma::zeros(18);
		ref(2) = 9.81;
        return ref;};
arma::vec unom(double t){
        return arma::zeros(4);};

int main()
{   arma::arma_rng::set_seed(50);//set seed for reproducibility
	//arma::arma_rng::set_seed_random();
 
	ofstream myfile;
    myfile.open ("test.csv");
 
 	QuadBasis basisobj;
 
 	KoopSys<QuadBasis> systK (0.01,&basisobj);
 
    QuadRotor syst1 (1./200.);
    syst1.Ucurr = {0.,0.,0.,0.}; 
 	arma::mat Rinit = arma::randu<arma::mat>(3,3);
 	arma::vec pinit = {1.,1.,1.,1.};
 	arma::vec Twistinit = arma::randu<arma::vec>(6);
 	arma::mat hinit; hinit.zeros(4,4);
 	hinit.submat(0,0,2,2) = Rinit; hinit.submat(0,3,3,3);
    syst1.Xcurr=arma::join_cols(hinit.as_col(),Twistinit);
 	systK.Xcurr = basisobj.zx(syst1.Xcurr);
 
 	arma::mat R = arma::eye(syst1.Ucurr.n_rows,syst1.Ucurr.n_rows);
 	arma::mat Qk = arma::zeros(basisobj.xdim,basisobj.xdim);
	arma::vec Qvec = {1,1,1,1,1,1,5,5,5};
 	Qk.submat(0,0,8,8)=arma::diagmat(Qvec);
    arma::vec umax = {20};
 /*
    arma::vec xwrap,zwrap,uK;
    
 	systK.Ucurr = {0.0}; 
    systK.Xcurr = basisobj.zx(syst1.Xcurr);
 	uK = systK.Ucurr;
 	zwrap=syst1.proj_func(systK.Xcurr);
    errorcost<KoopSys<CPBASIS>> costK (Qk,R,xdk,&systK);
    sac<KoopSys<CPBASIS>,errorcost<KoopSys<CPBASIS>>> sacsysK (&systK,&costK,0.,1.0,umax,unom);
    
       
    myfile<<"time,theta,thetadot,x,xdot,u,xK,thetaK\n";
 
    while (syst1.tcurr<60.0){
    myfile<<syst1.tcurr<<",";
    xwrap = syst1.proj_func(syst1.Xcurr); 
    myfile<<xwrap(0)<<","<<xwrap(1)<<",";
    myfile<<xwrap(2)<<","<<xwrap(3)<<",";
    myfile<<syst1.Ucurr(0)<<","<<zwrap(2)<<","<<zwrap(0)<<"\n";
	syst1.step();
	systK.update_XU(syst1.Xcurr,syst1.Ucurr);
	systK.calc_K();
	systK.step();zwrap=systK.proj_func(systK.Xcurr);
	sacsysK.SAC_calc();
	uK=sacsysK.ulist.col(0); 
	syst1.Ucurr = uK; 
    sacsysK.unom_shift();
    if(fmod(syst1.tcurr,5)<syst1.dt)cout<<"Time: "<<syst1.tcurr<<"\n";
    } 
       
    myfile.close();
 
 ofstream coeff;
 coeff.open("CP-koopman.csv");
 systK.K.save(coeff,arma::csv_ascii);
 coeff.close();*/
}

