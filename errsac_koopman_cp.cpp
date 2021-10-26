#include <iostream>

#include <fstream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"src/cartpend.hpp"
#include"src/error_cost.hpp"
#include"src/SAC.hpp"
#include"src/rk4_int.hpp"

#include"src/koopsys.hpp"
#include"src/cpbasis.hpp"

arma::vec xdk(double t){
        return arma::zeros(6);};
arma::vec unom(double t){
        return arma::zeros(1);};

int main()
{   ofstream myfile;
    myfile.open ("test.csv");
 
 	CPBASIS basisobj;
 
 	KoopSys<CPBASIS> systK (0.01,&basisobj);
 
    CartPend syst1 (0.1,0.1,9.81,2.0,0.01);
    arma::mat R = 0.3*arma::eye(1,1);
 	arma::mat Qk = {
        {200,0.,0.,0.,0.,0.},
        {0., 0.,0.,0.,0.,0.},
        {0.,0.,20.,0.,0.,0.},
        {0.,0.,0.,1.,0.,0.},
		{0.,0.,0.,0.,0.,0.},
		{0.,0.,0.,0.,0.,0.}};
    arma::vec umax = {20};
 
    arma::vec xwrap,zwrap,uK;
    syst1.Ucurr = {0.0}; 
    syst1.Xcurr = {-3., 0.0,0.0,0.0};
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
 coeff.close();
}

