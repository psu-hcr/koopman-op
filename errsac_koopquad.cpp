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
		ref(2) = -9.81;
        return ref;};
arma::vec unom(double t){
        return arma::randn(4);};

int main()
{   //arma::arma_rng::set_seed(50);//set seed for reproducibility
	arma::arma_rng::set_seed_random();
 
	ofstream myfile;
    myfile.open ("test.csv");
 	double DT = 1./200.;
 	QuadBasis basisobj;
 
 	KoopSys<QuadBasis> systK (DT,&basisobj);
 
    QuadRotor syst1 (DT);
    syst1.Ucurr = {0.,0.,0.,0.}; systK.Ucurr = syst1.Ucurr;
 	arma::mat Rinit = arma::normalise(arma::randn<arma::mat>(3,3));
 	arma::vec pinit = {1.,1.,1.,1.};
 	arma::vec Twistinit = arma::randn<arma::vec>(6);
 	arma::mat hinit; hinit.zeros(4,4);
 	hinit.submat(0,0,2,2) = Rinit; hinit.submat(0,3,3,3)=pinit;
    syst1.Xcurr=arma::join_cols(hinit.as_col(),Twistinit);
 	systK.Xcurr = basisobj.zx(syst1.get_measurement(syst1.Xcurr));
 
 	arma::mat R = arma::eye(syst1.Ucurr.n_rows,syst1.Ucurr.n_rows);
 	arma::mat Qk = arma::zeros(basisobj.xdim,basisobj.xdim);
	arma::vec Qvec = {1,1,1,1,1,1,5,5,5};
 	Qk.submat(0,0,8,8)=arma::diagmat(Qvec);
    arma::vec umax(size(syst1.Ucurr)); umax.fill(6);
 	errorcost<KoopSys<QuadBasis>> costK (Qk,R,xdk,&systK);
    sac<KoopSys<QuadBasis>,errorcost<KoopSys<QuadBasis>>> sacsysK (&systK,&costK,0.,1.0,umax,unom);
 
 myfile<<"time,ag1,ag2,ag3,u1,u2,u3,ag3K\n";
 arma::vec measure,agK;
 
 
	while (syst1.tcurr<10.0){
    myfile<<syst1.tcurr<<",";
    measure = syst1.get_measurement(syst1.Xcurr);
	agK=systK.Xcurr.subvec(0,2);
    myfile<<measure(0)<<","<<measure(1)<<",";
    myfile<<measure(2)<<","<<syst1.Ucurr(0)<<",";
    myfile<<syst1.Ucurr(1)<<","<<syst1.Ucurr(2)<<","<<agK(2)<<"\n";
	syst1.step();
	systK.update_XU(measure,syst1.Ucurr);
	systK.calc_K();
	systK.step();
	sacsysK.SAC_calc();
	syst1.Ucurr = sacsysK.ulist.col(0); 
    sacsysK.unom_shift();
    if(fmod(syst1.tcurr,5)<syst1.dt)cout<<"Time: "<<syst1.tcurr<<"\n";
    } 
       
    myfile.close();
 
 ofstream coeff;
 coeff.open("CP-koopman.csv");
 systK.K.save(coeff,arma::csv_ascii);
 coeff.close();
}

