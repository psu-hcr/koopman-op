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

arma::vec xd(double t){
        return arma::zeros(4);};
arma::vec unom(double t){
        return arma::zeros(1);};

int main()
{   ofstream myfile;
    myfile.open ("test.csv");
    CartPend syst1 (0.1,0.1,9.81,2.0,0.01);
    arma::mat Q = {
        {200,0.,0.,0.},
        {0., 0.,0.,0.},
        {0.,0.,20.,0.},
        {0.,0.,0.,1.}};
    arma::mat R = 0.3*arma::eye(1,1);
    arma::vec umax = {20};
 
    arma::vec xwrap;
    syst1.Ucurr = {0.0}; 
    syst1.Xcurr = {-3.1, 0.0,0.0,0.0};
    errorcost<CartPend> cost (Q,R,xd,&syst1);
    sac<CartPend,errorcost<CartPend>> sacsys (&syst1,&cost,0.,1.0,umax,unom);
    //arma::mat unom = arma::zeros<arma::mat>(1,sacsys.T_index);
       
    myfile<<"time,theta,thetadot,x,xdot,u\n";
 
    while (syst1.tcurr<30.0){
    myfile<<syst1.tcurr<<",";
    xwrap = syst1.proj_func(syst1.Xcurr); 
    myfile<<xwrap(0)<<","<<xwrap(1)<<",";//myfile<<syst1.Xcurr(0)<<","<<syst1.Xcurr(1)<<",";
    myfile<<xwrap(2)<<","<<xwrap(3)<<",";//myfile<<syst1.Xcurr(2)<<","<<syst1.Xcurr(3)<<",";
    myfile<<syst1.Ucurr(0)<<"\n";
    syst1.step();
    sacsys.SAC_calc();
    syst1.Ucurr = sacsys.ulist.col(0); 
    sacsys.unom_shift();
    if(fmod(syst1.tcurr,5)<syst1.dt)cout<<"Time: "<<syst1.tcurr<<"\n";
    } 
       
    myfile.close();
}

