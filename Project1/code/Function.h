/**
 * @file   Function.h
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Tue Apr 20 21:47:41 2021
 * 
 * @brief  Store the Right-Hand function f(u,t) and its Jacobian Matrix
 * 
 * 
 */
#ifndef _HSMATH_FUNCTION
#define _HSMATH_FUNCTION
#include"Point.h"
#include<Eigen/Dense>
using namespace Eigen;
Point func(Point &P,double mu){
    Point P1;
    double oldpos[3];
    double oldvel[3];
    double newpos[3];
    double newvel[3];
    for(int i=0;i<3;i++){
	oldpos[i]=P.getPosition()[i];
	oldvel[i]=P.getvelocity()[i];
    }
    for(int i=0;i<3;i++){
	newpos[i]=oldvel[i];
    }
    newvel[0]=2*oldvel[1]+oldpos[0]-mu*(oldpos[0]+mu-1)
	/pow(pow(oldpos[1],2)+pow(oldpos[2],2)+pow(oldpos[0]+mu-1,2),1.5)
	-(1-mu)*(oldpos[0]+mu)/pow(pow(oldpos[1],2)+pow(oldpos[2],2)+pow(oldpos[0]+mu,2),1.5);
    newvel[1]=-2*oldvel[0]+oldpos[1]-mu*oldpos[1]
	/pow(pow(oldpos[1],2)+pow(oldpos[2],2)+pow(oldpos[0]+mu-1,2),1.5)
	-(1-mu)*oldpos[1]/pow(pow(oldpos[1],2)+pow(oldpos[2],2)+pow(oldpos[0]+mu,2),1.5);
    newvel[2]=-mu*oldpos[2]/pow(pow(oldpos[1],2)+pow(oldpos[2],2)+pow(oldpos[0]+mu-1,2),1.5)
	-(1-mu)*oldpos[2]/pow(pow(oldpos[1],2)+pow(oldpos[2],2)+pow(oldpos[0]+mu,2),1.5);
    P1.setPosition(newpos);
    P1.setvelocity(newvel);
    return P1;
}
Matrix<double,6,6> Jacobian(Point &P,double mu){
    Matrix<double,6,6> J;
    for(int i=0;i<6;i++){
	for(int j=0;j<6;j++){
	    J(i,j)=0;//Initialize
	}
    }
    double oldpos[3];
    double oldvel[3];
    for(int i=0;i<3;i++){
	oldpos[i]=P.getPosition()[i];
	oldvel[i]=P.getvelocity()[i];
    }
    J(0,3)=1;
    J(1,4)=1;
    J(2,5)=1;
    J(3,4)=2;
    J(4,3)=-2;
    double u1=oldpos[0];
    double u2=oldpos[1];
    double u3=oldpos[2];
    double denominator1=pow(pow(u1+mu,2)+pow(u2,2)+pow(u3,2),2.5);
    double denominator2=pow(pow(u1+mu-1,2)+pow(u2,2)+pow(u3,2),2.5);
    J(3,0)=1+(1-mu)*(2*pow((u1+mu),2)-pow(u2,2)-pow(u3,2))/denominator1
	+mu*(2*pow(u1+mu-1,2)-pow(u2,2)-pow(u3,2))/denominator2;
    J(3,1)=3*mu*u2*(mu+u1-1)/denominator2-3*(mu+u1)*u2*(mu-1)/denominator1;
    J(3,2)=3*mu*u3*(mu+u1-1)/denominator2-3*(mu+u1)*u3*(mu-1)/denominator1;
    J(4,0)=3*mu*u2*(mu+u1-1)/denominator2-3*u2*(mu+u1)*(mu-1)/denominator1;
    J(4,1)=(1-mu)*(2*pow(u2,2)-pow(mu+u1,2)-pow(u3,2))/denominator1+mu*(2*pow(u2,2)-pow(mu+u1-1,2)-pow(u3,2))/denominator2+1;
    J(4,2)=3*mu*u2*u3/denominator2-3*u2*u3*(mu-1)/denominator1;
    J(5,0)=3*mu*u3*(mu+u1-1)/denominator2-3*u3*(mu+u1)*(mu-1)/denominator1;
    J(5,1)=3*mu*u2*u3/denominator2-3*u2*u3*(mu-1)/denominator1;
    J(5,2)=(1-mu)*(2*pow(u3,2)-pow(mu+u1,2)-pow(u2,2))/denominator1+mu*(2*pow(u3,2)-pow(mu+u1-1,2)-pow(u2,2))/denominator2;
    return J;
}
#else
//NOTHING
#endif
