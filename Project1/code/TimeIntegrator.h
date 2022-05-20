/**
 * @file   TimeIntegrator.h
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Sun Mar 28 00:41:45 2021
 * 
 * @brief  The time integrator to solve the ODEs numerically
 * 
 * 
 */ 

#ifndef HSMATH_TIMEINTEGRATOR
#define HSMATH_TIMEINTEGRATOR
#include "Function.h"
#include <string>
#include <iomanip>
#include <vector>
#include <stdexcept>
class TimeIntegrator{
protected:
    vector<Point> Pts;//Point in the integrator
    double mu;//parameter
    double T;//we need to solve the position on T
    double k;//time step
    int p;//Precision
    int s;//Steps(for LMM)
    vector<double> alpha;//the \alpha parameter of the LMMs
    vector<double> beta;//the \beta parameter of the LMMs
    int type;//Numerical Methods
public:
    TimeIntegrator();
    void SetInfo(double* Input,double step,double _T,double _mu,int p);
    void ForwardEuleronestep();
    virtual void onestepiterator()=0;//one_step_iterator
    virtual void solveODE()=0;//ODE_Solver!
    void PrintSolution();//Output the numerical solution.
    void DrawFigure();//Plot the track
    double calculateRelativeError(){//For the two test cases, we know that the solution satisfies the periodical property
	Point delta=Pts[Pts.size()-1]-Pts[0];
	return delta.norm()/Pts[0].norm();
    }
    double calculateAbsoluteError(){//Absolute Error
	Point delta=Pts[Pts.size()-1]-Pts[0];
	return delta.norm();
    }
    double calculateMaxNormError(){//Max-Norm
	Point delta=Pts[Pts.size()-1]-Pts[0];
	return delta.maxnorm();
    }
    Point getSolution(){
	return Pts[Pts.size()-1];
    }
};
void TimeIntegrator::SetInfo(double* Input,double step,double _T,double _mu,int p){
    Point _P;
    double pos[3];
    double vel[3];
    for(int i=0;i<3;i++){
	pos[i]=Input[i];
	vel[i]=Input[i+3];
    }
    _P.setPosition(pos);
    _P.setvelocity(vel);
    Pts.push_back(_P);//The element in P
    k=step;
    T=_T;
    mu=_mu;
    this->p=p;
}
TimeIntegrator::TimeIntegrator(){
    //Nothing
}
void TimeIntegrator::ForwardEuleronestep(){
    Point Pnew;
    Point tmp;
    tmp=func(Pts[Pts.size()-1],mu);
    tmp=k*tmp;
    Pnew=tmp+Pts[Pts.size()-1];
    Pts.push_back(Pnew);
    return;
}
void TimeIntegrator::PrintSolution(){
    cout<<Pts[Pts.size()-1]<<endl;
}
/** 
 * Output a .m file to draw the figure of track
 * 
 */
void TimeIntegrator::DrawFigure(){
    ofstream os;
    os.open("Trajectory.m");
    if(1){//os is empty!
	os<<"x = [";
	for(int i=0;i<Pts.size();i++){
	    os<<Pts[i].getPosition()[0]<<",";
	}
	os<<"];\n";
	os<<"y = [";
	for(int i=0;i<Pts.size();i++){
	    os<<Pts[i].getPosition()[1]<<",";
	}
	os<<"];\n";
    }
    os<<"plot(x,y);";
    os.close();
}

#else
//Do Nothing!
#endif
