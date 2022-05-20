#ifndef _HSMATH_CLASSICAL_RK
#define _HSMATH_CLASSICAL_RK
#include "TimeIntegrator.h"
class RungeKutta:public TimeIntegrator{
private:
    vector<double> RK_coefficient;
public:
    RungeKutta(){
        RK_coefficient.push_back(1.0/6);
	RK_coefficient.push_back(2.0/6);
	RK_coefficient.push_back(2.0/6);
	RK_coefficient.push_back(1.0/6);
    }
    void onestepiterator(){
	vector<Point> RKPoints;
	Point P=Pts[Pts.size()-1];
	Point start=P;
	Point P1=func(P,mu);
	RKPoints.push_back(P1);//it means y1
	for(int i=0;i<3;i++){
	    Point tmp=start+(i/2+1)*k/2*P1;
	    P1=func(tmp,mu);//y2,y3,y4
	    RKPoints.push_back(P1);
	}
        Point final=start;
	for(int i=0;i<4;i++){
	    final=final+k*RK_coefficient[i]*RKPoints[i];
	}
	Pts.push_back(final);
	RKPoints.clear();
    }
    void fiveorderRK(){
	vector<Point> RKPoints;
	Point P=Pts[Pts.size()-1];
	Point start=P;
	Point P1=func(P,mu);
	RKPoints.push_back(P1);//it means y1
	double b[6]={13.0/160,0,2375.0/5984,5.0/16,12.0/85,3.0/44};
	double c[6]={0,1.0/6,4.0/15,2.0/3,5.0/6,1};
	for(int i=1;i<6;i++){
	    Point tmp=start+c[i]*k*P1;
	    P1=func(tmp,mu);//y2,y3,y4
	    RKPoints.push_back(P1);
	}
        Point final=start;
	for(int i=0;i<6;i++){
	    final=final+k*b[i]*RKPoints[i];
	}
	Pts.push_back(final);
	RKPoints.clear();
    }
    void solveODE(){
	for(int i=0;i<T/k;i++){
	    onestepiterator();
	}
    }
    ~RungeKutta(){
    }
};

#else
//nothing
#endif
