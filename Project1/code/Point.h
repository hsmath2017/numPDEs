/**
 * @file   Point.h
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Sun Mar 28 03:43:54 2021
 * 
 * @brief  Store the information of the Position and Velocity vector of a point in the field.
 * 
 * 
 */
#ifndef _HSMATH_POINT
#define _HSMATH_POINT
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
class Point{
private:
    double Position[3];
    double velocity[3];
public:
    Point(){
    }
    ~Point(){
    }
    void setPosition(double (&Pos)[3]);
    double* getPosition();
    void setvelocity(double (&vel)[3]);
    double* getvelocity();
    const Point& operator=(const Point &_p);//get the value
    friend ostream& operator<<(ostream &_os,const Point &_p);
    friend Point operator+(const Point &_p1,const Point &_p2);//the addition of two vectors.
    friend Point operator*(double k,const Point &_p);//multiple a vector with a real number.
    friend Point operator-(const Point &_p1,const Point &_p2);
    double norm(){
	double n=0;
	for(int i=0;i<3;i++){
	    n=n+pow(Position[i],2);
	    n=n+pow(velocity[i],2);
	}
	n=pow(n,0.5);
	return n;
    }
    double maxnorm(){
	double n=0;
	for(int i=0;i<3;i++){
	    if(abs(Position[i])>n){
		n=abs(Position[i]);
	    }
	    if(abs(velocity[i])>n){
		n=abs(velocity[i]);
	    }
	}
	return n;
    }
};
Point operator+(const Point &_p1,const Point &_p2);//the addition of two vectors.
Point operator*(double k,const Point &_p);//multiple a vector with a real number.
Point operator-(const Point &_p1,const Point &_p2);
#else
//DO NOTHING
#endif
