#include "Point.h"
void Point::setPosition(double (&Pos)[3]){
    for(int i=0;i<3;i++){
	Position[i]=Pos[i];
    }
}
double* Point::getPosition(){
    double* pt=Position;
    return pt;
}
void Point::setvelocity(double (&vel)[3]){
    for(int i=0;i<3;i++){
	velocity[i]=vel[i];
    }
}
double* Point::getvelocity(){
    double* pt=velocity;
    return pt;
}
const Point& Point::operator=(const Point &_p){
    if(this!=&_p){
	for(int i=0;i<3;i++){
	    Position[i]=_p.Position[i];
	    velocity[i]=_p.velocity[i];
	}
    }
    return *this;
}
ostream& operator<<(ostream &_os,const Point &_p){
    _os<<"position:["<<_p.Position[0]<<","<<_p.Position[1]<<","<<_p.Position[2]<<"],";
    _os<<"velocity:["<<_p.velocity[0]<<","<<_p.velocity[1]<<","<<_p.velocity[2]<<"]";
    return _os;			
}
Point operator+(const Point &_p1,const Point &_p2){
    Point P;
    double Pos[3];
    double vel[3];
    for(int i=0;i<3;i++){
	Pos[i]=_p1.Position[i]+_p2.Position[i];
	vel[i]=_p1.velocity[i]+_p2.velocity[i];
    }
    P.setPosition(Pos);
    P.setvelocity(vel);
    return P;
}
Point operator*(double k,const Point &_p){
    Point P;
    double Pos[3];
    double vel[3];
    for(int i=0;i<3;i++){
	Pos[i]=k*_p.Position[i];
	vel[i]=k*_p.velocity[i];
    }
    P.setPosition(Pos);
    P.setvelocity(vel);
    return P;
}
Point operator-(const Point &_p1,const Point &_p2){
    Point P;
    double Pos[3];
    double vel[3];
    for(int i=0;i<3;i++){
	Pos[i]=_p1.Position[i]-_p2.Position[i];
	vel[i]=_p1.velocity[i]-_p2.velocity[i];
    }
    P.setPosition(Pos);
    P.setvelocity(vel);
    return P;
}