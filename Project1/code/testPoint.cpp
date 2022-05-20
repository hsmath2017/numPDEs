#include "Point.h"
int main(){
    Point p1,p2;
    double pos1[3]={0.5,3.5,4};
    double pos2[3]={1,2.4,3};
    double vel1[3]={0.6,3.5,10};
    double vel2[3]={4.5,5.5,10};
    double k=0.5;
    p1.setPosition(pos1);
    p1.setvelocity(vel1);
    p2.setPosition(pos2);
    p2.setvelocity(vel2);
    for(int i=0;i<3;i++){
	cout<<p1.getPosition()[i]<<",";
    }
    cout<<endl;
    cout<<p1<<endl;
    cout<<p2<<endl;
    Point add=p1+p2;
    Point minus=p1-p2;
    Point product=k*p1;
    cout<<add<<endl;
    cout<<product<<endl;
    cout<<minus<<endl;
}
