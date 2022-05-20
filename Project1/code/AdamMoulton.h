#ifndef _HSMATH_ADAMMOULTON
#define _HSMATH_ADAMMOULTON
#include <Eigen/LU>
#include "RungeKutta.h"
#include <map>
class AdamMoulton:public RungeKutta{
private:
    map<int,vector<double>> coefficient_table;
public:   
    AdamMoulton(){
    }
    void read_table();//Read the coefficient table of vector \alpha and \beta
    void generate_initial();//Generate the initial_values
    void onestepiterator();//The one_step iterator for A-B method
    void solveODE();//The ode_solver for A-B method
    Point Newton_Iterator();//Newton's Iterator method to solve unlinear equation
    Point Hfunc(Point& P);//Iterator function
    ~AdamMoulton(){
    }
};
Point AdamMoulton::Hfunc(Point &P){
    Point H;
    H=P-Pts[Pts.size()-1]-k*beta[0]*func(P,mu);
    for(int i=1;i<=s;i++){
	H=H-k*beta[i]*func(Pts[Pts.size()-i],mu);
    }
    return H;
}
void AdamMoulton::read_table(){//read the coefficient table,similar with A-B method!
    string tmp;
    ifstream in("AdamMoultonTable");
    getline(in,tmp);//the head of table
    int flag=0;
    while(!in.eof()&&flag<p){
	int s,p;//s:step num,p:precision,t:beta_s,always zero
	in>>s>>p;
	vector<double> vec;
	for(int i=0;i<p;i++){
	    string fraction;//read the coefficient:fraction frame
	    in>>fraction;
	    int pos;
	    for(int i=0;i<fraction.length();i++){
		if(fraction[i]=='/'){
		    pos=i;
		    break;
		}
	    }
	    int numerator,denominator;
	    numerator=stoi(fraction.substr(0,pos));
	    denominator=stoi(fraction.substr(pos+1));
	    double num;
	    num=numerator*1.0/denominator;
	    vec.push_back(num);
	}
	coefficient_table[p]=vec;//precition=p means beta=vec
	flag++;
    }
}
void AdamMoulton::generate_initial(){
    int init_num=p;//Need the number of initial value!
    if(p==1){
	s=p;
    }else{
	s=p-1;
    }
    for(int i=0;i<s-1;i++){
        RungeKutta::onestepiterator();
    }//s_steps Forward-Euler method to generate the initial_value!
}
void AdamMoulton::onestepiterator(){//Go forward one-step
    Point Pnew;
    Pnew=Newton_Iterator();
    Pts.push_back(Pnew);//one-step forward NOW!
}
void AdamMoulton::solveODE(){//Solve the ODE
    read_table();
    generate_initial();
    for(int i=0;i<coefficient_table[p].size();i++){
	beta.push_back(coefficient_table[p][i]);//Input the beta vector
    }
    for(int i=p;i<T/k;i++){
	onestepiterator();
    }
}
Point AdamMoulton::Newton_Iterator(){
    Matrix<double,6,6> Jacob;
    Matrix<double,6,6> I;//Identity Matrix
    I.setIdentity(6,6);//identity matrix
    Point sol=Pts[Pts.size()-1];//Guess:The solution is near the point "sol"
    Point H=Hfunc(sol);
    int looptime=0;
    while(H.norm()>1.0e-8&&looptime<4){//Newton Iterator!
	Jacob=I-k*beta[0]*Jacobian(sol,mu);
	Matrix<double,6,1> vec;//6*1 vector:column vec;
	for(int i=0;i<3;i++){
	    vec(i)=sol.getPosition()[i];
	    vec(i+3)=sol.getvelocity()[i];
	}//Write the right-hand vector!
	vec=Jacob*vec;
	for(int i=0;i<3;i++){
	    vec(i)=vec(i)-H.getPosition()[i];
	    vec(i+3)=vec(i+3)-H.getvelocity()[i];
	}
	Matrix<double,6,1> Sol=Jacob.lu().solve(vec);//get the vector
	double pos[3];
	double vel[3];
	for(int i=0;i<3;i++){
	    pos[i]=Sol(i);
	    vel[i]=Sol(i+3);
	}
	sol.setPosition(pos);
	sol.setvelocity(vel);
	H=Hfunc(sol);
	looptime++;
    }
    return sol;
}
#else
//do nothing
#endif
