#ifndef _HSMATH_BDF
#define _HSMATH_BDF
#include <Eigen/LU>
#include "RungeKutta.h"
#include <map>
class BDF:public RungeKutta{
private:
    map<int,vector<double>> coefficient_table;
public:   
    BDF(){
        //Just initialize the base class
    }
    void read_table();//Read the coefficient table of vector \alpha and \beta
    void generate_initial();//Generate the initial_values
    void onestepiterator();//The one_step iterator for A-B method
    void solveODE();//The ode_solver for A-B method
    Point Newton_Iterator();//Newton's Iterator method to solve unlinear equation
    Point Hfunc(Point& P);//Iterator function
    ~BDF(){
    }
};
Point BDF::Hfunc(Point &P){
    Point H;
    H=P-k*beta[0]*func(P,mu);
    for(int i=1;i<=p;i++){
	H=H+alpha[i-1]*Pts[Pts.size()-i];
    }
    return H;
}
void BDF::read_table(){//read the coefficient table,similar with A-B method!
    string tmp;
    ifstream in("BDFTable");
    getline(in,tmp);//the head of table
    int flag=0;
    while(!in.eof()&&flag<p){
	int s,p;//s:step num,p:precision,t:beta_s,always zero
	in>>s>>p;
	vector<double> vec;
	for(int i=0;i<p+1;i++){//read the value of alpha and beta_s
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
void BDF::generate_initial(){
    int init_num=p;//Need the number of initial value!
    for(int i=0;i<p-1;i++){
        RungeKutta::onestepiterator();
    }//s_steps Runge_Kutta method to generate the initial_value!
}
void BDF::onestepiterator(){//Go forward one-step
    Point Pnew;
    Pnew=Newton_Iterator();
    Pts.push_back(Pnew);//one-step forward NOW!
}
void BDF::solveODE(){//Solve the ODE
    read_table();
    generate_initial();
    for(int i=0;i<coefficient_table[p].size()-1;i++){
	alpha.push_back(coefficient_table[p][i]);//Input the alpha vector
    }
    beta.push_back(coefficient_table[p][coefficient_table[p].size()-1]);//beta value
    for(int i=p;i<T/k;i++){
	onestepiterator();
    }
}
Point BDF::Newton_Iterator(){
    Matrix<double,6,6> Jacob;
    Matrix<double,6,6> I;//Identity Matrix
    I.setIdentity(6,6);//identity matrix
    Point sol=Pts[Pts.size()-1];//Guess:The solution is near the point "sol"
    Point H=Hfunc(sol);
    int looptime=0;
    while(H.norm()>1.0e-8&&looptime<=4){//Newton Iterator!
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
