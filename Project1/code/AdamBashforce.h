#ifndef _HSMATH_ADAMBASHFORCE
#define _HSMATH_ADAMBASHFORCE
#include "RungeKutta.h"
#include <map>
class AdamBashforce:public RungeKutta{
private:
    map<int,vector<double>> coefficient_table;
public:   
    AdamBashforce()
    {
        //Just initialize the base class
    }
    void read_table();//Read the coefficient table of vector \alpha and \beta
    void generate_initial();//Generate the initial_values
    void onestepiterator();//The one_step iterator for A-B method
    void solveODE();//The ode_solver for A-B method
    ~AdamBashforce(){
    }
};
void AdamBashforce::read_table(){
    string tmp;
    ifstream in("AdamBashforceTable");
    getline(in,tmp);//the head of table
    int flag=0;
    while(!in.eof()&&flag<p){
	int s,p,t;//s:step num,p:precision,t:beta_s,always zero
	in>>s>>p>>t;
	vector<double> vec;
	for(int i=0;i<s;i++){
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
void AdamBashforce::generate_initial(){
    int init_num=p;//Need the number of initial value!
    for(int i=0;i<p-1;i++){
	RungeKutta::onestepiterator();
    }//p_steps Forward-Euler method to generate the initial_value!
}
void AdamBashforce::onestepiterator(){//Go forward one-step
    Point Pnew;
    Pnew=Pts[Pts.size()-1];//The position of point
    for(int i=0;i<beta.size();i++){//go step_forward!
	Pnew=Pnew+k*beta[i]*func(Pts[Pts.size()-1-i],mu);
    }
    Pts.push_back(Pnew);//one-step forward NOW!
}
void AdamBashforce::solveODE(){//Solve the ODE
    read_table();
    generate_initial();
    for(int i=0;i<coefficient_table[p].size();i++){
	beta.push_back(coefficient_table[p][i]);//Input the beta vector
    }
    for(int i=p;i<T/k;i++){
	onestepiterator();
    }
}

#else
//do nothing
#endif
