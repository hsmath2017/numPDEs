#include "AdamBashforce.h"
int main(){
    double input[6]={0.87978,0,0,0,-0.3797,0};
    double T1=19.14045706162071;
    int N0[4]={80000,16000,8000,2000};
    int N;
    double k=T1/N;
    double mu=1/81.45;
    int p=4;
    Point sol[3];
    for(int p=1;p<=4;p++){
	N=N0[p-1];
	k=T1/N;
	for(int i=0;i<3;i++){
	    AdamBashforce AB;
	    AB.SetInfo(input,k,T1,mu,p);
	    AB.solveODE();
	    sol[i]=AB.getSolution();
	    AB.PrintSolution();
	    N=N*2;
	    k=T1/N;
	}
	Point err1=sol[0]-sol[2];
	Point err2=sol[1]-sol[2];
	double precision=log(err1.norm()/err2.norm()-1)/log(2);
	cout<<precision<<endl;
    }
}
