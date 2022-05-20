#include "AdamMoulton.h"
int main(){
    double input[6]={0.87978,0,0,0,-0.3797,0};
    double T1=19.14045706162071;
    int N0[4]={1000,850,700,180};
    int N;
    double k=T1/N;
    double mu=1/81.45;
    int p=4;
    Point sol[3];
    for(int p=2;p<=5;p++){
	N=N0[p-2];
	k=T1/N;
	for(int i=0;i<3;i++){
	    AdamMoulton AM;
	    AM.SetInfo(input,k,T1,mu,p);
	    AM.solveODE();
	    sol[i]=AM.getSolution();
	    AM.PrintSolution();
	    N=N*2;
	    k=T1/N;
	}
	Point err1=sol[0]-sol[2];
	Point err2=sol[1]-sol[2];
	double precision=log(err1.norm()/err2.norm()-1)/log(2);
	cout<<precision<<endl;
    }
}
