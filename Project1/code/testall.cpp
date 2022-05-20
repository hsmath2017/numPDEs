#include "TimeIntegratorFactory.h"
int main(){
    int N=100000;
    int p=1;//precision
    int i;
    double input[6]={0.994,0,0,0,-2.0015851063790825224,0};
    double T1=17.06521656015796;
    double k=T1/N;
    double mu=1/81.45;
    for(i=1;i<=3;i++){
	for(p=2;p<=5;p++){
	    string Method="AdamMoulton";
	    TimeIntegrator* pIntegrator = TimeIntegratorFactory::Instance()
		->CreateTimeIntegrator(Method);
	    pIntegrator->SetInfo(input,k,T1,mu,p);
	    pIntegrator->solveODE();
	    pIntegrator->PrintSolution();
	    cout<<"Relative Error:"<<pIntegrator->calculateRelativeError()<<endl;
	}
	N=N*2;
	k=T1/N;
    }
}
