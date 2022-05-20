#include "TimeIntegratorFactory.h"
#include <ctime>
#include <chrono>
#include <ratio>
using namespace chrono;
int main(){
    ifstream in("InputFile");
    string tmp;
    getline(in,tmp);
    string method;
    int order;
    double Te;
    double input[6];
    double mu=1/81.45;
    int stepnum;
    in>>method>>order>>Te;
    for(int i=0;i<6;i++){
	in>>input[i];
    }
    in>>stepnum;
    double k=Te/stepnum;
    auto start = system_clock::now();
    TimeIntegrator* pIntegrator=TimeIntegratorFactory::Instance()->
	CreateTimeIntegrator(method);
    pIntegrator->SetInfo(input,k,Te,mu,order);
    pIntegrator->solveODE();
    auto end = system_clock::now();
    auto duration = duration_cast<microseconds>(end-start);
    pIntegrator->PrintSolution();
    pIntegrator->DrawFigure();
    cout<<pIntegrator->calculateMaxNormError()<<endl;
    cout<<"Time:"<<double(duration.count())*microseconds::period::num/microseconds::period::den<<" seconds."<<endl;
}
