#include "RungeKutta.h"
int main(){
    double input[6]={0.87978,0,0,0,-0.3797,0};
    double T1=19.14045706162071;
    int N=100;
    double k=T1/N;
    double mu=1/81.45;
    int p=4;
    Point sol[3];
    for(int i=0;i<3;i++){
        RungeKutta RK;
        RK.SetInfo(input,k,T1,mu,p);
        RK.solveODE();
        sol[i]=RK.getSolution();
        RK.PrintSolution();
        N=N*2;
        k=T1/N;
    }
    Point err1=sol[0]-sol[2];
    Point err2=sol[1]-sol[2];
    double precision=log(err1.norm()/err2.norm()-1)/log(2);
    cout<<precision<<endl;
}
