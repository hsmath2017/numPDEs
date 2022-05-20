/**
 * @file   VCycle.h
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Thu Jun  3 04:56:54 2021
 * 
 * @brief  Class VCycle.h:Achieve VCycle algorithm.
 * 
 * 
 */
#ifndef _HSMATH_VCYCLE
#define _HSMATH_VCYCLE
#include "Meshtype1.h"
#include "Function.h"
#include<fstream>
using namespace std;
class VCycle{
protected:
    Meshtype1 M;
    int mu1;//Relax times on 1st guess
    int mu2;//Relax times on correlated guess
    double LeftBC;
    double RightBC;
    valarray<double> Righthand;
    double maxerror;
public:
    VCycle(int sz,int mu1,int mu2):M(sz)
    {
	this->mu1=mu1;
	this->mu2=mu2;
	Righthand=generateRH(M);
    }
    void seterr(double e){
	maxerror=e;
    }
    void setmesh(Meshtype1 m){
	M=m;
    }
    Meshtype1 getmesh(){
	return M;
    }
    Meshtype1 WeightedJacob(int mu,Meshtype1& M,valarray<double>& Righthand);//Weighted-Jacob Iteration for mu times.
    void setBC(double L,double R){
	LeftBC=L;
	RightBC=R;
    }
    valarray<double> generateRH(Meshtype1& M){
	int sz=M.getsize();
	int n=sz+1;
	valarray<double> Righthand(sz);//Righthand vector 'f' denoted by function 'func'
	for(int i=0;i<sz;i++){
	    double d=0;
	    d=func(1.0*(i+1)/n)/(n*n);
	    Righthand[i]=d;
	}
	Righthand[0]=Righthand[0]+LeftBC;
	Righthand[sz-1]=Righthand[sz-1]+RightBC;
	return Righthand;
    }
    valarray<double> Residure(valarray<double>& f1,valarray<double>& v1);
    Meshtype1 cycle(Meshtype1& M,valarray<double>& f);//V-Cycle Recursive
    void getresult();
    double generateerror(){
	valarray<double> val=M.getinfo();
	valarray<double> f=generateRH(M);
	valarray<double> res=Residure(f,val);
	double maxnorm1=0;
	double maxnorm2=0;
	for(int i=0;i<res.size();i++){
	    if(abs(res[i])>maxnorm1){
		maxnorm1=abs(res[i]);
	    }
	    if(abs(f[i])>maxnorm2){
		maxnorm2=abs(f[i]);
	    }
	}
	double rel=maxnorm1/maxnorm2;
	return rel;
    }
    bool quit(){//quit the VCycle or not!
	valarray<double> val=M.getinfo();
	valarray<double> f=generateRH(M);
	valarray<double> res=Residure(f,val);
	double maxnorm1=0;
	double maxnorm2=0;
	for(int i=0;i<res.size();i++){
	    if(abs(res[i])>maxnorm1){
		maxnorm1=abs(res[i]);
	    }
	    if(abs(f[i])>maxnorm2){
		maxnorm2=abs(f[i]);
	    }
	}
	double rel=maxnorm1/maxnorm2;
	if(rel<maxerror){
	    return true;
	}else{
	    return false;
	}
    }
};
Meshtype1 VCycle::WeightedJacob(int mu,Meshtype1& M,valarray<double>& Righthand){
    int sz=M.getsize();
    int n=sz+1;
    for(int i=0;i<mu;i++){
	valarray<double> Oldone=M.getinfo();
	valarray<double> c=0.5*Righthand;//The right-hand vector c in Weighted-Jacobi Mode (1/2*f)
	valarray<double> vec(sz);//-D^{-1}(L+U)*u
	vec[0]=0.5*Oldone[1];
	vec[sz-1]=0.5*Oldone[sz-2];
	for(int i=1;i<sz-1;i++){
	    vec[i]=0.5*Oldone[i-1]+0.5*Oldone[i+1];
	}
	valarray<double> ustar(sz);
	ustar=vec+c;
	valarray<double> onetimeiterator(sz);
	onetimeiterator=(2.0/3)*ustar+(1.0/3)*Oldone;
	M.setinfo(onetimeiterator);
    }
    Meshtype1 m(sz);
    m=M;
    return m;
}
void VCycle::getresult(){
    valarray<double> v=M.getinfo();
    ofstream os;
    int sz=M.getsize();
    int n=sz+1;
    os.open("Result.m");
    os<<"x=[0:1/"<<n<<":1];\n";
    os<<"y=["<<LeftBC<<",";
    for(int i=0;i<sz;i++){
	os<<v[i]<<",";
    }
    os<<RightBC<<"];\n";
    os<<"plot(x,y)\n";
    os.close();
}
valarray<double> VCycle::Residure(valarray<double>& f,valarray<double>& v){ 
    int len=v.size();
    valarray<double> Av(len);
    Av[0]=2*v[0]-v[1];
    for(int i=1;i<len-1;i++){
	Av[i]=2*v[i]-v[i-1]-v[i+1];
    }
    Av[len-1]=2*v[len-1]-v[len-2];
    valarray<double> right=f-Av;//grid
    return right;
}
Meshtype1 VCycle::cycle(Meshtype1& M,valarray<double>& f){
    /*exit*/
    if(M.getsize()==3){
	M=WeightedJacob(50,M,f);
	return M;
    }
    Meshtype1 M1=WeightedJacob(mu1,M,f);//First of all:relax mu1 times
    valarray<double> v=M1.getinfo();
    int len=v.size();
    valarray<double> right=Residure(f,v);
    right=(4.0)*right;
    Meshtype1 M2(len+1);
    M2.setinfo(right);
    M2.tocoarse();
    valarray<double> g(len+1);
    g=M2.getinfo();
    M2.setzero();
    M2=cycle(M2,g);
    M2.torefine();
    M1=M1+M2;
    M1=WeightedJacob(mu2,M1,f);
    valarray<double> inf=M1.getinfo();
    valarray<double> Res=Residure(f,inf);
    return M1;
}
#else
//NOTHING
#endif
