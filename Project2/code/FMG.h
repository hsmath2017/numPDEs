/**
 * @file   FMG.h
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Sat Jun  5 23:19:36 2021
 * 
 * @brief  Achieve FMG Algorithm (A better multigrid solver)
 * 
 * 
 */
#ifndef _HSMATH_FMG
#define _HSMATH_FMG
#include "VCycle.h"
using namespace std;
class FMG:public VCycle{
public:
    FMG(int sz,int mu1,int mu2):VCycle(sz,mu1,mu2){firstin=true;}
    Meshtype1 FMGcycle(valarray<double>& f,int vcycletimes);
    void setfirstin(){
	firstin=true;
    }
private:
    bool firstin;
};
Meshtype1 FMG::FMGcycle(valarray<double> &f,int vcycletimes){
    static Meshtype1 m(4);
    if(firstin==true){
	m.setsize(3);
	m.setzero();
    }
    firstin=false;
    if(f.size()==3){//exit!
	m=cycle(m,f);
	return m;
    }
    int len=f.size();
    valarray<double> g;
    g=f;
    Meshtype1 M1(len+1);
    M1.setinfo(g);
    M1.tocoarse();
    g=M1.getinfo();
    m=FMGcycle(g,vcycletimes);
    m.torefine();
    g=f;
    for(int i=0;i<vcycletimes;i++){
	m=cycle(m,g);
    }
    return m;
}
#else
//NOTHING
#endif
