#ifndef _Test_Function
#define _Test_Function
#include "MultiGrid/ScalarFunction.h"
class Sol:public ScalarFunction<2>{
public:
    const Real operator()(const Vec<Real,2>& pt) const{
        Real res=exp(pt[0]+pt[1])*cos(pt[0])*sin(pt[1]);
        return res;
    };
};
class RightHand:public ScalarFunction<2>{
public:
    const Real operator()(const Vec<Real,2>& pt) const{
        Real res=-2*exp(pt[0]+pt[1])*cos(pt[0]+pt[1]);
        return res;
    }
};

class BdryCond:public ScalarFunction<2>{
public:
    const Real operator()(const Vec<Real,2>& pt) const{
        Real res=exp(pt[0]+pt[1])*cos(pt[0])*sin(pt[1]);
        return res;
    }
};
#else
#endif