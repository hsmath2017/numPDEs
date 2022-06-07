#ifndef _Scalar_Func
#define _Scalar_Func
#include "Core/Tensor.h"
template<int Dim>
class ScalarFunction{
public:
    virtual const Real operator()(const Vec<Real,Dim>& pt) const =0;
};
#else
#endif