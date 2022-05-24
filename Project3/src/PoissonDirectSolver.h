#ifndef _Direct_Solver
#define _Direct_Solver
#include "Core/RectDomain.h"
#include "Core/Tensor.h"
#include "Core/TensorExpr.h"
template<int Dim>
class PoissonDirectSolver{
public:
    void solve(Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs);
};
#else
//nothing
#endif