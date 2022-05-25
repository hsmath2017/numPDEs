#ifndef _MGSolver
#define _MGSolver
#include "FuncFiller.h"
#include "BoundaryFiller.h"
#include "Interpolation.h"
#include "Restrictor.h"
#include "Laplacian.h"

template<int Dim>
class MGSolver
{
public:
    struct MGParam{
        int numPreIter;
        int numPostIter;
        int numBottomIter;
        Real reltol;
        int maxIter;
    };
    template<class T>
    using Vector=std::vector<T>;
    using VPR=Vector<Restrictor<Dim>*>;
    using VPI=Vector<Interpolator<Dim>*>;
protected:
    Vector<RectDomain<Dim>> vDomain;
    Vector<BoundaryFiller<Dim>> vBdryFiller;
    Vector<Laplacian<Dim>> vLaplacian;
    VPR vpRestriction;
    VPI vpInterpolation;
    PoissonDirectSolver<Dim> bottomSolver;
    MGParam param;
public:
    MGSolver(const Vector<RectDomain<Dim>>& avDomain,const VPR& avpRestriction,const VPI& avpInterpolation);

    void setParam(const MGParam& aparam);

    void VCycle(int depth,Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs) const;

    void FMVCycle(int depth,Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs) const;

    void solve(Tensor<Real,Dim>& phi,ScalarFunction<Dim>* RightFunc,ScalarFunction<Dim>* BdryFunc,bool useFMVCycle=0) const; 
};

#else
//nothing
#endif