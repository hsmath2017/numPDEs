#ifndef _MGSolver
#define _MGSolver
#include "FuncFiller.h"
#include "BoundaryFiller.h"
#include "Interpolation.h"
#include "Restrictor.h"
#include "Laplacian.h"
#include "WeightedJacobi.h"
struct MGParam{
    int numPreIter;
    int numPostIter;
    int numBottomIter;
    Real reltol;
    int maxIter;
};
template<int Dim>
class MGSolver
{
public:
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
    MGParam param;
public:
    MGSolver(const Vector<RectDomain<Dim>>& avDomain,const VPR& avpRestriction,const VPI& avpInterpolation);

    void setParam(const MGParam& aparam);

    void VCycle(int depth,Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs) const;

    void FMVCycle(int depth,Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs) const;

    void solve(Tensor<Real,Dim>& phi,ScalarFunction<Dim>* RightFunc,ScalarFunction<Dim>* BdryFunc,bool useFMVCycle) const; 

    Real getError(const Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& sol){
        Real maxError=0;
        Box<Dim> bx=phi.box();
        loop_box_2(bx,i,j){
            if(std::abs(phi(i,j)-sol(i,j))>maxError){
                maxError+=std::abs(phi(i,j)-sol(i,j));
            }
        }
        auto sz=bx.size();
        maxError=maxError/(sz[0]*sz[1]);
        return maxError;
    }
};

#else
//nothing
#endif