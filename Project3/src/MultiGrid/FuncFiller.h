#ifndef _Func_Filler
#define _Func_Filler
#include "ScalarFunction.h"
#include "Core/RectDomain.h"
template<int Dim>
class FuncFiller{
protected:
    RectDomain<Dim> domain;
public:
    FuncFiller(const RectDomain<Dim>& adomain):domain(adomain){};
    //fill the inside tensor with the rhs function f.
    inline void fill(Tensor<Real,Dim>& res,const ScalarFunction<Dim>* pfunc) const;
};
template<int Dim>
void FuncFiller<Dim>::fill(Tensor<Real,Dim>& res,const ScalarFunction<Dim>* pfunc) const{
    Box<Dim> bx=res.box();
    assert(domain.contain(bx));
    assert(domain.getCentering()==NodeCentered);
    loop_box_2(bx,i,j){
        auto pt=domain.spacing()*Vec<Real,2>({i,j});
        res(i,j)=pfunc->operator()(pt);
    }
} 
#else
#endif