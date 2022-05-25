#ifndef _Smoother
#define _Smoother
#include "Core/RectDomain.h"
#include "Core/Tensor.h"
template<int Dim>
class Smoother{
protected:
    RectDomain<Dim> domain;
public:
    //smooth for one time.
    Smoother(const RectDomain<Dim>& adomain):domain(adomain){};
    virtual void apply(const Tensor<Real,Dim>& phi, const Tensor<Real,Dim>& rhs, Tensor<Real,Dim>& res) const = 0;
};
#else
//nothing
#endif