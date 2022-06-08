#ifndef _Weighted_Jacob
#define _Weighted_Jacob
#include "Smoother.h"
template<int Dim>
class WeightedJacobi:public Smoother<Dim>{
private:
    const Real weight=2.0/3; //The weight of the smoother
public:
    WeightedJacobi(const RectDomain<Dim>& adomain):Smoother<Dim>(adomain){};
    void apply(const Tensor<Real,Dim>& phi, const Tensor<Real,Dim>& rhs, Tensor<Real,Dim>& res) const;
};
template<int Dim>
void WeightedJacobi<Dim>::apply(const Tensor<Real,Dim>& phi, const Tensor<Real,Dim>& rhs, Tensor<Real,Dim>& res) const{
    if(Dim==2){
        auto dx=(this->domain).spacing();
        Real area=prod(dx);
        Tensor<Real,Dim> tmp;
        Box<Dim> bx=phi.box();
        Box<Dim> bx_rhs=rhs.box();
        Box<Dim> bx_res=res.box();
        assert(bx.contain(bx_res)&&bx_res.contain(bx));
        assert(bx.contain(bx_rhs)&&bx_rhs.contain(bx));
        tmp.resize(bx);
        loop_box_2(bx,i,j){
            tmp(i,j)=rhs(i,j);
        }
        auto newbx=bx.inflate(-1);
        loop_box_2(newbx,i,j){
            tmp(i,j)=0;//initialize
        }
        loop_box_2(newbx,i,j){
            tmp(i,j)+=area*rhs(i,j)/4.0;
            tmp(i,j)=tmp(i,j)+(1/4.0)*(phi(i-1,j)+phi(i+1,j)+phi(i,j-1)+phi(i,j+1));
        }
        res=tmp*weight+phi*(1-weight);
    }
};
#else
//nothing
#endif