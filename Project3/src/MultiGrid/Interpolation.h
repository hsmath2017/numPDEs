#ifndef _Interpolator
#define _Interpolator
#include "Core/RectDomain.h"
#include "Core/TensorExpr.h"
#include <vector>
template<int Dim>
class Interpolator{
protected:
    RectDomain<Dim> coarseDomain;
    RectDomain<Dim> fineDomain;
public:
    Interpolator(const RectDomain<Dim>& acoarse,const RectDomain<Dim>& afine):coarseDomain(acoarse),fineDomain(afine){};

    virtual void apply(const Tensor<Real,Dim>& data,Tensor<Real,Dim>& res) = 0;
};
template<int Dim>
class LinearInterpolator:public Interpolator<Dim>{
public:
    LinearInterpolator(const RectDomain<Dim>& acoarse,const RectDomain<Dim>& afine):Interpolator<Dim>(acoarse,afine){};

    void apply(const Tensor<Real,Dim>& data,Tensor<Real,Dim>& res){
        if(Dim==2){
            Tensor<Real,Dim> tmp;
            Box<Dim> bx=data.box();
            tmp.resize(bx);
            loop_box_2(bx,i,j){
                tmp(i,j)=0;
            }
            bx=bx.inflate(-1);
            loop_box_2(bx,i,j){
                tmp(i,j)=data(i,j);
            }
            loop_box_2(res.box(),i,j){
                if(i%2==0&&j%2==0){
                    res(i,j)=tmp(i/2,j/2);
                }else if(i%2==0){
                    res(i,j)=(tmp(i/2,j/2)+tmp(i/2,j/2+1))*1.0/2;
                }else if(j%2==0){
                    res(i,j)=(tmp(i/2,j/2)+tmp(i/2+1,j/2))*1.0/2;
                }else{
                    res(i,j)=(tmp(i/2,j/2)+tmp(i/2,j/2+1)+tmp(i/2+1,j/2)+tmp(i/2+1,j/2+1))*1.0/4;
                }
            }
        }
    }
};
#else
#endif