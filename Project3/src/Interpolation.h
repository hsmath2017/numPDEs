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
    virtual void apply(const Tensor<Real,Dim>& data,Tensor<Real,Dim>& res) = 0;
};
template<int Dim>
class LinearInterpolator:public Interpolator<Dim>{
public:
    void apply(const Tensor<Real,Dim>& data,Tensor<Real,Dim>& res){
        if(Dim==2){
            loop_box_2(res.box(),i,j){
                if(i%2==0&&j%2==0){
                    res(i,j)=data(i/2,j/2);
                }else if(i%2==0){
                    res(i,j)=(data(i/2,j/2)+data(i/2,j/2+1))*1.0/2;
                }else if(j%2==0){
                    res(i,j)=(data(i/2,j/2)+data(i/2+1,j/2))*1.0/2;
                }else{
                    res(i,j)=(data(i/2,j/2)+data(i/2,j/2+1)+data(i/2+1,j/2)+data(i/2+1,j/2+1))*1.0/4;
                }
            }
        }
    }
};
#else
#endif