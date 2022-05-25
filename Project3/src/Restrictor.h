#ifndef _Restrictor
#define _Restrictor
#include "Core/RectDomain.h"
#include "Core/TensorExpr.h"
#include <vector>
//Achieve the restrict operator.
template<int Dim>
class Restrictor{
protected:
    RectDomain<Dim> coarseDomain;
    RectDomain<Dim> fineDomain;
public:
    Restrictor(const RectDomain<Dim>& D1,const RectDomain<Dim>& D2):coarseDomain(D1),fineDomain(D2){};
    virtual void apply(const Tensor<Real,Dim>& data,Tensor<Real,Dim>& res) const = 0;
};
template<int Dim>
class Injection:public Restrictor<Dim>{
public:
    Injection(const RectDomain<Dim>& acoarse,const RectDomain<Dim>& afine):Restrictor<Dim>(acoarse,afine){};
    void apply(const Tensor<Real,Dim>& data,Tensor<Real,Dim>& res) const{
        if(Dim==2){
            Box<Dim> bx_res=res.box();
            loop_box_2(bx_res,i,j){
                int nx=i*2;
                int ny=j*2;
                res(i,j)=data(nx,ny);
            }
        }
    } 
};
template<int Dim>
class FullWeighting:public Restrictor<Dim>{
public:
    FullWeighting(const RectDomain<Dim>& acoarse,const RectDomain<Dim>& afine):Restrictor<Dim>(acoarse,afine){};
    void apply(const Tensor<Real,Dim>& data,Tensor<Real,Dim>& res) const{
        if(Dim==2){
            Box<Dim> bx_res=res.box();
            loop_box_2(bx_res,i,j){
                std::vector<Vec<int,2>> Dir={{-1,0},{1,0},{0,-1},{0,1}};
                Vec<int,2> origin{i*2,j*2};
                bool isBdry=false;
                for(auto it:Dir){
                    auto tmp=origin+it;
                    if(data.box().contain(tmp)==false){
                        isBdry=true;
                        break;
                    }
                }
                if(isBdry){
                    res(i,j)=data(origin);
                }else{
                    Real val=0;
                    val+=data(origin)*1.0/2;
                    for(auto it:Dir){
                        auto tmp=origin+it;
                        val+=data(tmp)*1.0/8;
                    }
                    res(i,j)=val;
                }
            }
        }
    };
};
#else
//nothing
#endif