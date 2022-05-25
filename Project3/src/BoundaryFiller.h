#ifndef _Boundary_Filler
#define _Boundary_Filler
#include "ScalarFunction.h"
#include "Core/RectDomain.h"
template<int Dim>
class BoundaryFiller{
protected:
    RectDomain<Dim> domain;
public:
    BoundaryFiller(const RectDomain<Dim>& adomain):domain(adomain){};
    //Mark the normal vector of the streight edges
    //Left/Down:Mark as low,Right/Up:Mark as high
    enum normal{low=-1,high=1};
    //Fill the value of boundary function on one side.
    void fillOneSide(Tensor<Real,Dim>& res,int dim,normal N,const ScalarFunction<Dim>* pfunc=nullptr) const;
    //Fill ALL sides.
    void fillAllSides(Tensor<Real,Dim>& res,ScalarFunction<Dim>* pfunc=nullptr) const;
};
template<int Dim>
void BoundaryFiller<Dim>::fillOneSide(Tensor<Real,Dim>& res,int dim,normal N,const ScalarFunction<Dim>* pfunc=nullptr) const{
    Vec<int,Dim> start;
    Vec<int,Dim> end;
    Vec<int,Dim> info={dim,N};
    if(Dim==2){
    switch (info)
    {
    case {0,-1}:
        start=domain.lo();
        end=start;
        end[0]=domain.hi()[0];
        break;
    case {0,1}:
        end=domain.hi();
        start=end;
        start[0]=domain.lo()[0];
        break;
    case {1,-1}:
        start=domain.lo();
        end=start;
        end[1]=domain.hi()[1];
        break;
    case {1,1}:
        end=domain.hi();
        start=end;
        start[1]=domain.lo()[1];
        break;
    }
    }
    Box<2> side(start,end);
    loop_box_2(side,i,j){
        auto pt=domain.spacing()*Vec<Real,2>{i,j};
        res(i,j)=pfunc->operator()(pt);
    }
};
template<int Dim>
void BoundaryFiller<Dim>::fillAllSides(Tensor<Real,Dim>& res,ScalarFunction<Dim>* pfunc=nullptr) const{
    auto bx=res.box();
    assert(domain.contain(bx));
    if(Dim==2){
        for(dim=0;dim<2;dim++){
            fillOneSide(res,dim,low,pfunc);
            fillOneSide(res,dim,high,pfunc);
        }
    }
};
#else
#endif