//Handle the laplacian operator
#ifndef _Laplacian
#define _Laplacian
#include "Smoother.h"
#include "Core/Tensor.h"
#include "Core/TensorExpr.h"
#include <vector>
template<int Dim>
class Laplacian{
protected:
    RectDomain<Dim> domain;
    Smoother<Dim>* psmoother;
public:
    Laplacian(const RectDomain<Dim>& adomain,Smoother<Dim>* apsmoother):domain(adomain),psmoother(apsmoother){};
    void apply(const Tensor<Real,Dim>& phi,Tensor<Real,Dim>& res) const;
    void smooth(const Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs, Tensor<Real,Dim>& res) const;
    void computeResidual(const Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs,Tensor<Real,Dim>& res) const;
    Vec<Real,2> getsize(){
        return domain.spacing();
    }
};
template<int Dim>
void Laplacian<Dim>::apply(const Tensor<Real,Dim>& phi,Tensor<Real,Dim>& res) const{
    Box<Dim> bx_phi=phi.box();
    Box<Dim> bx_res=res.box();
    assert(bx_phi.contain(bx_res)&&bx_res.contain(bx_phi));
    assert(domain.contain(bx_res));
    std::vector<Vec<int,Dim>> Directions;
    if(Dim==2){
        Directions={{1,0},{-1,0},{0,1},{0,-1}};
        loop_box_2(bx_res,i,j){
            std::vector<Vec<int,Dim>> neighbors;
            Vec<int,Dim> origin{i,j};
            for(auto it:Directions){
                auto tmp=origin+it;
                neighbors.push_back(tmp);
            }
            bool isBdry=false;
            for(auto it:neighbors){
                if(domain.contain(it)==false){
                    isBdry=true;
                    break;
                }
            }
            if(isBdry){
                res(origin)=phi(origin);
            }else{
                Real val=0;
                for(auto x:neighbors){
                    bool flag=true;
                    for(auto it:Directions){
                        auto tmp=it+x;
                        if(domain.contain(tmp)==false){
                            flag=false;
                            break;
                        }
                    }
                    if(flag){
                        val=val-phi(x);
                    }
                }
                val=val+4*phi(origin);
                val=val/prod(domain.spacing());
                res(origin)=val;
            }
        }
    }
};

template<int Dim>
void Laplacian<Dim>::smooth(const Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs, Tensor<Real,Dim>& res) const{
    psmoother->apply(phi,rhs,res);
};

template<int Dim>
void Laplacian<Dim>::computeResidual(const Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs,Tensor<Real,Dim>& res) const{
    Box<Dim> bx_phi=phi.box();
    Box<Dim> bx_rhs=rhs.box();
    Box<Dim> bx_res=res.box();
    assert(bx_phi.contain(bx_rhs)&&bx_rhs.contain(bx_phi));
    assert(bx_phi.contain(bx_res)&&bx_res.contain(bx_phi));
    Tensor<Real,Dim> LapPhi;
    LapPhi.resize(phi.box());
    apply(phi,LapPhi);
    Box<2> bx=LapPhi.box();
    bx=bx.inflate(-1);
    std::cout<<"LapPhi = "<<LapPhi<<std::endl;
    loop_box_2(bx,i,j){
        res(i,j)=rhs(i,j)-LapPhi(i,j);
    }
};
#else
//nothing
#endif
