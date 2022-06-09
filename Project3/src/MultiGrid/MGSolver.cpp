#include "MGSolver.h"
template<int Dim>
MGSolver<Dim>::MGSolver(const Vector<RectDomain<Dim>>& avDomain,const VPR& avpRestriction,const VPI& avpInterpolation):vDomain(avDomain),vpRestriction(avpRestriction),vpInterpolation(avpInterpolation){
    for(auto d:avDomain){
        Smoother<Dim>* psmoother=new WeightedJacobi<Dim>(d);
        Laplacian<Dim> L(d,psmoother);
        vLaplacian.push_back(L);
    }
};

template<int Dim>
void MGSolver<Dim>::setParam(const MGParam& aparam){
    param.maxIter=aparam.maxIter;
    param.numBottomIter=aparam.numBottomIter;
    param.numPostIter=aparam.numPostIter;
    param.numPreIter=aparam.numPreIter;
    param.reltol=aparam.reltol;
};

template<int Dim>
void MGSolver<Dim>::VCycle(int depth,Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs) const{
    std::cout<<"Depth = "<<depth<<std::endl;
    Laplacian<Dim> LevelOp=vLaplacian[depth];
    std::cout<<"size = "<<LevelOp.getsize()<<std::endl;
    Tensor<Real,Dim> res;
    res.resize(rhs.box());
    // loop_box_2(rhs.box(),i,j){
    //     res(i,j)=rhs(i,j);
    // }
    for(int i=0;i<param.numPreIter;i++){
        LevelOp.smooth(phi,rhs,res);
        loop_box_2(rhs.box(),i,j){
            phi(i,j)=res(i,j);
        }
    }
    if(depth==param.maxIter){
        for(int i=0;i<param.numBottomIter;i++){
            LevelOp.smooth(phi,rhs,res);
            loop_box_2(rhs.box(),i,j){
                phi(i,j)=res(i,j);
            }
        }
    }else{
        std::cout<<"phi = "<<phi<<std::endl;
        LevelOp.computeResidual(phi,rhs,res);
        std::cout<<"res = "<<res<<std::endl;
        Tensor<Real,Dim> newrhs;
        Box<2> bx=res.box();
        auto vechi=bx.hi()/2;
        auto veclo=bx.lo();
        Box<Dim> newbx(veclo,vechi);
        newrhs.resize(newbx);
        vpRestriction[depth]->apply(res,newrhs);
        std::cout<<"newrhs = "<<newrhs<<std::endl;
        Tensor<Real,Dim> newres;
        newres.resize(newbx);
        VCycle(depth+1,newres,newrhs);
        std::cout<<"newres = "<<newres<<std::endl;
        Tensor<Real,Dim> residual;
        residual.resize(bx);
        vpInterpolation[depth]->apply(newres,residual);
        std::cout<<"residual = "<<residual<<std::endl;
        Box<2> parasbx=bx.inflate(-1);
        loop_box_2(parasbx,i,j){
            phi(i,j)=phi(i,j)+residual(i,j);
        }
    }
    std::cout<<"phi = "<<phi<<std::endl;
    for(int i=0;i<param.numPostIter;i++){
        LevelOp.smooth(phi,rhs,res);
        loop_box_2(rhs.box(),i,j){
            phi(i,j)=res(i,j);
        }
    }
    return;
}

template<int Dim>
void MGSolver<Dim>::FMVCycle(int depth,Tensor<Real,Dim>& phi,const Tensor<Real,Dim>& rhs) const{
    if(depth==param.maxIter){
        phi=0;
    }else{
        Tensor<Real,Dim> newrhs,newphi;
        Box<Dim> bx=rhs.box();
        auto veclo=bx.lo()/2;
        auto vechi=bx.hi()/2;
        Box<Dim> newbx(veclo,vechi);
        newrhs.resize(newbx);
        newphi.resize(newbx);
        vpRestriction[depth]->apply(rhs,newrhs);
        FMVCycle(depth+1,newphi,newrhs);
        vpInterpolation[depth]->apply(newphi,phi);
    }
    VCycle(depth,phi,rhs);
}; 

//RightFunc f:-\Delta u=f, BdryFunc g:u|_{\partial\Omega}=g.
template<int Dim>
void MGSolver<Dim>::solve(Tensor<Real,Dim>& phi,ScalarFunction<Dim>* RightFunc,ScalarFunction<Dim>* BdryFunc,bool useFMVCycle) const{
    //Generate RHS
    std::vector<Vec<int,2>> Dir={{1,0},{-1,0},{0,1},{0,-1}};
    Tensor<Real,Dim> rhs;
    Box<Dim> bx=vDomain[0];
    rhs.resize(bx);
    bx=bx.inflate(-1);
    Tensor<Real,Dim> tmp(bx);
    FuncFiller<Dim> FF(vDomain[0]);
    FF.fill(tmp,RightFunc);
    loop_box_2(bx,i,j){
        rhs(i,j)=tmp(i,j);
    }
    BoundaryFiller<Dim> BF(vDomain[0]);
    //BF.fillAllSides(rhs,BdryFunc);
    Real area=prod(vDomain[0].spacing());
    //remote the RHS on the boundary of domain.
    bx=bx.inflate(1);
    Box<Dim> innerbx=bx.inflate(-1);
    loop_box_2(innerbx,i,j){
        Vec<int,2> origin{i,j};
        for(auto c:Dir){
            auto tmp=origin+c;
            bool flag=false;//is 'tmp' a boundary cell?
            for(auto d:Dir){
                auto test=tmp+d;
                if(bx.contain(test)==false){
                    flag=true;
                    break;
                }
            }
            if(flag){
                Vec<Real,2> dx=vDomain[0].spacing();
                dx[0]=dx[0]*tmp[0];
                dx[1]=dx[1]*tmp[1];
                rhs(origin)+=BdryFunc->operator()(dx)/area;
            }
        }
    }
    //BF.fillAllSides(phi,BdryFunc);
    //Solve
    if(useFMVCycle){
        FMVCycle(0,phi,rhs);
    }else{
        for(int i=0;i<100;i++){
        VCycle(0,phi,rhs);
        }
    }
    BF.fillAllSides(phi,BdryFunc);
    //Generate Error
};

//==================================================
template class MGSolver<2>;