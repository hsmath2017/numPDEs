#include "MGSolver.h"
template<int Dim>
MGSolver<Dim>::MGSolver(const Vector<RectDomain<Dim>>& avDomain,const VPR& avpRestriction,const VPI& avpInterpolation):vDomain(avDomain),vpRestriction(avpRestriction),vpInterpolation(avpInterpolation){
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
    Laplacian<Dim> LevelOp=Laplacian[depth];
    Tensor<Real,Dim> res=rhs;
    for(int i=0;i<param.numPreIter;i++){
        LevelOp.smooth(phi,rhs,res);
        phi=res;
    }
    if(depth==param.maxIter){
        //...
    }else{
        LevelOp.computeResidual(phi,rhs,res);
        Tensor<Real,Dim> newrhs;
        auto bx=res.box();
        Box<Dim> newbx;
        newbx.hi()=bx.hi()/2;
        newbx.lo()=bx.lo();
        newrhs.resize(newbx);
        vpRestriction[depth]->apply(res,newrhs);
        Tensor<Real,Dim> newres;
        newres.resize(newbx);
        VCycle(depth+1,newres,newrhs);
        Tensor<Real,Dim> residual;
        residual.resize(bx);
        vpInterpolation[depth]->apply(newres,residual);
        phi=phi+residual;
    }
    for(int i=0;i<param.numPostIter;i++){
        LevelOp.smooth(phi,rhs,res);
        phi=res;
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
        Box<Dim> newbx;
        newbx.lo()=bx.lo()/2;
        newbx.hi()=bx.hi()/2;
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
void MGSolver<Dim>::solve(Tensor<Real,Dim>& phi,ScalarFunction<Dim>* RightFunc,ScalarFunction<Dim>* BdryFunc,bool useFMVCycle=0) const{
    //Generate RHS
    Tensor<Real,Dim> rhs;
    Box<Dim> bx=vDomain[0];
    rhs.resize(bx);
    FuncFiller<Dim> FF(vDomain[0]);
    FF.fill(rhs,RightFunc);
    BoundaryFiller<Dim> BF(vDomain[0]);
    BF.fillAllSides(rhs,BdryFunc);
    //Solve
    if(useFMVCycle){
        FMVCycle(0,phi,rhs);
    }else{
        VCycle(0,phi,rhs);
    }
    //Generate Error
};

//==================================================
template class MGSolver<2>;