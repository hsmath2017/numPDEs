#include "../src/MultiGrid/MGSolver.h"
#include "Functions.h"
template<typename T>
using Vector=std::vector<T>;   
template<int Dim> 
using VPR=Vector<Restrictor<Dim>*>;
template<int Dim>
using VPI=Vector<Interpolator<Dim>*>;
int main(){
    Vector<RectDomain<2>> GridVec;
    VPR<2> VecRestrict;
    VPI<2> VecInterp;
    for(int i=4;i>=2;i--){
        Box<2> bx({0,0},{(int)std::pow(2,i),(int)std::pow(2,i)});
        int gridnum=std::pow(2,i);
        Vec<Real,2> dx={1.0/gridnum,1.0/gridnum};
        RectDomain<2> RD(bx,dx,NodeCentered,0);
        GridVec.push_back(RD);
    }
    for(int i=0;i<GridVec.size()-1;i++){
        FullWeighting<2>* pRes=new FullWeighting<2>(GridVec[i],GridVec[i+1]);
        VecRestrict.push_back(pRes);
        LinearInterpolator<2>* pInterp=new LinearInterpolator<2>(GridVec[i],GridVec[i+1]);
        VecInterp.push_back(pInterp);
    }
    MGSolver<2> MGS(GridVec,VecRestrict,VecInterp);
    MGParam param;
    param.maxIter=0;
    param.numPreIter=1;
    param.numPostIter=1;
    param.numBottomIter=1;
    param.reltol=0.01;
    MGS.setParam(param);
    ScalarFunction<2>* RhsFunc=new RightHand;
    ScalarFunction<2>* BdryFunc=new BdryCond;
    Tensor<Real,2> phi;
    phi.resize(GridVec[0]);
    MGS.solve(phi,RhsFunc,BdryFunc,false);
    FuncFiller<2> FF(GridVec[0]);
    ScalarFunction<2>* SolFunc=new Sol;
    Tensor<Real,2> sol;
    sol.resize(phi.box());
    FF.fill(sol,SolFunc);
    std::cout<<phi<<std::endl;
    std::cout<<sol<<std::endl;
    Real maxErr=MGS.getError(phi,sol);
    std::cout<<maxErr<<std::endl;
}