#ifndef _Test_MG
#define _Test_MG
#include "MGFactory.h"
template<int Dim>
class TestMultigrid{
protected:
    std::unique_ptr<MGSolver<Dim>> pMG;
    std::string multigridID;
public:
    TestMultigrid(const std::string& jsonfile);
    Tensor<Real,Dim> solve(const ScalarFunction<Dim>* pfunc,const ScalarFunction<Dim>* BdryFunc,bool useFMVCycle=0) const;
    Real computeError(const Tensor<Real,Dim>& res,const ScalarFunction<Dim>* pfunc,const int p=0) const;
    void plot(const Tensor<Real,Dim>& res,const std::string &file) const;
    void test(const ScalarFunction<Dim>* pfunc,const ScalarFunction<Dim>* BdryFunc,const int numEncryp,bool useFMVCycle=0) const;
};
#else
//nothing
#endif