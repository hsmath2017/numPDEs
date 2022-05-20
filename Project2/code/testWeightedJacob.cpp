/**
 * @file   testWeightedJacob.cpp
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Thu Jun  3 06:06:35 2021
 * 
 * @brief  test the function WeightedJacobian
 * 
 * 
 */
#include "VCycle.h"
using namespace std;
int main(){
    VCycle V(128,50000,50000);
    V.setBC(1,exp(sin(1)));
    Meshtype1 M=V.getmesh();
    valarray<double> arr=V.generateRH(M);
    M=V.WeightedJacob(20000,M,arr);
    V.setmesh(M);
    V.getresult();
}
