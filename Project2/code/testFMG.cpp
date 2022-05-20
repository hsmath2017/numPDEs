/**
 * @file   testFMG.cpp
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Sun Jun  6 00:01:29 2021
 * 
 * @brief  Designed to test the file FMG.h
 * 
 * 
 */
#include "FMG.h"
int main(){
    FMG F(1024,1,1);
    F.setBC(1,exp(sin(1)));
    Meshtype1 M=F.getmesh();
    valarray<double> arr=F.generateRH(M);
    M=F.FMGcycle(arr,5);
    F.setmesh(M);
    F.getresult();
}
