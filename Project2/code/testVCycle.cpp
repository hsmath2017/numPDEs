/**
 * @file   testVCycle.cpp
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Thu Jun  3 06:06:35 2021
 * 
 * @brief  test the function VCycle
 * 
 * 
 */
#include "VCycle.h"
#include <vector>
using namespace std;
int main(){
    VCycle V(1024,1,1);
    V.setBC(0,0);
    Meshtype1 M=V.getmesh();
    valarray<double> arr=V.generateRH(M);
    int n=1;
    vector<double> errors;
    for(n=1;n<=32;n++){
	for(int i=0;i<n;i++){
	    M=V.cycle(M,arr);
	    V.setmesh(M);
	}
	double err=V.generateerror();
	errors.push_back(err);
	if(n<32){
	M.setzero();
	}
    }
    ofstream os;
    os.open("VCycleerror1024.txt");
    os<<"\\hline"<<endl;
    os<<"cycletimes & error & ratio"<<endl;
    for(int i=0;i<32;i++){
	n=i+1;
	os<<n<<" & "<<errors[i];
	if(i>=1){
	    double ratio=errors[i-1]/errors[i];
	    os<<"& "<<ratio;
	}
	os<<"\\\\ "<<"\\hline";
	os<<endl;
    }
    os.close();
    V.getresult();
}
