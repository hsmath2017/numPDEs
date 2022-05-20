/**
 * @file   main.cpp
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Tue Jun  8 22:39:28 2021
 * 
 * @brief  Main function:Solve the 1D BVP
 * 
 * 
 */
#include "Mesh.h"
#include "Meshtype1.h"
#include "VCycle.h"
#include "FMG.h"
int main(){
    cout<<"Press the type of cycle, 1 means VCycle and 2 means FMG:";
    int type;
    cin>>type;
    cout<<"Press the boundary condition,format:Left Right:";
    double LBC,RBC;
    cin>>LBC>>RBC;
    cout<<"Press the type of stopping criteria, 1 means iteration time, 2 means relative error:";
    int stoptype;
    cin>>stoptype;
    cout<<"Press the number of grids:";
    int n;
    cin>>n;
    double relative_err=0;
    int maxloops=0;
    if(stoptype==2){
    cout<<"Press the maximum relative error:";
    cin>>relative_err;
    }else{
	cout<<"Press the maximum number of loops";
	cin>>maxloops;
    }
    //Default initial_guess:zeros!
    int looptime=1;
    VCycle VC(n,1,1);
    VC.setBC(LBC,RBC);
    VC.seterr(relative_err);
    FMG FC(n,1,1);
    FC.setBC(LBC,RBC);
    FC.seterr(relative_err);
    if(type==1){
	while(1){
	     if(VC.quit()==true&&stoptype==2){
        	break;
	     }else{
        	Meshtype1 M=VC.getmesh();
	 	valarray<double> arr=VC.generateRH(M);
                M=VC.cycle(M,arr);
	 	VC.setmesh(M);
	 	looptime++;
		if(looptime>maxloops&&stoptype==1){
		    break;
		}
	     }
	     if(looptime>=100){
        	cerr<<"Runtime Error"<<endl;
		break;
	     }
	}
	VC.getresult();
    }else{
	while(1){
	    if(FC.quit()==true&&stoptype==2){
		break;
	    }else{
		Meshtype1 M=FC.getmesh();
		valarray<double> arr=FC.generateRH(M);
		FC.setfirstin();
		FC.setfirstin();
	        M=FC.FMGcycle(arr,looptime);
       	        FC.setmesh(M);
	        double err=FC.generateerror();
		looptime++;
		if(looptime>maxloops&&stoptype==1){
		    break;
		}
	    }
	    if(looptime>=100){
		cerr<<"Runtime Error"<<endl;
		break;
	    }
	}
	FC.getresult();
    }
}
