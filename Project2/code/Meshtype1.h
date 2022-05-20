/**
 * @file   Meshtype1.h
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Wed Jun  2 23:48:31 2021
 * 
 * @brief  Class Meshtype1:the operators will be defined just like the lecture note.
 * 
 * 
 */
#ifndef _HSMATH_MESHTYPE1
#define _HSMATH_MESHTYPE1
#include"Mesh.h"
using namespace std;
class Meshtype1:public Mesh{
public:
    Meshtype1(int sz):Mesh(sz){};
    void tocoarse();
    void torefine();
    Meshtype1 operator=(Meshtype1 M1){
	size=M1.size;
	n=M1.n;
	nodes=M1.nodes;
	return *this;
    }
    friend Meshtype1 operator+(const Meshtype1& M1,const Meshtype1& M2){
	if(M1.size!=M2.size){
	    cerr<<"Length is different!"<<endl;
	}
	Meshtype1 ret(M1.size+1);
	for(int i=0;i<M1.size;i++){
	    ret.nodes[i]=M1.nodes[i]+M2.nodes[i];
	}
	return ret;
    }
};
/** 
 * Function tocoarse:change refine mesh to coarse mesh.
 * 
 */
void Meshtype1::tocoarse(){
    n=n/2;
    size=n-1;
    valarray<double> coarsenode(size);//Coarse Mesh
    for(int i=0;i<size;i++){
	coarsenode[i]=nodes[i*2]/4+nodes[i*2+1]/2+nodes[i*2+2]/4;
    }
    nodes=coarsenode;
}
/** 
 * Function torefine:change coarse mesh to refine mesh
 * 
 */
void Meshtype1::torefine(){
    n=n*2;
    size=n-1;
    valarray<double> refinenode(size);//Refine Mesh
    refinenode[0]=(nodes[0])/2;
    for(int i=1;i<size;i=i+2){
	refinenode[i]=nodes[i/2];
    }
    for(int i=2;i<size-1;i=i+2){
	refinenode[i]=(nodes[i/2]+nodes[i/2-1])/2;
    }
    refinenode[size-1]=(nodes[(size-1)/2-1])/2;
    nodes=refinenode;
}
#else
//NOTHING
#endif
