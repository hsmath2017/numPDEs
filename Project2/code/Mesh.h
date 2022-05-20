/*-*- C++ -*-*/
/**
 * @file   Mesh.h
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Wed Jun  2 23:09:59 2021
 * 
 * @brief  Class Mesh:Store the Mesh in the Multigrid model.
 * Virtual base class, the operators I will be defined later.
 * 
 */
#ifndef _HSMATH_MESH
#define _HSMATH_MESH
#include<iostream>
#include<cmath>
#include<valarray>
using namespace std;
class Mesh{
protected:
    int n;//The number of intervals
    valarray<double> nodes;//The information of vector v
    int size;//size of the vector 'nodes'
public:
    Mesh(int sz);//Constructor,create the mesh of BVP with given size
    void setsize(int sz){
	size=sz;
	n=size+1;
	nodes.resize(size);
    }
    void setzero(){
	for(int i=0;i<size;i++){
	    nodes[i]=0;
	}
    }
    int getsize(){
	return size;
    }
    valarray<double> getinfo(){
	return nodes;
    }
    void setinfo(valarray<double> vec){
	nodes=vec;
	size=vec.size();
	n=size+1;
    }
    virtual void tocoarse()=0;//Refine mesh change into coarse mesh
    virtual void torefine()=0;//Coarse mesh change into refine mesh
    void printMesh(){
	cout<<nodes[0];
	for(int i=1;i<size;i++){
	    cout<<" "<<nodes[i];
	}
	cout<<endl;
    }
};
Mesh::Mesh(int sz){
    n=sz;
    size=sz-1;
    nodes.resize(size);//x_{1} to x_{n-1}
    for(int i=0;i<size;i++){//Initial guess:zeroes!
	nodes[i]=0;
    }
}

#else
//NOTHING
#endif
