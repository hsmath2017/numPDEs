/**
 * @file   testMesh.cpp
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Thu Jun  3 01:20:57 2021
 * 
 * @brief  Test the generalization of Mesh, and two refining operators.
 * 
 * 
 */
#include "Meshtype1.h"
using namespace std;
int main(){
    Meshtype1 Mesh(4);
    Mesh.printMesh();
    Mesh.torefine();
    Mesh.printMesh();
    Mesh.tocoarse();
    Mesh.printMesh();
}
