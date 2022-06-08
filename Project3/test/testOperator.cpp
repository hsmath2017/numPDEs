#include "../src/MultiGrid/Restrictor.h"
#include "../src/MultiGrid/Interpolation.h"
int main(){
    Box<2> bx2({0,0},{4,4});
    RectDomain<2> RD2(bx2,{1.0/4,1.0/4},NodeCentered,0);
    Box<2> bx1({0,0},{2,2});
    RectDomain<2> RD1(bx1,{1.0/2,1.0/2},NodeCentered,0);
    FullWeighting<2> FW(RD1,RD2);
    Tensor<Real,2> T2(bx2);
    Tensor<Real,2> T1(bx1);
    loop_box_2(bx2,i,j){
        T2(i,j)=i+j;
    }
    FW.apply(T2,T1);
    std::cout<<"T1 = "<<T1<<std::endl;
    std::cout<<"T2 = "<<T2<<std::endl;
    LinearInterpolator<2> LI(RD1,RD2);
    LI.apply(T1,T2);
    std::cout<<"T1 = "<<T1<<std::endl;
    std::cout<<"T2 = "<<T2<<std::endl;    
}